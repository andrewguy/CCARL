import logging
import pickle

import click
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import auc, roc_curve

from ccarl.ccarl import CCARLClassifier, _calculate_binders, _log_rfu_values
from ccarl.glycan_graph_methods import digraph_to_glycan_string
from ccarl.plotting.metrics import plot_kfold_test_training_roc, plot_test_training_roc
from ccarl.plotting.microarray import plot_log_rfu_histogram
from ccarl.plotting.utils import remove_top_right_borders
from ccarl.utils import _setup_logging, generate_cv_folds
from ccarl.validate_cfg_structures import (ArrayMismatchError,
                                           validate_cfg_structures)


@click.command()
@click.argument('input', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output', type=click.Path(exists=False, dir_okay=False), required=True)
@click.option('--cfg_version', default='',
              help='CFG array version number. If present, will check that glycan strings provided roughly match array version.')
@click.option('--levenshtein_threshold', type=float, default=10,
              help='Total Levenshtein distance threshold for calling a match.')
@click.option('-v', '--verbose', is_flag=True,
              help='Verbose output. Will print number of mismatches to stdout.')
def validate_structures(input, output, cfg_version, levenshtein_threshold, verbose):
    '''Validate CFG glycan structures. Output corrected glycan strings.

    Takes a CSV file as input, and writes a similar CSV output with corrected glycan strings.
    Input CSV file can contain additional columns such as RFU values etc.
    '''
    logger = logging.getLogger(__name__)
    _setup_logging(verbose)
    if verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
    csv_data = pd.read_csv(input)
    try:
        out_df = validate_cfg_structures(csv_data, cfg_version, levenshtein_threshold, log=logger)
    except ArrayMismatchError:
        return
    out_df.to_csv(output, index=False)
    return


@click.command()
@click.argument('input', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output', type=click.Path(exists=False, dir_okay=False), required=True)
@click.option('--zscore_low', default=1.5,
              help='Z-score threshold for negative binders')
@click.option('--zscore_high', default=3.5,
              help='Z-score threshold for positive binders')
@click.option('--histogram', type=click.Path(exists=False, dir_okay=False), required=False,
              help='Filename for histogram of binding levels.')
def identify_binders(input, output, zscore_low, zscore_high, histogram):
    '''Identify binding and non-binding glycans from glycan microarray data.

    Takes a CSV file as input. CSV input file must contain a column labelled 'RFU'.
    Writes an output CSV file with an additional 'Binding' column with binary 0 or 1 indicators of binding.
    Intermediate/indeterminate binding glycans will not be written to the output CSV file.
    Set zscore_low = zscore_high if you want to ignore intermediate
    '''
    logger = logging.getLogger(__name__)
    _setup_logging()
    csv_data = pd.read_csv(input)
    thresholds = (zscore_low, zscore_high)
    if zscore_low > zscore_high:
        logger.error("zscore_low can't be higher than zscore_high. Aborting.")
        return
    csv_data['log_rfu'] = _log_rfu_values(csv_data.RFU)
    csv_data['ternary_binding_class'] = _calculate_binders(csv_data.log_rfu, thresholds=thresholds)

    csv_data_subset = csv_data[csv_data.ternary_binding_class != 1].copy()
    csv_data_subset['Binding'] = csv_data_subset['ternary_binding_class'].astype(bool).astype(int)
    csv_data_subset.to_csv(output, index=False)
    if histogram:
        fig, ax = plt.subplots()
        remove_top_right_borders(ax)
        plot_log_rfu_histogram(csv_data.log_rfu, csv_data.ternary_binding_class, ax, legend=True)
        fig.savefig(histogram)
    return


@click.command()
@click.argument('input', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output_prefix', type=click.Path(exists=False, dir_okay=False), required=True)
@click.option('--support_positive', default=0.4,
              help='Support threshold for frequent subtree mining the positive binding set')
@click.option('--support_all', default=0.05,
              help='Support threshold for frequent subtree mining all glycans')
@click.option('--format', default='CFG',
              help='Glycan string format for input file. Currently only CFG is supported.')
@click.option('--cross_validation', is_flag=True,
              help='Run 5-fold cross-validation to assess model performance.')
@click.option('--plot_roc', is_flag=True,
              help='Plot ROC curves')
@click.option('--save_model', is_flag=True,
              help='Save model as a Pickle file')
@click.option('-v', '--verbose', is_flag=True,
              help='Verbose output.')
def identify_motifs(input, output_prefix, support_positive, support_all, format, cross_validation,
                    plot_roc, verbose, save_model):
    '''Extract glycan motifs from glycan microarray data.

    Requires a CSV file containing 'Structure' and 'Binding' columns.
    'Binding' column should contain a binary indicator (1=binding, 0=non-binder) for binding.
    Will output identified motifs to stdout, along with model performance metrics.
    '''
    csv_data = pd.read_csv(input)

    cf = CCARLClassifier(support_pos=support_positive, support_all=support_all, num_mrmr_features=10)
    cf.train(csv_data.Structure, csv_data.Binding.to_numpy(dtype=int), glycan_format=format, parse_linker=True)

    output_results_for_training(cf, csv_data, test=None, title="----Training results on full dataset:----\n")
    if cross_validation:
        folds = generate_cv_folds(csv_data)
        models = []
        for i, (train, test) in enumerate(folds):
            cf = CCARLClassifier(support_pos=support_positive, support_all=support_all, num_mrmr_features=10)
            cf.train(train.Structure, train.Binding.to_numpy(dtype=int), glycan_format=format, parse_linker=True)
            output_results_for_training(cf, train, test, title=f"----Training results for fold {i+1}----\n")
            models.append(cf)
    if plot_roc and cross_validation:
        fig, ax = plt.subplots(figsize=(4, 3))
        plot_kfold_test_training_roc(ax, models, *zip(*folds))
        fig.savefig(f"{output_prefix}_ROC_curves_CV.svg")
    if plot_roc:
        fig, ax = plt.subplots(figsize=(4, 3))
        plot_test_training_roc(ax, cf, csv_data, test_df=None)
        fig.savefig(f"{output_prefix}_ROC_curve_full_data.svg")
    if save_model:
        with open(f"{output_prefix}_model.pkl", mode='w') as f:
            pickle.dump(f, cf)


def output_results_for_training(cf, train, test=None, title=''):
    print(title)
    coefs = cf._model.coef_[0]
    fpr, tpr, _ = roc_curve(train.Binding, cf.predict_proba(train.Structure)[:, 1], drop_intermediate=False)
    auc_value = auc(fpr, tpr)
    sorted_coefs, sorted_features = zip(*sorted(zip(coefs, cf.subtree_features), key=lambda x: -x[0]))
    print("Features:\n")
    for i, feature in enumerate(sorted_features):
        feature_str = digraph_to_glycan_string(feature)
        print(f"Feature {i+1}: {feature_str}")
    print("\nModel coefficients:\n")
    for i, coef in enumerate(sorted_coefs):
        print(f"Feature {i+1} model coefficient: {coef: 0.3f}")
    print(f"\nAUC Value (Train): {auc_value: 0.4f}")
    if test is not None:
        fpr, tpr, _ = roc_curve(test.Binding, cf.predict_proba(test.Structure)[:, 1], drop_intermediate=False)
        auc_value = auc(fpr, tpr)
        print(f"AUC Value (Test): {auc_value: 0.4f}")
    print("\n")
    return
