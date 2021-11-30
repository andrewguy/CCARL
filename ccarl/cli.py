import logging
import pickle

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import auc, roc_curve

from ccarl import CCARLClassifier
from ccarl.ccarl import _calculate_binders, _log_rfu_values
from ccarl.glycan_graph_methods import digraph_to_glycan_string
from ccarl.plotting.features import render_features_pdf, render_glycan_list_pdf
from ccarl.plotting.metrics import plot_kfold_test_training_roc, plot_test_training_roc
from ccarl.plotting.microarray import plot_log_rfu_histogram
from ccarl.plotting.utils import remove_top_right_borders
from ccarl.utils import _setup_logging, generate_cv_folds, OptionEatAll
from ccarl.validate_cfg_structures import (ArrayMismatchError,
                                           validate_cfg_structures)


@click.group()
def cli():
    '''A collection of utilities for processing glycan microarray data.

    If you are starting with CFG data, it is recommended to first clean
    your data with `ccarl validate-structures`, and then identify binding
    glycans using median absolute deviation values using `ccarl identify-binders`.
    From here you can identify motifs (and build a predictive model) using
    `ccarl identify-motifs`.
    Finally, you can predict the binding of unknown glycans using `ccarl predict-binding`.
    '''
    pass


@cli.command()
@click.argument('input', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output', type=click.Path(exists=False, dir_okay=False), required=True)
@click.option('--cfg-version', default='',
              help='CFG array version number. Will check that glycan strings provided roughly match array version.')
@click.option('--levenshtein-threshold', type=float, default=10,
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


@cli.command()
@click.argument('input', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output', type=click.Path(exists=False, dir_okay=False), required=True)
@click.option('--zscore-low', default=1.5,
              help='Z-score threshold for negative binders (default: 1.5)')
@click.option('--zscore-high', default=3.5,
              help='Z-score threshold for positive binders (default: 3.5)')
@click.option('--histogram', type=click.Path(exists=False, dir_okay=False), required=False,
              help='Filename for histogram of binding levels. Format is determined by file extension.')
def identify_binders(input, output, zscore_low, zscore_high, histogram):
    '''Identify binding and non-binding glycans from glycan microarray data.

    Takes a CSV file as input. CSV input file must contain a column labelled 'RFU'.
    Writes an output CSV file with an additional 'Binding' column with binary 0 or 1 indicators of binding.

    Glycan binding/non-binding is calculated using Median Absolute Deviation.

    Glycans that fall between zscore_low and zscore-high will be classed as "intermediate" binders. These
    will not be written to the output CSV file, and will be ignored for subsequent analyses.
    If you don't want this behaviour, set zscore-high and zscore-low to the same value.
    '''
    logger = logging.getLogger(__name__)
    _setup_logging()
    csv_data = pd.read_csv(input)
    thresholds = (zscore_low, zscore_high)
    if zscore_low > zscore_high:
        logger.error("zscore-low can't be higher than zscore-high. Aborting.")
        return
    csv_data['log_rfu'] = _log_rfu_values(csv_data.RFU)
    csv_data['ternary_binding_class'] = _calculate_binders(csv_data.log_rfu, thresholds=thresholds)

    csv_data_subset = csv_data[csv_data.ternary_binding_class != 1].copy()
    csv_data_subset['Binding'] = csv_data_subset['ternary_binding_class'].astype(bool).astype(int)
    csv_data_subset.drop(columns=['ternary_binding_class'], inplace=True)
    csv_data_subset.to_csv(output, index=False)
    if histogram:
        fig, ax = plt.subplots()
        remove_top_right_borders(ax)
        plot_log_rfu_histogram(csv_data.log_rfu, csv_data.ternary_binding_class, ax, legend=True)
        fig.savefig(histogram)
    return


@cli.command()
@click.argument('input', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output_prefix', type=click.Path(exists=False, dir_okay=False), required=True)
@click.option('--support-positive', default=0.4,
              help='Support threshold for frequent subtree mining the positive binding set (default=0.4)')
@click.option('--support-all', default=0.05,
              help='Support threshold for frequent subtree mining all glycans (default=0.05)')
@click.option('--format', default='CFG',
              help='Glycan string format for input file. Currently only CFG is supported. (default=CFG)')
@click.option('--cross-validation', is_flag=True,
              help='Run 5-fold cross-validation to assess model performance.')
@click.option('--plot-roc', is_flag=True,
              help='Plot ROC curves')
@click.option('--save-model', is_flag=True,
              help='Save model as a Pickle file')
@click.option('--render-motifs', is_flag=True,
              help='Save PDF of motifs as SNFG-style drawings')
@click.option('--render-motifs-cv', is_flag=True,
              help='Save PDF of motifs as SNFG-style drawings for all cross-validation folds')
def identify_motifs(input, output_prefix, support_positive, support_all, format, cross_validation,
                    plot_roc, save_model, render_motifs, render_motifs_cv):
    '''Extract glycan motifs from glycan microarray data.

    Requires a CSV file containing 'Structure' and 'Binding' columns.
    'Binding' column should contain a binary indicator (1=binding, 0=non-binder) for binding.
    Will output identified motifs to stdout, along with model performance metrics.
    '''
    logger = logging.getLogger(__name__)
    _setup_logging()
    if render_motifs_cv and not cross_validation:
        logger.error("--render-motifs-cv was set, but not --cross-validation. Did you mean to set --cross-validation? "
                     "Run aborted.")
        return
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
            if render_motifs_cv:
                render_features_pdf(cf, f"{output_prefix}.cv_fold_{i+1}_features.pdf")
    if plot_roc and cross_validation:
        fig, ax = plt.subplots(figsize=(4, 3))
        plot_kfold_test_training_roc(ax, models, *zip(*folds))
        fig.savefig(f"{output_prefix}.ROC_curves_CV.svg", bbox_inches="tight")
    if plot_roc:
        fig, ax = plt.subplots(figsize=(4, 3))
        plot_test_training_roc(ax, cf, csv_data, test_df=None)
        fig.savefig(f"{output_prefix}.ROC_curve_full_data.svg", bbox_inches="tight")
    if save_model:
        with open(f"{output_prefix}.model.pkl", mode='wb') as f:
            pickle.dump(cf, f)
    if render_motifs:
        render_features_pdf(cf, f"{output_prefix}.features.pdf")


def output_results_for_training(cf, train=None, test=None, title=''):
    print(title)
    coefs = cf._model.coef_[0]
    sorted_coefs, sorted_features = zip(*sorted(zip(coefs, cf.subtree_features), key=lambda x: -x[0]))
    print("Features:\n")
    for i, feature in enumerate(sorted_features):
        feature_str = digraph_to_glycan_string(feature)
        print(f"Feature {i+1}: {feature_str}")
    print("\nModel coefficients:\n")
    for i, coef in enumerate(sorted_coefs):
        print(f"Feature {i+1} model coefficient: {coef: 0.3f}")
    if train is not None:
        fpr, tpr, _ = roc_curve(train.Binding, cf.predict_proba(train.Structure)[:, 1], drop_intermediate=False)
        auc_value = auc(fpr, tpr)
        print(f"\nAUC Value (Train): {auc_value: 0.4f}")
    if test is not None:
        fpr, tpr, _ = roc_curve(test.Binding, cf.predict_proba(test.Structure)[:, 1], drop_intermediate=False)
        auc_value = auc(fpr, tpr)
        print(f"AUC Value (Test): {auc_value: 0.4f}")
    print("\n")
    return


@cli.command()
@click.argument('input', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('model', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output', type=click.Path(exists=False, dir_okay=False), required=True)
@click.option('--format', default='CFG',
              help='Glycan string format for input file. Currently only CFG is supported (default=CFG).')
@click.option('-v', '--verbose', is_flag=True,
              help='Verbose output. Will print model details to stdout.')
def predict_binding(input, model, output, format, verbose):
    '''Predict binding of unknown glycans from model trained on glycan microarray data.

    Requires a previously trained model file (a Python pickle object).
    Use ccarl identify-motifs with the --save-model flag to generate a trained model.

    Requires a CSV file containing a 'Structure' column.
    Will output binding probability to a csv file.
    '''
    csv_data = pd.read_csv(input)
    with open(model, 'rb') as f:
        cf = pickle.load(f)
    preds = cf.predict_proba(csv_data.Structure, glycan_format=format)
    csv_data['Binding_Probability'] = preds[:, 1]
    add_features_to_df(csv_data, cf, glycan_format=format)
    csv_data.to_csv(output, index=False, float_format='%1.3g')
    if verbose:
        output_results_for_training(cf, title="----Model Details----\n")
        print(f"Binding predictions saved to {output}")
    return


def add_features_to_df(csv_data, model, glycan_format='CFG', parse_linker=True):
    '''Add feature columns to DataFrame'''
    features = model._generate_features(csv_data.Structure, glycan_format=glycan_format, parse_linker=parse_linker)
    features = np.array(features)
    for i, idx in enumerate(np.argsort(-model._model.coef_[0])):
        csv_data[f'Feature_{i+1}'] = features[:, idx]
    return


@cli.command()
@click.argument('input', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output', type=click.Path(exists=False, dir_okay=False), required=True)
@click.option('--format', default='CFG',
              help='Glycan string format for input file. Currently only CFG is supported (default=CFG).')
def render_glycans(input, output, format):
    '''Render a list of glycans as SNFG-like symbols

    Does not currently support rendering of glycan strings which contain restricted linkage information.

    Requires a CSV file containing a 'Structure' column.
    Will save glycans as a PDF file, with one glycan per page. Programs such as Inkscape can be
    used to load individual pages from this PDF file and save in other formats.
    '''
    csv_data = pd.read_csv(input)
    render_glycan_list_pdf(csv_data.Structure, output, format=format)
    return


@cli.command()
@click.argument('input', type=click.Path(exists=True, dir_okay=False), required=True)
@click.option("--models", cls=OptionEatAll, help='Previously trained models (accepts multiple arguments)', type=tuple)
@click.option('--format', default='CFG',
              help='Glycan string format for input file. Currently only CFG is supported (default=CFG).')
def binding_overlap(input, models, format):
    '''Cross-tabulate predicted binding from multiple models and some other binary
    separation of glycans.

    For example, you may have a set of glycans that are present on Cell Type 1, and a set of glycans that are
    present on Cell Type 2, and you wish to find a lectin which will allow you to adequately distinguish
    the two cell types. Using a selection of previously trained models, you can use this tool to
    assess which of the trained models adequately separates the pre-defined classes (e.g. Cell Type 1
    and Cell Type 2).

    The input CSV file should contain a `Structure` column of glycan strings and a
    `Class` column with the different glycan classes. A glycan can be present in multiple classes, and
    in that case it should be present on multiple rows.

    Use `ccarl identify-motifs` with the `--save-model` flag to generate trained models.
    '''
    csv_data = pd.read_csv(input)
    model_results = []
    for model in models:
        print(f"----Predicted binding cross-tab for model {model}----\n")
        with open(model, 'rb') as f:
            cf = pickle.load(f)
        preds = cf.predict_proba(csv_data.Structure, glycan_format=format)[:, 1]
        model_results.append(preds > 0.5)
        crosstab = pd.crosstab(csv_data['Class'], preds > 0.5, margins=False, colnames=['Predicted Binder'])
        print(crosstab)
        print()
    if len(models) > 1:
        crosstab = pd.crosstab(csv_data['Class'], model_results, margins=False, colnames=models)
        print("----Predicted binding cross-tab for all models----\n")
        print(crosstab)
    return
