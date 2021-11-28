from sklearn.metrics import roc_curve, auc
from matplotlib.lines import Line2D
import numpy as np


def plot_test_training_roc(ax, model, train_df, test_df=None):
    '''Plot test and training ROC curves.

    Args:
        ax [matplotlib.axes.Axes]: Matplotlib axis.
        model [CCarlClassifier]: CCARLClassifier object
        test_df: DataFrame containing test set
        train_df: DataFrame containing training set
    Returns:
        matplotlib.axes.Axes
    '''
    ax.set_title('')
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    clf = model
    glycan_train = train_df.Structure
    y_train = train_df.Binding
    if test_df is not None:
        glycan_test = test_df.Structure
        y_test = test_df.Binding
        fpr, tpr, _ = roc_curve(y_test, clf.predict_proba(glycan_test)[:, 1], drop_intermediate=False)
        ax.plot(fpr, tpr, color=(0, 0.45, 0.70), label=f'Test, AUC:{auc(fpr, tpr): 2.2f}')
    fpr_train, tpr_train, _ = roc_curve(y_train, clf.predict_proba(glycan_train)[:, 1], drop_intermediate=False)
    ax.plot(fpr_train, tpr_train, color=(0.8, 0.4, 0), label=f'Training, AUC:{auc(fpr_train, tpr_train): 2.2f}')
    ax.plot([0, 1], [0, 1], linestyle='--', color='grey', linewidth=0.8, dashes=(5, 10))
    ax.legend()
    return ax


def plot_kfold_test_training_roc(ax, models, train_dfs, test_dfs):
    '''Plot multiple ROC curves.

    Args:
        ax: Matplotlib Axes
        models: List of CCARLClassifier objects
        test_dfs: List of DataFrames containing test sets
        train_dfs: List of DataFrames containing training sets
    Returns:
        matplotlib.axes.Axes
    '''
    ax.set_title('')
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    aucs = []
    aucs_tr = []
    for model, train_df, test_df in zip(models, train_dfs, test_dfs):
        glycan_test = test_df.Structure
        y_test = test_df.Binding
        glycan_train = train_df.Structure
        y_train = train_df.Binding
        clf = model
        fpr, tpr, _ = roc_curve(y_test, clf.predict_proba(glycan_test)[:, 1], drop_intermediate=False)
        fpr_train, tpr_train, _ = roc_curve(y_train, clf.predict_proba(glycan_train)[:, 1], drop_intermediate=False)
        auc_ = auc(fpr, tpr)
        auc_tr = auc(fpr_train, tpr_train)
        aucs.append(auc_)
        aucs_tr.append(auc_tr)
        ax.plot(fpr, tpr, color=(0, 0.45, 0.70), alpha=0.5, label=f'Test, AUC:{auc_: 2.2f}')
        ax.plot(fpr_train, tpr_train, color=(0.8, 0.4, 0), alpha=0.5, label=f'Training, AUC:{auc_tr: 2.2f}')

    ax.plot([0, 1], [0, 1], linestyle='--', color='grey', linewidth=0.8, dashes=(5, 10))
    custom_lines = [Line2D([0], [0], color=(0, 0.45, 0.70), alpha=0.5),
                    Line2D([0], [0], color=(0.8, 0.4, 0), alpha=0.5)]
    ax.legend(custom_lines, [f'Test, AUC: {np.mean(aucs):2.2f} ({np.std(aucs):2.2f})',
                             f'Training, AUC:{np.mean(aucs_tr):2.2f} ({np.std(aucs_tr):2.2f})'])
    return ax
