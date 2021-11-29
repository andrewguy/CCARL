import numpy as np


def plot_log_rfu_histogram(x, binding_classes, ax, title='', legend=False):
    '''Plot a histogram of counts and show positive and intermediate binders.

    Args:
        x: log RFU values
        ax (AxesSubplot): Matplotlib figure axis.
        title (str, optional): Title for figure.
        thresholds (tuple, optional): Modified z-score thresholds for intermediate and
            positive binders.

    Returns:
        AxesSubplot: Matplotlib figure axis, to allow method chaining.
    '''
    binwidth = 0.05
    data_neg = x[binding_classes == 0]
    data_pos = x[binding_classes == 2]
    data_int = x[binding_classes == 1]
    ax.hist(data_neg, bins=np.arange(min(data_neg), max(data_neg) + binwidth, binwidth))
    # Sometimes there aren't any intermediates or positives...
    legend = ['Negative', 'Positive']
    if len(data_int) != 0:
        bins = np.arange(min(data_int), max(data_int) + binwidth, binwidth)
        ax.hist(data_int, bins=bins, color='orange')
        legend = ['Negative', 'Intermediate', 'Positive']
    if len(data_pos) != 0:
        ax.hist(data_pos, bins=np.arange(min(data_pos), max(data_pos) + binwidth, binwidth),
                color='red')
    ax.set_title(title)
    ax.set_xlabel('log(RFU)')
    ax.set_ylabel('Counts')
    if legend:
        ax.legend(legend)
    return ax
