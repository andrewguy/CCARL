import numpy as np


def mad_based_outlier(points, thresh=3.5, two_sided=False):
    '''Identifies outliers using Median Absolute Deviation.

    By default, returns outliers who are greater than the main distribution.

    Args:
        points (numpy.array): Data points.
        thresh (float): Modified Z-score threshold.
        two_sided (bool, default=False): If True, returns outliers both greater
            than and less than the main distribution.
    Returns:
        numpy.array: A boolean array of outliers.
    '''
    if len(points.shape) == 1:
        points = np.array(points)[:, None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
    signed_diff = np.sum((points - median), axis=-1)
    modified_z_score = 0.6745 * signed_diff / med_abs_deviation
    return modified_z_score > thresh
