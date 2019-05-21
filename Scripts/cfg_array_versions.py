'''A set of data for each CFG array version.'''

import os
import sys
import pandas as pd
import numpy as np
from itertools import zip_longest
from pyxdameraulevenshtein import damerau_levenshtein_distance, normalized_damerau_levenshtein_distance

cfg_array_versions = {}
DATA_DIR = '../Data/CFG_Array_Versions/'
for filename in os.listdir(DATA_DIR):
    if filename.endswith('.csv'):
        version = filename.replace('.csv', '').replace('Array_', '')
        cfg_array_versions[version] = pd.read_csv(os.path.join(DATA_DIR, filename), header=None)


def get_likely_cfg_array_version(glycan_list, mismatch_threshold=5):
    '''Get the most likely CFG glycan array given a list of glycans.

    Uses a scaled Levenshtein distance to compute similarity between glycan strings,
    and returns the array with the minimum sum of scaled levenshtein distances
    for each pair of glycans in the glycan list and corresponding reference array.

    We need to do this because sometimes the array version is not provided, and 
    there are slight spelling errors in the provided glycan names. It is easier to
    match to a reference list of glycans for a particular array version, with all 
    errors corrected.

    Args:
        glycan_list (list): A list of glycan strings ordered by index.
    Returns:
        CFG glycan list (list), most likely array version (string), number of mismatches (int)
    '''
    for i, glycan in enumerate(glycan_list):
        # Handle odd characters in some excel files. Nonbreaking spaces, greek letters etc.
        glycan_list[i] = glycan.replace('–', '-').replace('α', 'a') \
                            .replace('β', 'b').replace('[', '(') \
                            .replace(']', ')').replace(' ', '').replace(u"\u00A0", '')
    likely_array = None
    likely_array_mismatches = None
    scaled_levenshtein_total = 0
    for key, value in cfg_array_versions.items():
        # Take into account glycans which are almost the same.
        array_version = [x.replace(' ', '') for x in value[1]]
        scaled_levenshtein_sum = np.sum([normalized_damerau_levenshtein_distance(x, y) for x, y in zip_longest(glycan_list, array_version, fillvalue='')])
        non_matches = len([x for x in zip(glycan_list, array_version) if x[0] != x[1]])
        #if not likely_array or non_matches < likely_array_mismatches:
        if not likely_array or scaled_levenshtein_sum < scaled_levenshtein_total:
            likely_array = key
            likely_array_mismatches = non_matches
            scaled_levenshtein_total = scaled_levenshtein_sum
    if likely_array_mismatches > mismatch_threshold:
        raise ValueError("Glycan list does not match to known array versions.")
    return list(cfg_array_versions[likely_array][1]), likely_array, likely_array_mismatches, scaled_levenshtein_total
