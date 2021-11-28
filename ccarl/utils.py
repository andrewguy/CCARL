'''General utility functions'''

import logging
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold


def _setup_logging(verbose=False):
    logging.basicConfig(
        format='%(levelname)s:%(message)s',
        stream=sys.stdout,
        level=logging.DEBUG if verbose else logging.WARNING
    )


def generate_cv_folds(df, n_splits=5, random_seed=None):
    if random_seed:
        np.random.seed(random_seed)
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True)
    folds = []
    binding = df.Binding.to_numpy(dtype=int)
    for train, test in skf.split(df.Structure, binding):
        glycan_list_train = [df.Structure[x] for x in train]
        glycan_list_test = [df.Structure[x] for x in test]
        df_train = pd.DataFrame({'Structure': glycan_list_train, 'Binding': binding[train]})
        df_test = pd.DataFrame({'Structure': glycan_list_test, 'Binding': binding[test]})
        folds.append((df_train, df_test))
    return folds
