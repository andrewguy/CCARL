'''General utility functions'''

import logging
import sys
import numpy as np
import pandas as pd
import click
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


class OptionEatAll(click.Option):
    '''Click Option class that allows for multiple arguments, similar to *args.

    In some versions of Click, needs to be used with type=tuple:

        @click.option("--arg", cls=OptionEatAll, type=tuple)

    From https://stackoverflow.com/questions/48391777/nargs-equivalent-for-options-in-click
    '''
    def __init__(self, *args, **kwargs):
        self.save_other_options = kwargs.pop('save_other_options', True)
        nargs = kwargs.pop('nargs', -1)
        assert nargs == -1, 'nargs, if set, must be -1 not {}'.format(nargs)
        super(OptionEatAll, self).__init__(*args, **kwargs)
        self._previous_parser_process = None
        self._eat_all_parser = None

    def add_to_parser(self, parser, ctx):

        def parser_process(value, state):
            # method to hook to the parser.process
            done = False
            value = [value]
            if self.save_other_options:
                # grab everything up to the next option
                while state.rargs and not done:
                    for prefix in self._eat_all_parser.prefixes:
                        if state.rargs[0].startswith(prefix):
                            done = True
                    if not done:
                        value.append(state.rargs.pop(0))
            else:
                # grab everything remaining
                value += state.rargs
                state.rargs[:] = []
            value = tuple(value)

            # call the actual process
            self._previous_parser_process(value, state)

        retval = super(OptionEatAll, self).add_to_parser(parser, ctx)
        for name in self.opts:
            our_parser = parser._long_opt.get(name) or parser._short_opt.get(name)
            if our_parser:
                self._eat_all_parser = our_parser
                self._previous_parser_process = our_parser.process
                our_parser.process = parser_process
                break
        return retval
