'''General utility functions'''

import logging
import sys


def _setup_logging(verbose=False):
    logging.basicConfig(
        format='%(levelname)s:%(message)s',
        stream=sys.stdout,
        level=logging.DEBUG if verbose else logging.WARNING
    )
