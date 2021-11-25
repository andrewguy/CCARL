"""A module for correcting errors in CFG glycan strings.

Most data files from the CFG do not strictly conform to the modified-CNFG format that
glycan data should be provided as. Problems include the end of long glycan strings
being cut off due to field size limits in Excel, missing symbols such as opening
or closing brackets, use of non-breaking space characters instead of standard spaces,
occasional use of greek alpha and beta instead of 'a' and 'b', etc.

If using data from the CFG, this module should first be used to fix any problems in the
glycan strings that would then cause issues with the CCARL CFG parsing module (which
cannot fix syntax errors in glycan strings).

This module uses a string-matching approach using Levenshtein distances to compare
glycan strings between provided data and corrected data for each CFG array version.
If the CFG array version is provided, this module will verify that this is the closest
matching array version, as well as verifying that the provided data matches within some
tolerance limit (set by the `levenshtein_threshold` parameter). If no CFG array version
is provided, this module will identify the closest matching array version (within a
tolerance limit).In both cases this module returns the idealised glycan strings for
that array version.
"""

import logging

from .glycan_parsers.cfg_array_versions import get_likely_cfg_array_version

logger = logging.getLogger(__name__)


def validate_cfg_structures(input_df, cfg_version='', levenshtein_threshold=2.0, log=logger):
    """Find closest matching set of CFG structures.

    Args:
        input_df (pandas.DataFrame): CFG Data. Glycans should be in a column labelled 'Structure'.
        cfg_version (str, optional): CFG version to match. Will match any version if not supplied. Defaults to ''.
        levenshtein_threshold (float, optional): Total Levenshtein distance for matching glycan names. Defaults to 2.
        log (Log, optional): Allows setting verbosity of output. Defaults to Log().

    Raises:
        ValueError: Raised if no matching CFG version is identified.

    Returns:
        pandas.DataFrame: A copy of the input dataframe with Structure column replaced with correct glycan strings.
    """
    structures = list(input_df.Structure)
    threshold = levenshtein_threshold
    glycans, array_ver, mismatches, sum_lev = get_likely_cfg_array_version(structures, distance_threshold=threshold)
    if cfg_version and (cfg_version != array_ver):
        msg = f"Supplied CFG version ({cfg_version}) does not match likely CFG version ({array_ver})."
        log.error(msg)
        raise ArrayMismatchError(msg)
    log.info(f"Matched array version is {array_ver} with {mismatches} mismatches and total levenshtein distance of {sum_lev}.")
    csv_data = input_df.copy()
    csv_data.Structure = glycans
    csv_data.attrs['cfg_version'] = array_ver
    return csv_data


class ArrayMismatchError(ValueError):
    '''Raised when there is a mismatch in expected vs closest CFG array versions'''
    pass
