from pytest import fixture
import pytest
import pandas as pd
from ccarl.validate_cfg_structures import validate_cfg_structures

TEST_CSV_DATA = './tests/test_data/MAL_I_0.1ug_13881_v5.0_DATA.csv'


@fixture
def csv_df():
    return pd.read_csv(TEST_CSV_DATA)


def test_validate_without_version(csv_df):
    out_df = validate_cfg_structures(csv_df)
    assert len(out_df.Structure) == 611
    struct57 = out_df[out_df['Chart Number'] == 57].Structure.iloc[0]
    assert struct57 == 'Neu5Aca2-6Galb1-4GlcNAcb1-2Mana1-6(Neu5Aca2-6Galb1-4GlcNAcb1-2Mana1-3)Manb1-4GlcNAcb1-4GlcNAcb-Sp21'


def test_validate_with_wrong_version(csv_df):
    for version in ["2.0", "2.1", "3.0", "3.1", "3.2", "4.0", "4.1", "4.2", "5.1", "5.2"]:
        with pytest.raises(ValueError):
            _ = validate_cfg_structures(csv_df, cfg_version=version)


def test_validate_with_right_version(csv_df):
    version = "5.0"
    out_df = validate_cfg_structures(csv_df, cfg_version=version)
    assert len(out_df.Structure) == 611
    struct57 = out_df[out_df['Chart Number'] == 57].Structure.iloc[0]
    assert struct57 == 'Neu5Aca2-6Galb1-4GlcNAcb1-2Mana1-6(Neu5Aca2-6Galb1-4GlcNAcb1-2Mana1-3)Manb1-4GlcNAcb1-4GlcNAcb-Sp21'


def test_validate_and_access_version(csv_df):
    version = "5.0"
    out_df = validate_cfg_structures(csv_df, cfg_version=version)
    assert out_df.attrs['cfg_version'] == version
