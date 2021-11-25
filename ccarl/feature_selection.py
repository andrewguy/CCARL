
from tempfile import NamedTemporaryFile
import subprocess
import os


def run_fast_mrmr(feature_df, mrmr_bin, mrmr_reader, num_features=10):
    '''Run fast-mRMR algorithm on a feature dataframe.

    Feature dataframe should contain class label in column 0, with features in
    every other column.

    Args:
        feature_df (DataFrame): A DataFrame containing class labels and features
        mrmr_bin_dir (str): Path to the mRMR binary.
        mrmr_reader (str): Path to the mRMR data reader binary.
        num_features (int, optional): The number of features to return from mRMR
            algorithm. Defaults to 10 features.

    Return:
        DataFrame: A filtered subset of the original DataFrame containing
            selected features.
    '''
    suffix = ".csv"
    with NamedTemporaryFile(mode='w', delete=False, suffix=suffix) as f:
        feature_df.to_csv(path_or_buf=f, index=False)
        name = f.name

    prefix = name[:-len(suffix)]
    mrmr_data = prefix + '.mrmr'

    subprocess.run(f'{mrmr_reader} -f {name} -o {mrmr_data}', shell=True)
    output = subprocess.check_output(f'{mrmr_bin} -f {mrmr_data} -a {num_features}', shell=True)
    os.remove(name)
    os.remove(mrmr_data)
    features = [int(x) for x in output.decode().strip().split(',')]
    return list(feature_df.iloc[:, features].columns.values)
