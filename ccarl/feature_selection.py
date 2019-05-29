
from tempfile import NamedTemporaryFile, _get_candidate_names, gettempdir
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
    tempdir = gettempdir()
    mrmr_data = os.path.join(tempdir, prefix + '.mrmr')

    os.system('cd {tempdir} && {mrmr_reader} {name}'.format(name=name, mrmr_reader=mrmr_reader, 
                                                            tempdir=tempdir))
    os.rename(os.path.join(tempdir, "data.mrmr"), mrmr_data)
    output = subprocess.check_output('{mrmr_bin} -f {mrmr_data} -a {n}'.format(
                 n=num_features + 1, mrmr_bin=mrmr_bin, mrmr_data=mrmr_data), shell=True)
    os.remove(name)
    os.remove(mrmr_data)
    features = [int(x) for x in output.decode().strip().split(',')]
    return list(feature_df.iloc[:, features].columns.values)