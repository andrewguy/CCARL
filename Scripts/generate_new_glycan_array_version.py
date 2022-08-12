import pandas as pd

from ccarl.glycan_parsers.cfg_array_versions import clean_glycan_string
from ccarl.glycan_parsers.cfg_parser import CFGGlycanParser

"""
This is an example script for adding new glycan array versions to the CCARL tool.
"""

if __name__ == "__main__":
    # Add your path to cfg data in csv format here. Note: csv file must include 'Structure' and 'Chart Number' columns
    input_path = f"<path-to-cfg-data>.csv"
    # Add CFG array version number here
    new_array_version = "x.x"
    df = pd.read_csv(input_path)
    # This path can be left as is
    output_path = f"Data/CFG_Array_Versions/Array_{new_array_version}.csv"
    parser = CFGGlycanParser()
    failed, succeeded = [], []
    for structure in df["Structure"]:
        try:
            structure = clean_glycan_string(structure)
            _ = parser.string_to_graph(structure)
            succeeded.append(structure)
        except Exception as e:
            failed.append(structure)

    print(f"Failed to parse {len(failed)} structures. {len(succeeded)} parsed successfully ")
    print(f"Failed structures that will need to be fixed manually:")
    for structure in failed:
        print(structure)

    if not failed:  # Only write new glycan array file if all structures were successfully parsed.
        df['Structure'] = succeeded
        df.to_csv(output_path, columns=["Chart Number", "Structure"], index=False, header=False)
        print(f"Created new file: {output_path}")
