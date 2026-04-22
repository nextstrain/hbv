import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--metadata-in", required=True)
parser.add_argument("--mapping", required=True)
parser.add_argument("--metadata-out", required=True)
args = parser.parse_args()

metadata = pd.read_csv(args.metadata_in, sep="\t", low_memory=False)

mapping = pd.read_csv(
    args.mapping,
    sep=",",
    comment="#",
    dtype=str,
    keep_default_na=False,
)

mapping["old"] = mapping["old"].str.strip()
mapping["new"] = mapping["new"].str.strip()

mapping_dict = dict(zip(mapping["old"], mapping["new"]))

metadata["subgenotype_genbank"] = metadata["subgenotype_genbank"].replace(mapping_dict)

metadata.loc[metadata["subgenotype_genbank"] == "REMOVE", "subgenotype_genbank"] = ""

metadata.to_csv(args.metadata_out, sep="\t", index=False)
