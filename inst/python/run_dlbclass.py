"""Classify DLBClass GSM files with DLBclass-tool, writing each result where asked.

    python run_dlbclass.py --tool ~/git/DLBclass-tool gsm1.tsv out1.tsv [gsm2.tsv out2.tsv ...]
    python run_dlbclass.py --tool ~/git/DLBclass-tool --versions

Bundled with GAMBLR.predict and run by its run_dlbclass(); DLBclass-tool
(https://github.com/getzlab/DLBclass-tool) is installed separately.

Each GSM is in the DLBClass layout (a feature-name column, then one column
per sample), e.g. from GAMBLR.predict::dlbclass_gsm(..., output_format = "gsm").
Uses DLBclass-tool's own functions (construct_reduced_winning_version and
classify_samples_winning_model), so nothing is written inside the tool's
directory. The model paths in DLBclass-tool are relative, so it runs from
the tool's directory; GSM and output paths are resolved first.
"""
import argparse
import os
import sys
import warnings

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--tool", default=os.path.expanduser("~/git/DLBclass-tool"), help="DLBclass-tool checkout")
parser.add_argument("--versions", action="store_true", help="print Python, pandas and torch versions and exit")
parser.add_argument("pairs", nargs="*", help="GSM file and output file, repeated")
args = parser.parse_args()
if args.versions:
    import platform
    import pandas
    import torch
    print(f"python={platform.python_version()}\tpandas={pandas.__version__}\ttorch={torch.__version__}")
    sys.exit(0)
if not args.pairs or len(args.pairs) % 2:
    sys.exit("give GSM and output files in pairs")
pairs = [(os.path.abspath(os.path.expanduser(args.pairs[i])), os.path.abspath(os.path.expanduser(args.pairs[i + 1])))
         for i in range(0, len(args.pairs), 2)]

tool = os.path.abspath(os.path.expanduser(args.tool))
os.chdir(tool)
sys.path.insert(0, tool)
import pandas as pd  # noqa: E402
import src.classify_generic as cg  # noqa: E402
import src.format_data as fd  # noqa: E402

clusters = {1.0: "C1", 2.0: "C2", 3.0: "C3", 4.0: "C4", 5.0: "C5"}
for gsm_file, out_file in pairs:
    gsm = pd.read_csv(gsm_file, sep="\t", index_col=0)
    classified = cg.classify_samples_winning_model(fd.construct_reduced_winning_version(gsm))
    classified.index.name = "sample"
    classified["PredictedCluster"] = classified["PredictedCluster"].map(clusters)
    classified.to_csv(out_file, sep="\t", index=True)
    print(f"{gsm_file}: {classified.shape[0]} samples -> {out_file}")
