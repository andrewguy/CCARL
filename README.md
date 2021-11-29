# CCARL

A repository for code associated with the "Carbohydrate Classification Accounting for Restricted Linkages" (CCARL) tool.

This tool enables identification of key motifs from glycan microarray data.

Additionally, it can also be used to classify binding to glycans with unknown binding behaviour (given a set of glycan microarray data for a given glycan-binding protein).

## __Installation & Usage__

## Docker Container

The easiest way to run the CCARL tool is to use the provided Docker container which provides easy access to the CCARL command-line scripts:
```bash
docker pull andrewguy/ccarl:v1.0.0

docker run --rm -v `pwd`:/data andrewguy/ccarl:v1.0.0 --help
```

To validate CFG data (often needed because the CFG data usually contains errors which break the glycan parsing tool):

```bash
docker run --rm -v `pwd`:/data andrewguy/ccarl:v1.0.0 validate-structures /data/tests/test_data/ConA_13799-10ug_V5.0_DATA.csv /data/ConA.validated.csv -v
```

To identify binary thresholds for binding/non-binding glycans from RFU data, saving a histogram of binding:

```bash
docker run --rm -v `pwd`:/data andrewguy/ccarl:v1.0.0 identify-binders /data/ConA.validated.csv /data/ConA.binders.csv --histogram /data/ConA_hist.png
```

To identify binding motifs from binary binding data, performing 5-fold cross-validation, plotting ROC curves and saving the model trained on the entire dataset:
```bash
docker run --rm -v `pwd`:/data andrewguy/ccarl:v1.0.0 identify-motifs /data/ConA.binders.csv /data/ConA.results --cross_validation --plot_roc --save_model
```

To perform binding prediction on unknown glycans:

```bash
docker run --rm -v `pwd`:/data andrewguy/ccarl:v1.0.0 predict_binding /data/tests/test_data/test_unknowns.csv /data/ConA.results.model.pkl /data/ConA.predicted.csv
```

## Local Installation

Local installation is useful if you want to integrate the CCARL package into other Python code/scripts, or want to use aspects of the CCARL tool that aren't covered by the command-line scripts.

CCARL **requires** installation of both [`fast-mrmr`](https://github.com/andrewguy/fast-mRMR) and [`gBolt`](https://github.com/Jokeren/gBolt). The binaries for both of these tools should be put in your `$PATH` (otherwise you can specify the binary location when initialising the `CCARLClassifier` object).

Installation of CCARL can be performed by cloning the git repository and running `setup.py`. It is also recommended you install CCARL within a virtual environment:

```bash
git clone https://github.com/andrewguy/CCARL.git && cd CCARL
conda create -n ccarl python=3.8
conda activate ccarl
python setup.py install
```

Note that CCARL **does not** support Python 2.x. If you are still using Python 2.x for any reason, consider this a good time to make the change. :)

To run CCARL from a local installation, you can use the command-line scripts:

```bash
ccarl --help
```

You can also use it as a package from within a Python script/file:

```python
from ccarl import CCARLClassifier
import pandas as pd

csv_data = pd.read_csv("nicely_formatted_glycan_data.csv")

clf = CCARLClassifier()

cf.train(csv_data.Structure, csv_data.Binding.to_numpy(dtype=int))

# Do some more analysis with your trained model here, or run inference on other glycans.
```

## __Citing__

If you use this tool in any of your work, please cite:

Coff, L., Chan, J., Ramsland, P.A., Guy, A.J.  Identifying glycan motifs using a novel subtree mining approach. BMC Bioinformatics 21, 42 (2020). https://doi.org/10.1186/s12859-020-3374-4