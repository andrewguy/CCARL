"""CCARL is a Python package for analysis and identification of glycan binding motifs.

Carbohydrate Classifier Accounting for Restricted Linkages (CCARL) uses a frequent subtree mining
approach combined with feature selection using mRMR to identify key binding motifs from 
glycan microarray data. CCARL also utilises a glycan representation in which the *absence* of 
potential linkages is explicitly modelled, allowing the tool to distinguish between "terminal"
motifs and otherwise identical motifs within the core of a glycan.
"""

from ccarl.ccarl import CCARLClassifier


__version__ = '1.0.0'
__all__ = ["ccarl", "validate_cfg_structures"]
