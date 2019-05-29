import numpy as np
from sklearn.metrics import matthews_corrcoef, make_scorer
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression

from .glycan_graph_methods import generate_digraph_from_glycan_string, add_termini_nodes_to_graphs
from .glycan_features import extract_features_from_glycan_graphs, generate_features_from_subtrees
from .glycan_graph_methods import get_permitted_connections
from .stats_utils import mad_based_outlier


class CCARLClassifier:
    def __init__(self, support_pos=0.4, support_all=0.05, z_score_thresholds=(1.5, 3.5),
                 permitted_connections=None, mrmr_reader='mrmr-reader',
                 mrmr_bin='fast-mrmr', gbolt_bin='gbolt'):
        self._model = None
        self._features = None
        self._z_score_thresholds = z_score_thresholds
        self._support_all = support_all
        self._support_pos = support_pos
        self._log_rfu = None
        self._rfu = None
        self.glycan_graphs = None
        self.glycan_graphs_with_restriction = None
        self._binding_class = None
        self.subtree_features = None
        self._permitted_connections = permitted_connections
        self._mrmr_reader = mrmr_reader
        self._mrmr_bin = mrmr_bin
        self._gbolt_bin = gbolt_bin
        return

    def train(self, glycans, y, glycan_format="CFG", raw_rfu=True):
        # Optionally process RFU values into tertiary binding classes.
        if raw_rfu:
            self._rfu = y
            self._log_rfu = _log_rfu_values(y, shift_data=True)
            self._binding_class = _calculate_binders(self._log_rfu,
                                                     self._z_score_thresholds)
        else:
            self._binding_class = y
        self.glycan_graphs = [generate_digraph_from_glycan_string(x, parse_linker=True,
                                                                 format=glycan_format)
                              for x in glycans]
        if self._permitted_connections is None:
            self._permitted_connections = get_permitted_connections(self.glycan_graphs)
        self.glycan_graphs_with_restriction = add_termini_nodes_to_graphs(self.glycan_graphs, 
            permitted_connections=self._permitted_connections)
        feature_set_by_mrmr, freq_subtrees_subset, feature_df = extract_features_from_glycan_graphs(
            self.glycan_graphs_with_restriction,
            self._binding_class, self._gbolt_bin, self._mrmr_reader, self._mrmr_bin,
            support_all=self._support_all, support_pos=self._support_pos, 
            parent_edge_types=True)
        mcc_scorer = make_scorer(matthews_corrcoef)

        # Use features to train a classifier...
        # Now perform some additional feature selection with L1 regularisation to reduce unimportant features.
        X = feature_df.loc[:, feature_set_by_mrmr].values
        y = feature_df.iloc[:, 0].astype('bool').astype('int').values
        
        logistic_clf_lasso = LogisticRegressionCV(scoring=mcc_scorer, cv=5, penalty='l1',
                                                  solver='liblinear', class_weight='balanced',
                                                  Cs=100)

        logistic_clf_lasso.fit(X, y)
        coefs = logistic_clf_lasso.coef_[0]
        feature_set_reduced = [x for x, c in zip(feature_set_by_mrmr, coefs) if c != 0]
        freq_subtrees_reduced = [x for x, c in zip(freq_subtrees_subset, coefs) if c != 0]
        if len(feature_set_reduced) == 0:
            print("No features found following L1 regularization!")
            print("Reverting to original features identified by mRMR.")
            feature_set_reduced = feature_set_by_mrmr
            freq_subtrees_reduced = freq_subtrees_subset
        # Run logistic regression again on reduced set of features.
        X = feature_df.loc[:, feature_set_reduced].values
        y = feature_df.iloc[:, 0].astype('bool').astype('int').values

        # Set C to inf to obtain unpenalised logistic regression
        logistic_clf_reduced = LogisticRegression(penalty='l2', C=np.inf, solver='lbfgs',
                                                  class_weight='balanced')
        logistic_clf_reduced.fit(X, y)
        self._model = logistic_clf_reduced
        self._subtrees = [x['subtree'] for x in freq_subtrees_reduced]
        self._training_features = X
        self._training_classes = y
        return

    def predict(self, X, glycan_format="CFG"):
        glycan_graphs = [generate_digraph_from_glycan_string(x, parse_linker=True,
                                                             format=glycan_format)
                         for x in X]
        glycan_graphs_with_restriction = add_termini_nodes_to_graphs(glycan_graphs, 
            permitted_connections=self._permitted_connections)
        features = [generate_features_from_subtrees(self._subtrees, glycan) for 
                    glycan in glycan_graphs_with_restriction]
        return self._model.predict(features)

def _calculate_binders(x, thresholds=(2.0, 2.5)):
    '''Calculate positive and intermediate binders using log transformed data.
    
    Args:
        x (np.array): Log RFU counts data.
        pos_threshold (float): Modified z-score threshold for positive binders.
        intermediate_threshold (float): Modified z-score threshold for intermediate binders.
    Returns:
        Array: A numpy array of binding classes, where 0 = no binding,
            1 = intermediate binding and 2 = strong binding.
    '''
    neg_mask = ~mad_based_outlier(x, thresh=thresholds[0])
    pos_mask = mad_based_outlier(x, thresh=thresholds[1])
    int_mask = ~np.logical_or(neg_mask, pos_mask)
    binding_class_ternary = int_mask * 1 + pos_mask * 2
    return binding_class_ternary


def _log_rfu_values(y, shift_data=True):
    '''Return log transformed RFU values for glycan microarray.
    
    Args:
        shift_data (bool): If true, shifts data so that the minimum value is 1
            (applies a shift of min_value + 1 to each data point) before calculating log of data.
    Returns:
        Array: RFU data, log transformed and optionally scaled.
    '''
    if np.count_nonzero(y) == 0:
        raise ValueError(f"No non-zero values in RFU data.")
    if shift_data:
        y = -min(y) + y + 1
    log_y = np.log10(y)
    return log_y