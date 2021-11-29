import numpy as np
from sklearn.metrics import matthews_corrcoef, make_scorer
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from statsmodels.stats.outliers_influence import variance_inflation_factor
import warnings

from ccarl.glycan_graph_methods import generate_digraph_from_glycan_string, add_termini_nodes_to_graphs
from ccarl.glycan_features import extract_features_from_glycan_graphs, generate_features_from_subtrees
from ccarl.glycan_graph_methods import get_permitted_connections
from ccarl.stats_utils import mad_based_outlier


class CCARLClassifier:
    def __init__(self, support_pos=0.4, support_all=0.05, num_mrmr_features=10,
                 z_score_thresholds=(1.5, 3.5),
                 permitted_connections=None, mrmr_reader='mrmr-reader',
                 mrmr_bin='fast-mrmr', gbolt_bin='gbolt'):
        self._model = None
        self._features = None
        self._mrmr_subtrees = None
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
        self._num_mrmr_features = num_mrmr_features
        return

    def train(self, glycans, y, glycan_format="CFG", parse_linker=True):
        self._binding_class = y
        self.glycan_graphs = [generate_digraph_from_glycan_string(x, parse_linker=True,
                                                                  format=glycan_format)
                              for x in glycans]
        # Define the allowed linkages for each sugar type
        if self._permitted_connections is None:
            self._permitted_connections = get_permitted_connections(self.glycan_graphs)
        self.glycan_graphs_with_restriction = add_termini_nodes_to_graphs(
            self.glycan_graphs,
            permitted_connections=self._permitted_connections)
        mrmr_features, self._mrmr_subtrees, feature_df = extract_features_from_glycan_graphs(
            self.glycan_graphs_with_restriction,
            self._binding_class, self._gbolt_bin, self._mrmr_reader, self._mrmr_bin,
            support_all=self._support_all, support_pos=self._support_pos,
            parent_edge_types=True, num_mrmr_features=self._num_mrmr_features)

        mcc_scorer = make_scorer(matthews_corrcoef)
        # Use features to train a classifier...
        # Now perform some additional feature selection with L1 regularisation to reduce unimportant features.
        X = feature_df.loc[:, mrmr_features].values
        y = feature_df.iloc[:, 0].astype('bool').astype('int').values

        n_folds = 5
        logistic_clf_lasso = LogisticRegressionCV(scoring=mcc_scorer, cv=n_folds, penalty='l1',
                                                  solver='liblinear', class_weight='balanced',
                                                  Cs=100)
        # MCC function handles NaN values anyway, which is where this warning pops up
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="invalid value encountered in double_scalars")
            logistic_clf_lasso.fit(X, y)
        # For l1 regularisation, C needs to scale inversely with class size.
        # i.e. with more data samples, the best C value is proportionally less.
        best_C_scaled = logistic_clf_lasso.C_[0] * (1 - 1 / n_folds)
        logistic_clf_lasso_final = LogisticRegression(penalty='l1', solver='liblinear',
                                                      class_weight='balanced', C=best_C_scaled)
        logistic_clf_lasso_final.fit(X, y)
        coefs = logistic_clf_lasso_final.coef_[0]

        feature_set_reduced = [x for x, c in zip(mrmr_features, coefs) if c != 0]
        freq_subtrees_reduced = [x for x, c in zip(self._mrmr_subtrees, coefs) if c != 0]
        if len(feature_set_reduced) == 0:
            print("No features found following L1 regularization!")
            print("Reverting to original features identified by mRMR.")
            feature_set_reduced = mrmr_features
            freq_subtrees_reduced = self._mrmr_subtrees

        # Run logistic regression again on reduced set of features.
        X = feature_df.loc[:, feature_set_reduced].values
        y = feature_df.iloc[:, 0].astype('bool').astype('int').values
        if X.shape[1] > 1:
            # Calculate Variance Inflation Factors, and remove redundant features
            # Drop feature if VIF is inf, but do it stepwise.
            # We will handle divide by zero errors later
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="divide by zero encountered in double_scalars")
                vif = [variance_inflation_factor(X, i) for i in range(X.shape[1])]

            # Indices of features to keep.
            to_keep = list(range(len(vif)))

            for _ in range(len(vif)):
                if len(vif) == 1:
                    break
                if np.any(np.array(vif) == np.inf):
                    to_remove = len(vif) - vif[::-1].index(np.inf) - 1
                    del(to_keep[to_remove])
                    if len(to_keep) == 1:
                        break
                    with warnings.catch_warnings():
                        warnings.filterwarnings("ignore", message="divide by zero encountered in double_scalars")
                        vif = [variance_inflation_factor(X[:, np.array(to_keep)], i)
                               for i in range(X[:, np.array(to_keep)].shape[1])]
                else:
                    break

            X = X[:, np.array(to_keep)]
            freq_subtrees_reduced = [freq_subtrees_reduced[x] for x in to_keep]

        # Set C to inf to obtain unpenalised logistic regression
        logistic_clf_reduced = LogisticRegression(penalty='l2', C=100, solver='lbfgs',
                                                  class_weight='balanced')
        logistic_clf_reduced.fit(X, y)
        self._model = logistic_clf_reduced
        self.subtree_features = [x['subtree'] for x in freq_subtrees_reduced]
        self._training_features = X
        self._training_classes = y
        return

    def predict(self, glycans, glycan_format="CFG", parse_linker=True):
        features = self._generate_features(glycans, glycan_format, parse_linker=parse_linker)
        return self._model.predict(features)

    def predict_proba(self, glycans, glycan_format="CFG", parse_linker=True):
        features = self._generate_features(glycans, glycan_format, parse_linker=parse_linker)
        return self._model.predict_proba(features)

    def _generate_features(self, glycans, glycan_format="CFG", parse_linker=True):
        glycan_graphs = [generate_digraph_from_glycan_string(x, parse_linker=parse_linker,
                                                             format=glycan_format)
                         for x in glycans]
        glycan_graphs_with_restriction = add_termini_nodes_to_graphs(glycan_graphs,
                                                                     self._permitted_connections)
        features = [generate_features_from_subtrees(self.subtree_features, glycan) for
                    glycan in glycan_graphs_with_restriction]
        return features


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
        raise ValueError("No non-zero values in RFU data.")
    if shift_data:
        y = -min(y) + y + 1
    log_y = np.log10(y)
    return log_y
