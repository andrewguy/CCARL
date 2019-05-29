import pandas as pd
import numpy as np
import networkx as nx

from .glycan_graph_methods import generate_digraph_from_glycan_string
from .frequent_subtrees import get_frequent_subtrees
from .glycan_graph_methods import find_root_node, add_termini_nodes_to_graphs
from .feature_selection import run_fast_mrmr
from .glycan_graph_methods import graph_fingerprint

def extract_features_from_glycan_graphs(glycan_graphs, binding_class,
    gbolt_path, mrmr_reader, mrmr_bin, support_all=0.05,
    support_pos=0.4, parent_edge_types=False):
    '''Get list of features from an initial list of glycans.
    
    Uses MRMR to generate candidate list of motifs from a DiGraph
    representation with null nodes added.
    
    Extraction of frequent subtrees is performed using gBolt, a fast
    implementation of gSpan.
    
    Args:
        glycan_list (list): A list of glycans as networkX DiGraphs.
        binding_class (np.array): Output classes, with 0 indicating no binding
            and 1 indicating positive binding.
        gbolt_path (str): Path to gbolt executable.
        mrmr_reader_path (str): Path to fast-mrmr reader executable.
        mrmr_bin_dir (str): Path to directory containing fast-mrmr executable.
        support_all (float, optional): Minimum support threshold for finding 
            subtrees in all glycans.
        support_pos (float, optional): Minimum support threshold for finding
            subtrees in positive glycans.
        parent_edge_types (bool): Extend subtrees by adding extra features for
            each parent edge type.
    Returns:
        list: A final set of features selected by MRMR.
        list: A list of dictionaries containing information of a particular
            subtree/motif. Dictionary keys include `'parents'`, which contains
            the index of all parent glycans for that particular subtree,
            `'subtree'`, which contains the DiGraph implementation of that
            subtree, and optionally `'sibling_count'`, which contains the
            number of sibling nodes for the root node of that subtree (inclusive
            of the root node).
         DataFrame: A DataFrame containing all features extracted (including
            those not selected by MRMR)
    '''
    num_mrmr_features = 10
    # Get frequent subtrees using gBolt.
    # Support is set to 5%. Seems to be a reasonable threshold that doesn't miss important subtrees.
    freq_subtrees = get_all_frequent_subtrees(glycan_graphs, binding_class,
                                              support_all=support_all, support_pos=support_pos,
                                              gbolt_path=gbolt_path)
    # We now need to generate a set of input features based on frequent subtrees.
    number_of_glycans = len(glycan_graphs)
    number_of_features = len(freq_subtrees)
    features_ = np.zeros((number_of_glycans, number_of_features), dtype='int')

    for i, freq_subtree in enumerate(freq_subtrees):
        for parent in freq_subtree['parents']:
            features_[parent][i] = 1

    feature_df = pd.DataFrame(features_, columns=[str(i) for i in range(number_of_features)])
    feature_df.insert(0, 'class', binding_class)
    # Select top features using mRMR algorithm, using the fast-mrmr implementation.
    feature_set_by_mrmr = run_fast_mrmr(feature_df, mrmr_bin=mrmr_bin,
        mrmr_reader=mrmr_reader, num_features=num_mrmr_features)
    
    # Don't calculate parent edge types.
    if not parent_edge_types:
        freq_subtrees_subset = [freq_subtrees[int(x)] for x in feature_set_by_mrmr]
        return feature_set_by_mrmr, freq_subtrees_subset, feature_df
    
    feature_df_extended, freq_subtrees_extended = get_parent_edge_type_features(feature_df,
        feature_set_by_mrmr, freq_subtrees, glycan_graphs)
    feature_set_by_mrmr_extended = run_fast_mrmr(feature_df_extended, mrmr_bin=mrmr_bin,
        mrmr_reader=mrmr_reader, num_features=num_mrmr_features)

    freq_subtrees_subset = [freq_subtrees_extended[int(x)] for x in feature_set_by_mrmr_extended]
    return feature_set_by_mrmr_extended, freq_subtrees_subset, feature_df_extended

def generate_features_from_subtrees(frequent_subtree_features, glycan):
    '''Generate a feature set for an individual glycan based on a set of chosen subtree features.
    
    Args:
        frequent_subtree_features (list): A list of frequent subtree features.
        glycan (DiGraph): A candidate glycan in nx.DiGraph format (with null nodes).
    Returns:
        np.array: An array of features, corresponding to provided subtrees.
    '''
    main_graph = glycan
    features = []
    for subtree in frequent_subtree_features:
        nm = nx.isomorphism.categorical_node_match('label', None)
        em = nx.isomorphism.categorical_edge_match('label', None)

        link_types = {'alpha': 'a', 'beta': 'b', 'unknown_link': ''}

        parent_link = None
        # Check if parent edge type is in subtree.
        for node_id, label in nx.get_node_attributes(subtree, 'label').items():
            if label in link_types:
                parent_link = link_types[label]
                subtree = subtree.copy()
                subtree.remove_node(node_id)
                break

        root_node = find_root_node(subtree)

        is_subtree = 0
        GM = nx.isomorphism.DiGraphMatcher(main_graph, subtree, node_match=nm, edge_match=em)
        if parent_link is not None:
            for match in GM.subgraph_isomorphisms_iter():
                reverse_match = {v: k for k, v in match.items()}
                root_in_main_graph = reverse_match[root_node]
                parents = list(main_graph.predecessors(root_in_main_graph))
                if not parents and parent_link == link_types['unknown_link']:
                    is_subtree = 1
                    break
                for parent in parents:
                    parent_edge = nx.get_edge_attributes(main_graph, 'label')[(parent, root_in_main_graph)]
                    if parent_edge[0] == parent_link:
                        is_subtree = 1
                        break
        elif GM.subgraph_is_isomorphic():
            is_subtree = 1
        features.append(is_subtree)
    return np.array(features)


def get_parent_edge_type_features(feature_df, selected_features, frequent_subtrees, glycan_graphs):
    '''Get a list of new features based on the edge type of parent connections.
    
    Edge type will typically be 'alpha', 'beta', or 'unknown_link'.
    
    Args:
        feature_df (DataFrame): A feature dataframe that contains all extracted features.
        selected_features (list): A list of feature column ids for currently selected features.
        frequent_subtrees (list): A list of frequent subtrees with indices corresponding to
            feature ids.
        glycan_graphs (list): A list of all original glycan DiGraphs.
    Returns:
        DataFrame: An updated feature DataFrame containing all old features plus new 
            features based on parent linkage type.
        list: A list of all frequent subtrees, where indexes match column labels in returned 
            feature DataFrame.
    '''
    feature_df_new = feature_df.copy()
    extended_frequent_subtrees = []

    nm = nx.isomorphism.categorical_node_match('label', None)
    em = nx.isomorphism.categorical_edge_match('label', None)
    
    link_types = ('alpha', 'beta', 'unknown_link')
    edge_labels = ('a', 'b', '')

    # Loop over each subtree that was selected by MRMR.
    for selected_feature in selected_features:
        G = frequent_subtrees[int(selected_feature)]['subtree']
        root_node = find_root_node(G)
        next_node_id = max(G.nodes()) + 1
        
        parent_feats = [np.zeros((1, len(glycan_graphs)), dtype='int') for _ in link_types]
        parent_graphs = [G.copy() for _ in link_types]
        
        for parent_graph, link_type in zip(parent_graphs, link_types):
            parent_graph.add_node(next_node_id, label=link_type)
            parent_graph.add_edge(next_node_id, root_node, label=('','',''))

        for main_graph_id in frequent_subtrees[int(selected_feature)]['parents']:
            main_graph = glycan_graphs[main_graph_id]
            GM = nx.isomorphism.DiGraphMatcher(main_graph, G, node_match=nm, edge_match=em)
            has_parent_types = [False for i in link_types]
            # Find location of all matching subtrees in glycan.
            # Collect parent linkage types for each match.
            for match in GM.subgraph_isomorphisms_iter():
                reverse_match = {v: k for k, v in match.items()}
                root_in_main_graph = reverse_match[root_node]
                parents = list(main_graph.predecessors(root_in_main_graph))
                if not parents:
                    # Unknown linkage is always last.
                    has_parent_types[-1] = True
                for parent in parents:
                    parent_edge = nx.get_edge_attributes(main_graph, 'label')[(parent, root_in_main_graph)]
                    for i, edge_label in enumerate(edge_labels):
                        if parent_edge[0] == edge_label:
                            has_parent_types[i] = True
            for parent_feat, has_parent_type in zip(parent_feats, has_parent_types):
                parent_feat[0, main_graph_id] = int(has_parent_type)
        for parent_feat, parent_graph in zip(parent_feats, parent_graphs):
            # Generate a new unique column label.
            new_column_label = str(int(feature_df_new.columns[-1]) + 1)
            # Need to only include rows that have the same index as the original feature_df.
            # Some rows were probably dropped when removing intermediate binders.
            feature_df_new[new_column_label] = parent_feat[0, feature_df_new.index]
            extended_frequent_subtrees.append({'parents': list(np.where(parent_feat == 1)[1]),
                                           'subtree': parent_graph})
    return feature_df_new, frequent_subtrees + extended_frequent_subtrees


def get_all_frequent_subtrees(glycan_graphs, binding_class, support_all=0.05,
    support_pos=0.3, gbolt_path='gbolt'):
    '''Gets all frequent subtrees from both all and only positive graphs.
    
    Uses different thresholds for each.
    
    Args:
        glycan_graphs (list): A list of glycan graphs as nx.DiGraph objects.
        binding_class (np.array): A numpy array of binding classes, where
            0 = no binding, 1 = strong binding.
        support_all (float): Minimum support for subtrees from all graphs.
        support_pos (float): Minimum support for subtrees from positive binding
            graphs.
        gbolt_path (str): Path to gBolt executable.
    Returns:
        list: A list of dictionaries with keys 'subtree' and 'parents'. Subtree
            objects are found under 'subtree' key, while the id of parent trees
            that contain that subgraph are given in a list under the 'parents'
            key.
        '''
    glycan_graphs_pos = [glycan for glycan, binding in 
                         zip(glycan_graphs, binding_class) if binding == 1]
    frequent_subtrees_pos = get_frequent_subtrees(glycan_graphs_pos, 
                                                  gbolt_path=gbolt_path, 
                                                  support=support_pos)
    positive_indices = np.where(binding_class == 1)[0]
    for frequent_subtree in frequent_subtrees_pos:
        frequent_subtree['parents'] = [positive_indices[x] for x in 
                                       frequent_subtree['parents']]
    # Get frequent subtrees using gBolt.
    # Default support is set to 5%.
    # Seems to be a reasonable threshold that doesn't miss important subtrees.
    frequent_subtrees = get_frequent_subtrees(glycan_graphs, 
                                              gbolt_path=gbolt_path, 
                                              support=support_all)
    #Initial filtering to remove duplicates
    graph_fingerprints_main = set(graph_fingerprint(x['subtree']) for x in frequent_subtrees)
    graph_fingerprints_positive = [graph_fingerprint(x['subtree']) for x in frequent_subtrees_pos]
    new_subtrees = []
    for i, fingerprint in enumerate(graph_fingerprints_positive):
        if fingerprint not in graph_fingerprints_main:
            new_subtrees.append(frequent_subtrees_pos[i])
            
    nm = nx.isomorphism.categorical_node_match('label', None)
    em = nx.isomorphism.categorical_edge_match('label', None)
    
    #For new subtrees, find all parents graphs 
    for new_subtree in new_subtrees:
        parents = new_subtree['parents']
        for i, main_graph in enumerate(glycan_graphs):
            if i in parents:
                continue
            GM = nx.isomorphism.DiGraphMatcher(main_graph, new_subtree['subtree'],
                                               node_match=nm, edge_match=em)
            if GM.subgraph_is_isomorphic():
                parents.append(i)
        parents.sort()
        new_subtree['parents'] = parents
    return frequent_subtrees + new_subtrees