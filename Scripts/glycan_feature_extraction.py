import networkx as nx
import numpy as np
import glycan_graph_methods

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
        root_node = glycan_graph_methods.find_root_node(G)
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




def create_new_features_for_root_crowdedness(feature_df, feature_set_by_mrmr, frequent_subtrees, cfg_graphs_null_nodes):
    '''Get a list of new features based on the number of siblings of root node in subtree.
    
    Features are split into a single sibling (1) or less than three siblings (<3).
    
    More than 2 siblings is implicitly captured by the other sibling features.
    
    In glycan terms, more than 2 siblings typically refers to a bisecting feature (i.e. bisecting GlcNAc).
    
    Args:
        feature_df (DataFrame): A feature dataframe that contains all extracted features.
        selected_features (list): A list of feature column ids for currently selected features.
        frequent_subtrees (list): A list of frequent subtrees with indices corresponding to
            feature ids.
    Returns:
        DataFrame: An updated feature DataFrame containing all old features plus new 
            features based on parent linkage type.
        list: A list of all frequent subtrees, where indexes match column labels in returned 
            feature DataFrame.
    '''
    feature_df_new = feature_df.copy()
    extended_features = []

    nm = nx.isomorphism.categorical_node_match('label', None)
    em = nx.isomorphism.categorical_edge_match('label', None)

    for selected_feature in feature_set_by_mrmr:
        new_siblings_less_than_two_feature = np.zeros((1, len(cfg_graphs_null_nodes)), dtype='int')
        new_siblings_less_than_three_feature = np.zeros((1, len(cfg_graphs_null_nodes)), dtype='int')

        G = frequent_subtrees[int(selected_feature)]['subtree']
        root_node = glycan_graph_methods.find_root_node(G)

        for main_graph_id in frequent_subtrees[int(selected_feature)]['parents']:
            main_graph = cfg_graphs_null_nodes[main_graph_id]

            GM = nx.isomorphism.DiGraphMatcher(main_graph, G, node_match=nm, edge_match=em)
            has_less_than_two_siblings, has_less_than_three_siblings = False, False
            node_labels = nx.get_node_attributes(main_graph, 'label')
            for match in GM.subgraph_isomorphisms_iter():
                reverse_match = {v: k for k, v in match.items()}
                root_in_main_graph = reverse_match[root_node]
                siblings = glycan_graph_methods.get_siblings(main_graph, root_in_main_graph)
                #Ignore null nodes when finding siblings.
                siblings = {x for x in siblings if node_labels[x] != 'null'}
                if len(siblings) < 2:
                    has_less_than_two_siblings = True
                if len(siblings) < 3:
                    has_less_than_three_siblings = True

            new_siblings_less_than_two_feature[0, main_graph_id] = int(has_less_than_two_siblings)
            new_siblings_less_than_three_feature[0, main_graph_id] = int(has_less_than_three_siblings)

        new_column_label = str(int(feature_df_new.columns[-1]) + 1)
        feature_df_new[new_column_label] = new_siblings_less_than_two_feature[0, feature_df_new.index]
        new_column_label = str(int(feature_df_new.columns[-1]) + 1)
        feature_df_new[new_column_label] = new_siblings_less_than_three_feature[0, feature_df_new.index]

        extended_features.append({'parents': list(np.where(new_siblings_less_than_two_feature == 1)[1]),
                                  'subtree': G,
                                  'sibling_count': '1'})
        extended_features.append({'parents': list(np.where(new_siblings_less_than_three_feature == 1)[1]),
                                  'subtree': G,
                                  'sibling_count': '<3'})
    return feature_df_new, frequent_subtrees + extended_features
