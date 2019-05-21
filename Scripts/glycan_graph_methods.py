import itertools
import networkx as nx
from collections import defaultdict, Counter
from tempfile import NamedTemporaryFile, _get_candidate_names, gettempdir
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import pandas as pd
import subprocess
import glob
import os

from cfg_parser import CFGGlycanParser
from glycan_plotting import draw_glycan_diagram

parser = CFGGlycanParser()

def generate_digraph_from_glycan_string(glycan_string, parse_linker=False):
    '''Generate a Networkx digraph object from a CFG glycan string.

    Nodes are renumbered with the root glycan encoded as 0, and node numbers incremented from there.

    Args:
        glycan_string (string): A CFG glycan string in 'extended IUPAC condensed' format.

    Returns:
        A networkx DiGraph object.
    '''
    _output_graph, node_labels, vertex_labels = parser.cfg_string_to_graph(glycan_string, parse_linker=parse_linker)
    max_node_num = len(node_labels) - 1
    g = nx.DiGraph()
    g.add_nodes_from([(max_node_num - key, {'label': value}) for key, value in node_labels.items()])
    g.add_edges_from([(max_node_num - key[1], max_node_num - key[0], {'label': value})
                      for key, value in vertex_labels.items()])
    # Set depth from root node.
    nx.set_node_attributes(g, nx.shortest_path_length(g, source=0), 'depth')
    # A manual hack to deal with Neu5,9Ac2 residues - convert to a Neu5Ac with 9Ac modification.
    if 'Neu5,9Ac2' in node_labels.values():
        for node in list(g.nodes):
            if nx.get_node_attributes(g, 'label')[node] == 'Neu5,9Ac2':
                new_node = max(g.nodes) + 1
                nx.set_node_attributes(g, {node: {'label': 'Neu5Ac'}})
                g.add_node(new_node, label='Ac')
                g.add_edge(node, new_node, label=('', '', '9'))
    return g


def graph_fingerprint(g):
    '''Generate a fingerprint for a DiGraph that is guaranteed to be the same for
    isomorphic graphs with identical node and edge labels (stored as a 'label' attribute).
    Non-isomorphic graphs may or may not have the same fingerprint.

    This function is used to speed up identification of common subtrees within a set of graphs
    by allowing grouping of potentially isomorphic graphs using the `itertools.groupby` function.

    Args:
        g: A Networkx DiGraph object.
    Returns:
        A tuple of properties.
    '''
    node_labels = nx.get_node_attributes(g, 'label')
    edge_labels = tuple(sorted(nx.get_edge_attributes(g, 'label').values()))
    properties = [(value, node_labels[key], len(list(g.successors(key)))) for key, value in g.degree()]
    properties.sort()
    return (tuple(properties), edge_labels)


def group_by_isomorphism(subtrees, fingerprint=graph_fingerprint, minimum_support=10):
    '''Group a list of graphs by isomorphism.

    Usually would take a list of all subtrees from a set of graphs,
    and returns a list of lists containing all subtrees in isomorphic groups.

    Args:
        subtrees: A list of graphs.

    Returns:
        A list of lists, grouping input graphs in isomorphic groups.
    '''
    sorted_graphs = sorted(subtrees, key=fingerprint)

    isomorphic_list = []
    for _key, group in itertools.groupby(sorted_graphs, key=fingerprint):
        similar_graphs = list(group)
        n = len(similar_graphs)
        if n < minimum_support:
            continue
        if n > 1:
            #Check all graphs against first graph. Those that are isomorphic to first can be
            # grouped and removed from later comparisons. Repeat for number 2 etc.
            while similar_graphs:
                to_pop = [0]
                isomorphic_graphs = []
                for i in range(1, len(similar_graphs)):
                    if subgraphs_are_equivalent(similar_graphs[0], similar_graphs[i]):
                        to_pop.append(i)
                for index in sorted(to_pop, reverse=True):
                    isomorphic_graphs.append(similar_graphs.pop(index))
                if len(isomorphic_graphs) >= minimum_support:
                    isomorphic_list.append(isomorphic_graphs)
        else:
            isomorphic_list.append(similar_graphs)
    return isomorphic_list

def generate_graph_from_glycan_string(glycan_string, parse_linker=False):
    '''Generate a Networkx graph object from a CFG glycan string.

    Nodes are renumbered with the root glycan encoded as 0, and node numbers incremented from there.

    Note:
        Should probably use `generate_digraph_from_glycan_string` instead, as glycans are best
        represented as a directional graph, and most of the methods in this module are built to
        function on DiGraph objects.

    Args:
        glycan_string (string): A CFG glycan string in 'extended IUPAC condensed' format.

    Returns:
        A networkx Graph object.
    '''
    _output_graph, node_labels, vertex_labels = parser.cfg_string_to_graph(glycan_string, parse_linker=parse_linker)
    max_node_num = len(node_labels) - 1
    g = nx.Graph()
    g.add_nodes_from([(max_node_num - key, {'label': value}) for key, value in node_labels.items()])
    g.add_edges_from([(max_node_num - key[0], max_node_num - key[1], {'label': value})
                      for key, value in vertex_labels.items()])
    # Set depth from root node.
    nx.set_node_attributes(g, nx.shortest_path_length(g, target=0), 'depth')
    return g


def subgraphs_are_equivalent(subgraph_1, subgraph_2):
    '''Return True if the subgraphs are isomorphic with the same node and edge labels.

    Compares node and edge 'label' attributes. Other attributes are ignored.

    Args:
        subgraph_1: A networkx Graph object.
        subgraph_2: A networkx Graph object.
    Returns:
        True if graphs are isomorphic with same edge and node labels. False otherwise.
    '''
    nm = nx.isomorphism.categorical_node_match('label', None)
    em = nx.isomorphism.categorical_edge_match('label', None)
    return nx.is_isomorphic(subgraph_1, subgraph_2, node_match=nm, edge_match=em)


def powerset(iterable):
    '''powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)'''
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))


def add_termini_nodes_to_graphs(graph_list):
    '''Adds termini nodes to glycan graphs. Returns a new list of graphs.

    This is used as often the absence of attached sugars at a particular carbon position is important
    for recognition by lectins or antibodies. Absence of attached residues at a position that would
    sometimes have an attached residue is indicated by a labelled 'null' node.

    Args:
        graph_list (list): A list of glycans in DiGraph form.
    Return:
        list of DiGraphs: A list of glycans in DiGraph form.
    '''
    # Create a list of possible connections from carbons within each residue type.
    saccharide_connections = defaultdict(Counter)
    for G in graph_list:
        for node in G.nodes():
            sugar = nx.get_node_attributes(G, 'label')[node]
            connections = [x[1]['label'][2] for x in G[node].items()]
            saccharide_connections[sugar].update(connections)
    # Add null connections for each graph in the graph list.
    # A null connection will be added if there exists a possible connection at that
    # carbon position, but doesn't exist for that residue.
    new_graph_list = []
    for graph in graph_list:
        G = graph.copy()
        for node in list(G.nodes()):
            sugar = nx.get_node_attributes(G, 'label')[node]
            connections = [x[1]['label'][2] for x in G[node].items()]
            null_connections = []
            for possible_connection in saccharide_connections[sugar].keys():
                if possible_connection not in connections:
                    null_connections.append(possible_connection)
            if null_connections:
                next_node_index = max(G.nodes()) + 1
                G.add_node(next_node_index, label='null')
                G.add_edge(node, next_node_index, label=('', '', ','.join((sorted(null_connections)))))
        new_graph_list.append(G)
    return new_graph_list

def convert_from_digraph_to_graph(digraph):
    '''Convert from a DiGraph to a 'pseudo directed graph' as a Graph object

    The returned Graph object contains a new node that replaces each edge in
    the original DiGraph, with node attributes containing both direction flag
    and original edge attributes:

    `'label' = ('_DirectionNode_', old_edge_attributes)`.

    The edges in the returned Graph object have edge attributes 'label',
    with `'label' = 0` indicating the closest original node is a parent node
    and `'label' = 1` indicating the closest original node is a child node.

    Args:
        digraph (DiGraph): A networkx DiGraph.

    Returns: A networkx Graph
    '''
    G = digraph
    G_und = nx.Graph()
    G_und.add_nodes_from(G.nodes(data=True))
    max_node_id = max(G_und.nodes())
    edge_attributes = nx.get_edge_attributes(G, 'label')
    for edge in G.edges():
        max_node_id += 1
        G_und.add_node(max_node_id, label=('_DirectionNode_', edge_attributes[edge]))
        G_und.add_edge(edge[0], max_node_id, label=0)
        G_und.add_edge(max_node_id, edge[1], label=1)
    return G_und

def find_root_node(g, guess=None):
    '''Find a root node for digraph.

    Args:
        g (Networkx Digraph): A digraph with a single root node.
        guess (int): An initial guess for the root node index/id.
    Returns:
        root node (int)'''
    if guess is None:
        guess = list(g.nodes())[0]
    if g.in_degree(guess) == 0:
        root_node = guess
    else:
        predecessor = next(g.predecessors(guess))
        root_node = find_root_node(g, predecessor)
    return root_node


def get_siblings(G, node):
    """
    Get siblings of a node in a DiGraph.
    Args:
        G (DiGraph): A NetworkX DiGraph
        node (int/str): A node id.
    Returns:
        list: A list of all sibling IDs for given node.
    """
    siblings = set()
    parents = list(G.predecessors(node))
    for parent in parents:
        siblings.update(list(G.successors(parent)))
    return siblings


def convert_from_graph_to_digraph(graph):
    '''Convert from a Graph to a DiGraph object

    See `convert_from_digraph_to_graph` for the required format for the input
    Graph object. This function is designed to be used in conjunction with
    the `convert_from_digraph_to_graph` function.

    Args:
        graph (Graph): A networkx Graph.

    Returns: A networkx DiGraph
    '''
    G = graph
    G_dir = nx.DiGraph()
    G_dir.add_nodes_from([x for x in G.nodes(data=True) if x[1]['label'][0] != '_DirectionNode_'])
    edge_nodes = [x for x in G.nodes(data=True) if x[1]['label'][0] == '_DirectionNode_']
    edge_directions = nx.get_edge_attributes(G, 'label')
    for edge in edge_nodes:
        neighbors = list(G.neighbors(edge[0]))
        parent = [neighbor for key, value in edge_directions.items()
                  for neighbor in neighbors if set(key) == {edge[0], neighbor} and value == 0][0]
        child = [neighbor for key, value in edge_directions.items()
                 for neighbor in neighbors if set(key) == {edge[0], neighbor} and value == 1][0]
        G_dir.add_edge(parent, child, label=edge[1]['label'][1])
    return G_dir

def hierarchy_pos(G, root, levels=None, width=1., height=1.):
    '''If there is a cycle that is reachable from root, then this will see infinite recursion.
       G: the graph
       root: the root node
       levels: a dictionary
               key: level number (starting from 0)
               value: number of nodes in this level
       width: horizontal space allocated for drawing
       height: vertical space allocated for drawing'''
    TOTAL = "total"
    CURRENT = "current"
    def make_levels(levels, node=root, currentLevel=0, parent=None):
        """Compute the number of nodes for each level
        """
        if not currentLevel in levels:
            levels[currentLevel] = {TOTAL : 0, CURRENT : 0}
        levels[currentLevel][TOTAL] += 1
        neighbors = list(G.neighbors(node))
        if parent is not None and parent in neighbors:
            neighbors.remove(parent)
        for neighbor in neighbors:
            levels =  make_levels(levels, neighbor, currentLevel + 1, node)
        return levels

    def make_pos(pos, node=root, currentLevel=0, parent=None, vert_loc=0):
        dx = 1/levels[currentLevel][TOTAL]
        left = dx/2
        pos[node] = ((left + dx*levels[currentLevel][CURRENT])*width, vert_loc)
        levels[currentLevel][CURRENT] += 1
        neighbors = list(G.neighbors(node))
        if parent is not None and parent in neighbors:
            neighbors.remove(parent)
        for neighbor in neighbors:
            pos = make_pos(pos, neighbor, currentLevel + 1, node, vert_loc-vert_gap)
        return pos
    if levels is None:
        levels = make_levels({})
    else:
        levels = {l:{TOTAL: levels[l], CURRENT:0} for l in levels}
    vert_gap = height / (max([l for l in levels])+1)
    return make_pos({})

def draw_glycan_graph(G):
    """Draw a simple representation of a glycan from a Digraph.

    Args:
        G (DiGraph): A networkx DiGraph object, with node and edge attributes under 'labels'.
    Returns:
        None
    """
    root_node = find_root_node(G)
    pos = hierarchy_pos(G, root_node)
    edge_labels = {k: v[0] + v[1] + '-' + v[2] for k, v in nx.get_edge_attributes(G, 'label').items()}
    nx.draw(G, pos=pos, with_labels=True, font_weight='bold', labels=nx.get_node_attributes(G, 'label'))
    nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels, rotate=False)
    return

def write_lg(graphs):
    '''Write a LineGraph string for use with gBolt.

    Args:
        graphs (list): A list of networkx Graph objects.
    Returns:
        str: A string to be written to a .lg file for gBolt.
        dict: Mapping of node id to label.
        dict: Mapping of edge id to label.
    '''
    output_lines = []
    node_labels = {x for graph in graphs for x in nx.get_node_attributes(graph, 'label').values()}
    node_label_dict = {label: i for i, label in enumerate(node_labels)}
    edge_labels = {x for graph in graphs for x in nx.get_edge_attributes(graph, 'label').values()}
    edge_label_dict = {label: i for i, label in enumerate(edge_labels)}

    for g_id, graph in enumerate(graphs):
        node_labels_by_nodeid = nx.get_node_attributes(graph, 'label')
        edge_labels_by_edgeid = nx.get_edge_attributes(graph, 'label')
        if "id" in graph.graph:
            output_lines.append("t # %s\n" % graph.graph['id'])
        else:
            output_lines.append("t # %s\n" % g_id)
        nodes_output = []
        for n in graph.nodes():
            if n == '':
                n = '[EMPTY_NODE]'
            nodes_output.append("v %i %s\n" % (n, node_label_dict[node_labels_by_nodeid[n]]))
        # Nodes appear to need to be sorted for gBolt to find all subtrees. Have filed a bug report with gBolt 31/1/2019.
        nodes_output.sort(key=lambda x: int(x.split()[1]))
        
        output_lines.extend(nodes_output)
        
        for edge in graph.edges():
            label = edge_label_dict[edge_labels_by_edgeid[edge]]
            if label is not None:
                output_lines.append("e %i %i %s\n" % (edge[0], edge[1], label))
            else:
                output_lines.append("e %i %i\n" % (edge[0], edge[1]))
    node_id_to_label = {v: k for k, v in node_label_dict.items()}
    edge_id_to_label = {v: k for k, v in edge_label_dict.items()}
    return ''.join(output_lines), node_id_to_label, edge_id_to_label

def remove_dangling_edges_from_lg(data_string, node_id_to_label, edge_id_to_label):
    '''Preprocessing for linegraph returned by gBolt. Removes all dangling edges.
    
    Used so we can perform frequent subtree mining on directed graphs.
    
    Args:
        data_string (str): A string containing data from LineGraph files
            output from gBolt.
        node_id_to_label (dict): A dictionary mapping node ids to labels.
        edge_id_to_label (dict): A dictionary mapping edge ids to labels.
    Returns:
        str: Filtered lg file data from gBolt.
    '''
    direction_vertices = set(str(key) for key, value in node_id_to_label.items() 
                             if value[0] == '_DirectionNode_')
    to_keep = []
    for graph in data_string.split('t')[1:]:
        lines = graph.split('\n')
        vertice_labels = [str(x.split()[-1]) for x in lines if x.startswith('v')]  
        if len(vertice_labels) == 2 * len([x for x in vertice_labels if x in direction_vertices]) + 1:
            to_keep.append(graph)
    return 't' + 't'.join(to_keep)

def read_lg(data_string, node_id_to_label, edge_id_to_label, remove_dangling=True):
    '''Convert a LineGraph file from gspan into a networkx Graph object.

    Args:
        data_string (str): A string containing data from LineGraph files
            output from gBolt.
        node_id_to_label (dict): A dictionary mapping node ids to labels.
        edge_id_to_label (dict): A dictionary mapping edge ids to labels.
    Returns:
        list: A list of dictionaries. Each dictionary contains info on a single subtree,
            with parent ids accessed by `'parent'` and subtree graph accessed by `'subtree'`.
    '''
    if remove_dangling:
        data_string = remove_dangling_edges_from_lg(data_string,
                                                    node_id_to_label,
                                                    edge_id_to_label)
    graph_map = {}
    node_map = {}
    counts = {}
    original_graphs = {}
    graph_id = 0
    for line in data_string.splitlines():
        line = line.strip()
        if line.startswith("t #"):
            node_map = {}
            graph_id += 1
            graph = nx.Graph(id=graph_id)
            graph_map[graph_id] = graph
            counts[graph_id] = line.split(" ")[4]
            node_id = 0
        elif line.startswith("v"):
            parts = line.split(" ")
            label_id = int(parts[2])
            graph_map[graph_id].add_node(node_id, label=node_id_to_label[label_id])
            node_map[parts[1]] = node_id
            node_id += 1
        elif line.startswith("e"):
            parts = line.split(" ")
            label_id = int(parts[3])
            graph_map[graph_id].add_edge(node_map[parts[1]], node_map[parts[2]], label=edge_id_to_label[label_id])
        elif line.startswith("x"):
            original_graphs[graph_id] = [int(x) for x in line.split()[1:]]

    graph_list = [{'subtree': graph_map[i], 'parents': original_graphs[i]} for i in graph_map.keys()]
    return graph_list

def has_dangling_edges(G):
    '''Check if graph has an unconnected edge.

    This is used when filtering results from gBolt, as we
    are representing DiGraphs as Graphs, and have 'edge nodes'
    that may be left on their own (without corresponding parent
    and child nodes attached).
    Args:
        G (Graph): A Networkx Graph objects
    Returns:
        bool: True if graph has dangling edges
    '''
    is_dangling = False
    for node in G:
        label = nx.get_node_attributes(G, 'label')[node]
        if (label[0] == '_DirectionNode_') and (G.degree(node) == 1):
            is_dangling = True
            continue
    return is_dangling

def run_gbolt(graphs, support=0.05, gbolt_path='gbolt'):
    '''Extract common subtrees using gBolt.

    Args:
        graphs (list): A list of networkx Graph objects.
        support (float, optional): minimum support for returned subtrees.
        gbolt_location (str): Path to gBolt executable.
    Returns:
        list: A list of dictionaries. Each dictionary contains info on a single subtree,
            with parent ids accessed by `'parent'` and subtree graph accessed by `'subtree'`.
    '''
    lg_str, node_id_to_label, edge_id_to_label = write_lg(graphs)
    with NamedTemporaryFile(mode='w', suffix='.lg') as f:
        f.write(lg_str)
        f.flush()
        input_filename = f.name
        output_name = os.path.join(gettempdir(), next(_get_candidate_names()))
        cmd = '{gbolt} -i {input_filename} -s {support} '\
              '-o {out}  -d '.format(gbolt=gbolt_path,
                                                              input_filename=input_filename,
                                                              support=str(support),
                                                              out=output_name
                                                              )
        output_file_pattern = output_name + ".t*"
        process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
        output, _ = process.communicate()
        output_files = glob.glob(output_file_pattern)
        output_data = []
        for file in output_files:
            with open(file, 'r') as f_out:
                output_data.append(f_out.read())
            os.remove(file)
    output = read_lg('\n'.join(output_data), node_id_to_label, edge_id_to_label)
    return output

def get_frequent_subtrees(digraph_list, support=0.05, gbolt_path='gbolt'):
    '''Extract frequent subtrees from a glycan DiGraph using gBolt.

    Args:
        digraph_list (list): A list of glycan DiGraph objects to extract subtrees from.
        support (float, optional): Minimum support for returned subtrees (between 0 and 1).
        gbolt_path (str): Path to gbolt executable.
    Returns:
        list: A list of dictionaries with keys 'subtree' and 'parents'. Subtree objects are
            found under 'subtree' key, while the id of parent trees that contain that
            subgraph are given in a list under the 'parents' key.
    '''
    graphs = [convert_from_digraph_to_graph(x) for x in digraph_list]
    gbolt_results = run_gbolt(graphs, support=support, gbolt_path=gbolt_path)
    for result in gbolt_results:
        result['subtree'] = convert_from_graph_to_digraph(result['subtree'])
    return gbolt_results

def run_fast_mrmr(feature_df, mrmr_bin_dir, mrmr_reader, num_features=10):
    '''Run fast-mRMR algorithm on a feature dataframe.

    Feature dataframe should contain class label in column 0, with features in every other column.

    Args:
        feature_df (DataFrame): A DataFrame containing class labels and features.
        mrmr_bin_dir (str): The full path to the directory containing the mRMR binary.
        mrmr_reader (str): The full path of the mRMR data reader binary.
        num_features (int, optional): The number of features to return from mRMR algorithm.
            Defaults to 10 features.

    Return:
        DataFrame: A filtered subset of the original DataFrame containing selected features.
    '''
    with NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
        feature_df.to_csv(path_or_buf=f, index=False)
        name = f.name
    os.system('cd /tmp && {mrmr_reader} {name}'.format(name=name, mrmr_reader=mrmr_reader))
    output = subprocess.check_output('cd {mrmr_dir} && ./fast-mrmr -f /tmp/data.mrmr -a {n}'.format(
                 n=num_features + 1, mrmr_dir=mrmr_bin_dir), shell=True)
    #os.remove(name)
    os.remove('/tmp/data.mrmr')
    features = [int(x) for x in output.decode().strip().split(',')]
    return list(feature_df.iloc[:, features].columns.values)

def plot_svc_coefficients_and_subtrees(classifier, feature_names, frequent_subtrees,
                                       top_features=5, positive_color='blue', negative_color='red'):
    '''Plot the top positive and negative coefficients from a Support Vector Classifier.

    Also plots relevant glycans for each position.

    Args:
        classifier: A SVC classifier object from sklearn.
        feature_names (np.array): A list of feature names for the SVC features.
        frequent_subtrees (dict): A dictionary of frequent subtrees returned by the `get_frequent_subtrees` method.
        top_features (int): Number of top features to plot for both positive and negative coefficents.
    Returns:
        matplotlib.Figure: A matplotlib Figure object.
    '''
    coef = classifier.coef_.ravel()
    top_positive_coefficients = np.argsort(coef)[-top_features:]
    top_negative_coefficients = np.argsort(coef)[:top_features]
    top_coefficients = np.hstack([top_negative_coefficients, top_positive_coefficients])
    fig, ax = plt.subplots(figsize=(30, 10))
    # Generate glycan plots
    for i, feature in enumerate(top_coefficients):
        if coef[feature] > 0:
            y_pos = -1.0
        else:
            y_pos = 1.0
        x_pos = i
        inset_ = inset_axes(ax,
                        width=2,                     # inch
                        height=2,                    # inch
                        bbox_transform=ax.transData, # relative axes coordinates
                        bbox_to_anchor=(x_pos,y_pos),    # relative axes coordinates
                        loc=10)                       # loc=lower left corner
        plt.sca(inset_)

        #draw_glycan_graph(frequent_subtrees[int(feature_names[feature])]['subtree'])
        draw_glycan_diagram(frequent_subtrees[int(feature_names[feature])]['subtree'], inset_,
                            draw_terminal_connection_labels=True)
    colors = [negative_color if c < 0 else positive_color for c in coef[top_coefficients]]
    ax.bar(np.arange(2 * top_features), coef[top_coefficients], color=colors)
    feature_names = np.array(feature_names)
    ax.set_xticks(np.arange(0, 2 * top_features))
    ax.set_xticklabels(feature_names[top_coefficients], rotation=90, fontsize='xx-large')
    ax.set_xlabel("Feature ID", fontsize='xx-large')
    ax.set_ylabel("SVC Coefficient", fontsize='xx-large')
    for label in ax.get_yticklabels():
        label.set_fontsize('xx-large')
    return fig

def mad_based_outlier(points, thresh=3.5, two_sided=False):
    '''Identifies outliers using Median Absolute Deviation.
    
    By default, returns outliers who are greater than the main distribution.
    
    Args:
        points (numpy.array): Data points.
        thresh (float): Modified Z-score threshold.
        two_sided (bool, default=False): If True, returns outliers both greater than and
            less than the main distribution.
    Returns:
        numpy.array: A boolean array of outliers.
    '''
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
    signed_diff = np.sum((points - median), axis=-1)
    modified_z_score = 0.6745 * signed_diff / med_abs_deviation
    return modified_z_score > thresh


class GlycanData(object):
    '''Glycan Data object that exposes several methods for plotting and analysing data.'''
    
    def __init__(self, csv_file):
        '''
        Args:
            csv_file (str): A csv file that contains cfg glycan data.
        '''
        self._csv_data = pd.read_csv(csv_file)
        self._bounds = (np.min(self._csv_data['RFU']), 
                        np.max(self._csv_data['RFU']))
        
    def plot_rfu_histogram(self, ax, xrange=None, title=''):
        '''Plot histogram of rfu values.
        '''
        if xrange is None:
            xrange = [self._bounds[0] - 100, self._bounds[1]]
        full_range = self._bounds[1] - self._bounds[0]
        nbins = int(100*full_range/(xrange[1]-xrange[0]))
        ax.hist(self._csv_data['RFU'], bins=nbins)
        ax.set_xlim(*xrange)
        ax.set_title(title)
        ax.set_xlabel('log(RFU)')
        ax.set_ylabel('Counts')
        return ax

    def rfu_values(self):
        '''Return raw RFU values for glycan microarray'''
        return self._csv_data['RFU']

    def log_rfu_values(self, shift_data=True):
        '''Return log transformed RFU values for glycan microarray.
        
        Args:
            shift_data (bool): If true, shifts data so that the minimum value is 1
                (applies a shift of min_value + 1 to each data point) before calculating log of data.
        Returns:
            Array: RFU data, log transformed and optionally scaled.
        '''
        x = self._csv_data['RFU']
        if np.count_nonzero(x) == 0:
            raise ValueError(f"No non-zero values in file.")
        if shift_data:
            x = -min(x) + x + 1
        log_x = np.log10(x)
        return log_x

    def calculate_binders(self, thresholds=(2.0, 2.5)):
        '''Calculate positive and intermediate binders using log transformed data.
        
        Args:
            x (np.array): RFU counts data.
            pos_threshold (float): Modified z-score threshold for positive binders.
            intermediate_threshold (float): Modified z-score threshold for intermediate binders.
        Returns:
            Array: A numpy array of binding classes, where 0 = no binding,
                1 = intermediate binding and 2 = strong binding.
        '''
        x = self.log_rfu_values()
        neg_mask = ~mad_based_outlier(x, thresh=thresholds[0])
        pos_mask = mad_based_outlier(x, thresh=thresholds[1])
        int_mask = ~np.logical_or(neg_mask, pos_mask)
        binding_class_ternary = int_mask * 1 + pos_mask * 2
        return binding_class_ternary

    def plot_log_rfu_histogram(self, ax, title='', thresholds=(2.0, 2.5)):
        '''Plot a histogram of counts and show positive and intermediate binders.
        
        Args:
            ax (AxesSubplot): Matplotlib figure axis.
            title (str, optional): Title for figure.
            thresholds (tuple, optional): Modified z-score thresholds for intermediate and
                positive binders.

        Returns:
            AxesSubplot: Matplotlib figure axis, to allow method chaining.
        '''
        x = self.log_rfu_values()
        binding_class_ternary = self.calculate_binders(thresholds=thresholds)
        binwidth = 0.05
        data_neg = x[binding_class_ternary == 0]
        data_pos = x[binding_class_ternary == 2]
        data_int = x[binding_class_ternary == 1]
        ax.hist(data_neg, bins=np.arange(min(data_neg), max(data_neg) + binwidth, binwidth))
        # Sometimes there aren't any intermediates or positives...
        if len(data_pos) != 0:
            ax.hist(data_pos, bins=np.arange(min(data_pos), max(data_pos) + binwidth, binwidth),
                    color='red')
        if len(data_int) != 0:
            bins = np.arange(min(data_int), max(data_int) + binwidth, binwidth)
            ax.hist(data_int, bins=bins, color='orange')
        ax.set_title(title)
        ax.set_xlabel('log(RFU)')
        ax.set_ylabel('Counts')
        return ax