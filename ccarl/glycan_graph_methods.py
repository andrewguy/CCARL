import networkx as nx
from collections import defaultdict

from .glycan_parsers.cfg_parser import CFGGlycanParser
from .glycan_parsers.gsl_parser import GSLGlycanParser

parsers = {
    'CFG': CFGGlycanParser(),
    'GSL': GSLGlycanParser()
}

PERMITTED_CONNECTIONS = {
    'Gal': {6, 3, 4, 2},
    'Glc': {4, 3, 6},
    'Man': {6, 3, 2, 4},
    'GalNAc': {6, 3, 4},
    'Fuc': {3},
    'Neu5Ac': {9, 8},
    'GlcNAc': {4, 6, 3},
    '-MDPLys': {4},
    'GlcA': {3},
    'Neu5Gc': {8},
    'GlcNA': {4}
}


def generate_digraph_from_glycan_string(glycan_string, parse_linker=False,
                                        format='CFG'):
    '''Generate a Networkx digraph object from a CFG glycan string.

    Nodes are renumbered with the root glycan encoded as 0, and node numbers
    incremented from there.

    Args:
        glycan_string (string): A CFG glycan string in 'extended IUPAC
            condensed' format.
        parse_linker(bool, optional): Set to True if linker should be included
            in final graph.
        format ('str'): Format for glycan string. Choose from ['CFG', 'GSL'].

    Returns:
        A networkx DiGraph object.
    '''
    parser = parsers[format]
    _, node_labels, vertex_labels = parser.string_to_graph(glycan_string, parse_linker=parse_linker)
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
    properties = [(value, node_labels[key], len(list(g.successors(key)))) for
                  key, value in g.degree()]
    properties.sort()
    return (tuple(properties), edge_labels)


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


def get_permitted_connections(graph_list):
    '''Return permitted connections for a set of glycan DiGraphs.

    Permitted connections are inferred based on connection types present in
    input glycan list.

    Args:
        graph_list (list): A list of glycans in DiGraph form.
    Return:
        dict: A dictionary of permitted connection sets for each sugar type.
    '''
    # Create a list of possible connections from carbons within each residue type.
    saccharide_connections = defaultdict(set)
    for G in graph_list:
        for node in G.nodes():
            sugar = nx.get_node_attributes(G, 'label')[node]
            connections = [x[1]['label'][2] for x in G[node].items()]
            saccharide_connections[sugar].update(connections)
    return saccharide_connections


def _add_null_connections(graph, permitted_connections):
    '''Add null connections to a graph based on a dictionary of permitted links.

    Args:
        graph (DiGraph): A glycan in DiGraph form.
        permitted_connections (dict, optional): A dictionary of permitted
            connections for each sugar type.
    Returns:
        DiGraph: A copy of the input DiGraph with restricted linkages added.
    '''
    G = graph.copy()
    for node in list(G.nodes()):
        sugar = nx.get_node_attributes(G, 'label')[node]
        connections = {x[1]['label'][2] for x in G[node].items()}
        null_connections = permitted_connections.get(sugar, set()) - connections
        if null_connections:
            next_node_index = max(G.nodes()) + 1
            G.add_node(next_node_index, label='null')
            G.add_edge(node, next_node_index,
                       label=('', '', ','.join((sorted(null_connections)))))
    return G


def add_termini_nodes_to_graphs(graph_list, permitted_connections=PERMITTED_CONNECTIONS):
    '''Adds termini nodes to glycan graphs. Returns a new list of graphs.

    This is used as often the absence of attached sugars at a particular carbon
    position is important for recognition by lectins or antibodies. Absence of
    attached residues at a position that would sometimes have an attached
    residue is indicated by a labelled 'null' node.

    Args:
        graph_list (list): A list of glycans in DiGraph form.
        permitted_connections (dict, optional): A dictionary of permitted
            connections for each sugar type. Set to None if permitted
            connections are to be inferred from connections present in input
            graph list.
    Return:
        list of DiGraphs: A list of glycans in DiGraph form.
    '''
    if permitted_connections is None:
        permitted_connections = get_permitted_connections(graph_list)
    # Add null connections for each graph in the graph list.
    # A null connection will be added if there exists a possible connection at that
    # carbon position, but doesn't exist for that residue.
    new_graph_list = [_add_null_connections(graph, permitted_connections)
                      for graph in graph_list]
    return new_graph_list


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


def digraph_to_glycan_string(G):
    '''Convert a Digraph object into a CFG-like string'''
    root_node = find_root_node(G)
    glycan_str = _digraph_to_string(G, root_node)
    return glycan_str


def _digraph_to_string(G, node):
    node_repr = G.nodes[node]['label']
    children = G.neighbors(node)
    children_reprs = []
    for child in children:
        edge_data = G.edges[(node, child)]['label']
        link_repr = f"{edge_data[0]}{edge_data[1]}-{edge_data[2]}"
        child_repr = _digraph_to_string(G, child)
        children_reprs.append(child_repr + link_repr)
    if children_reprs:
        full_child_repr = children_reprs[0] + ''.join([f"({x})" for x in children_reprs[1:]])
    else:
        full_child_repr = ''
    return f"{full_child_repr}{node_repr}"
