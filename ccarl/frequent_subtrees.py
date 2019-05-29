import networkx as nx
import numpy as np
import glob
import subprocess
from tempfile import NamedTemporaryFile, _get_candidate_names, gettempdir
import os

from .glycan_graph_methods import graph_fingerprint
from .glycan_graph_methods import convert_from_digraph_to_graph
from .glycan_graph_methods import convert_from_graph_to_digraph


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