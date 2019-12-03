import networkx as nx
import numpy as np

from ..glycan_graph_methods import generate_digraph_from_glycan_string
from ..glycan_plotting import set_x_positions, set_y_positions



def cfg_to_kcf(glycan, glycan_id=''):
    '''Convert a glycan string in CFG format to KCF format
    
    Args:
        glycan (str): A glycan string in CFG format.
    Returns:
        str: A glycan in KCF format.
    '''
    xscale = 5
    yscale = 10
    glycan_graph = generate_digraph_from_glycan_string(glycan, parse_linker=True)
    y_positions = set_y_positions(glycan_graph)
    x_positions = set_x_positions(glycan_graph)
    y_mid = np.max(list(y_positions.values())) - np.ptp(list(y_positions.values())) / 2
    node_dict = nx.get_node_attributes(glycan_graph, 'label')
    edge_dict = nx.get_edge_attributes(glycan_graph, 'label')
    node_str = '\n'.join([f"      {key}  {value}    {x_positions[key] * xscale:.0f}    {(y_positions[key] - y_mid) * yscale:.0f}" \
        for key, value in sorted(node_dict.items())])
    edge_str = '\n'.join([f"      {i}  {x[0][1]}:{''.join(x[1][0:2])}   {x[0][0]}:{x[1][2]}" for i, x in enumerate(edge_dict.items())])
    kcf_string = (f'ENTRY         Glycan\n'
                  f'NODE  {len(node_dict)}\n'
                  f'{node_str}\n'
                  f'EDGE   {len(edge_dict)}\n'
                  f'{edge_str}\n'
                  f'///\n'
                 )
    return kcf_string


                          
def kcf_to_digraph(kcf_glycan):
    '''Converts a glycan string in KCF format to DiGraph format.
    
    Args:
       kcf_glycan (str): A glycan in KCF string format.
    Returns:
        networkx.DiGraph: A glycan in DiGraph format.
    '''
    lines = kcf_glycan.split('\n')
    node_dict = {}
    edge_dict = {}
    for line in lines:
        line = line.strip()
        if line.startswith('ENTRY'):
            pass
        elif line.startswith('NODE'):
            node_flag = True
        elif line.startswith("EDGE"):
            node_flag = False
            edge_flag = True
        elif not line:
            pass
        elif line.startswith('///'):
            break
        elif node_flag:
            node_id, node_label, _, _ = line.split()
            node_dict[int(node_id)] = {'label': node_label}
        elif edge_flag:
            _, edge_child, edge_parent = line.split()
            child_id = int(edge_child.split(':')[0])
            parent_id = int(edge_parent.split(":")[0])
            try:
                parent_link = edge_parent.split(':')[1]
            except IndexError:
                parent_link = ''
            try:
                child_link = ''.join([x for x in edge_child.split(":")[1] if x.isnumeric()])
            except IndexError:
                child_link = ''
            try:
                child_anomer = ''.join([x for x in edge_child.split(":")[1] if x.isalpha()])
            except IndexError:
                child_anomer = ''
            edge_dict[(parent_id, child_id)] = {'label': (child_anomer, child_link, parent_link)}
    G = nx.DiGraph()
    G.add_nodes_from(node_dict.keys())
    G.add_edges_from(edge_dict.keys())
    nx.set_edge_attributes(G, edge_dict)
    nx.set_node_attributes(G, node_dict)
    return G
            