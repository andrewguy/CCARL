from collections import namedtuple, deque, defaultdict
from enum import Enum
import functools

import numpy as np
import matplotlib

from matplotlib.path import Path
from matplotlib.textpath import TextPath
import matplotlib.patches as patches
from matplotlib.transforms import Affine2D
from matplotlib.lines import Line2D

import glycan_graph_methods

import networkx as nx

from math import sqrt, sin, cos, pi

default_line_color = 'black'
default_scale_factor = 1.0

GLYCAN_CODES = [['Hexose', 'Glc', 'Man', 'Gal', 'Gul', 'Alt', 'All', 'Tal', 'Ido', ''],
                ['HexNAc', 'GlcNAc', 'ManNAc', 'GalNAc', 'GulNAc', 'AltNAc', 'AllNAc', 'TalNAc', 'IdoNAc', ''],
                ['Hexosamine', 'GlcN', 'ManN', 'GalN', 'GulN', 'AltN', 'AllN', 'TalN', 'IdoN', ''],
                ['Hexuronate', 'GlcA', 'ManA', 'GalA', 'GulA', 'AltA', 'AllA', 'TalA', 'IdoA', ''],
                ['Deoxyhexose', 'Qui', 'Rha', '', '6dGul', '6dAlt', '', '6dTal', '', 'Fuc'],
                ['DeoxyhexNAc', 'QuiNAc', 'RhaNAc', '', '', '6dAltNAc', '', '6dTalNAc', '', 'FucNAc'],
                ['Di-deoxyhexose', 'Oli', 'Tyv', '', 'Abe', 'Par', 'Dig', 'Col', '', ''],
                ['Pentose', '', 'Ara', 'Lyx', 'Xyl', 'Rib', '', '', '', ''],
                ['Deoxynonulosonate', '', 'Kdn', '', '', '', 'Neu5Ac', 'Neu5Gc', 'Neu', 'Sia '],
                ['Di-deoxynonulosonate', '', 'Pse', 'Leg', '', 'Aci', '', '4eLeg', '', ''],
                ['Unknown', 'Bac', 'LDmanHep', 'Kdo', 'Dha', 'DDmanHep', 'MurNAc', 'MurNGc', 'Mur', ''],
                ['Assigned', 'Api', 'Fru', 'Tag', 'Sor', 'Psi', '', '', '', '']
                ]

GLYCAN_CODE_CONVERSION = {
    "KDN": "Kdn",
    "GlcNA": "GlcNAc"
}

NULL_CODES = ['null']

MODIFICATION_CODES = ['S', 'P', 'Ac', '3S', '4S', '6S', '6P', '9Ac']

CONNECTION_CODES = {'alpha': '\u03B1',
                    'beta': '\u03B2',
                    'unknown_link': '\u03B1/\u03B2'}

def draw_glycan(glycan_code, ax, x, y, scale=0.1, line_weight=1.0, zorder=2):
    '''Draw a glycan symbol on a matplotlib axis.
    
    Returns the current axis.'''
    shape_function = get_glycan_shape(glycan_code)
    glycan_color = get_glycan_color(glycan_code)
    shape_function(ax, x, y, glycan_color, scale=scale, line_weight=line_weight, zorder=zorder)
    return ax

def get_glycan_color(glycan_code):
    '''Return an appropriate color given a glycan code.
    
    Args:
        glycan_code (str): A short glycan code (CFG notation).
    Returns:
        GlycanColor (Enum): A color value for the glycan code.
    '''
    transpose_codes = list(map(list, zip(*GLYCAN_CODES)))
    glycan_color_list = [GlycanColor.WHITE, GlycanColor.BLUE, GlycanColor.GREEN, GlycanColor.YELLOW,
                         GlycanColor.ORANGE, GlycanColor.PINK, GlycanColor.PURPLE, GlycanColor.LIGHT_BLUE,
                         GlycanColor.BROWN, GlycanColor.RED]
    if glycan_code in GLYCAN_CODE_CONVERSION:
        glycan_code = GLYCAN_CODE_CONVERSION[glycan_code]
    if not glycan_code:
        raise ValueError("Must provide a glycan code.")
    if glycan_code in NULL_CODES:
        return GlycanColor.BLACK
    for codes, color in zip(transpose_codes, glycan_color_list):
        if glycan_code in codes:
            return color
    return GlycanColor.WHITE
    #raise ValueError("Glycan code {} not found in standard list.".format(glycan_code))

def get_glycan_shape(glycan_code):
    '''Return an appropriate shape function given a glycan code.
    
    Args:
        glycan_code (str): A short glycan code (CFG notation).
    Returns:
        GlycanColor (Enum): A color value for the glycan code.
    '''
    shape_functions = [draw_circle, draw_square, draw_bisected_square, draw_bisected_diamond, draw_triangle, draw_bisected_triangle,
                       draw_rectangle, draw_star, draw_diamond, draw_flat_diamond, draw_flat_hexagon, draw_pentagon]
    if glycan_code in GLYCAN_CODE_CONVERSION:
        glycan_code = GLYCAN_CODE_CONVERSION[glycan_code]
    if not glycan_code:
        raise ValueError("Must provide a glycan code.")
    if glycan_code in NULL_CODES:
        return draw_cross
    for codes, shape in zip(GLYCAN_CODES, shape_functions):
        if glycan_code in codes:
            return shape
    if glycan_code in MODIFICATION_CODES:
        return draw_modification_text(glycan_code)
    if glycan_code in CONNECTION_CODES:
        return draw_modification_text(CONNECTION_CODES[glycan_code])
    return draw_pentagon
    #raise ValueError("Glycan code {} not found in standard list.".format(glycan_code))


class GlycanColor(Enum):
    WHITE      = (1, 1, 1, 1)
    BLUE       = (0, 144/255, 188/255, 1)
    GREEN      = (0, 166/255, 81/255, 1)
    YELLOW     = (1, 212/255, 0, 1)
    LIGHT_BLUE = (143/255, 204/255, 233/255, 1)
    PINK       = (246/255, 158/255, 161/255, 1)
    PURPLE     = (165/255, 67/255, 153/255, 1)
    BROWN      = (161/255, 122/255, 77/255, 1)
    ORANGE     = (244/255, 121/255, 32/255, 1)
    RED        = (237/255, 28/255, 36/255, 1)
    BLACK      = (0, 0, 0, 1)
    CLEAR      = "None"


def draw_shape(verts_list, ax, x, y, color_list, scale=0.1, line_weight=1.0, zorder=2, edgecolor=GlycanColor.BLACK):
    '''Draw a shape given by a list of vertice coordinate sets.'''
    codes_list = [[Path.MOVETO, *[Path.LINETO]*(len(verts)-2), Path.CLOSEPOLY] for verts in verts_list]
    path_list = [Path(verts * scale, codes) for verts, codes in zip(verts_list, codes_list)]
    trans_path_list = [path.transformed(matplotlib.transforms.Affine2D().translate(x, y)) for path in path_list]
    patch_list = [patches.PathPatch(t_path, facecolor=color.value, lw=line_weight, zorder=zorder,
                  edgecolor=edgecolor.value) for t_path, color in 
                  zip(trans_path_list, color_list)]
    for patch in patch_list:
        ax.add_patch(patch)
    return
    
    
def draw_circle(ax, x, y, color, scale=0.1, line_weight=1.0, zorder=2):
    path = Path(Path.unit_circle().vertices * scale, Path.unit_circle().codes)
    trans = matplotlib.transforms.Affine2D().translate(x, y)
    t_path = path.transformed(trans)
    patch = patches.PathPatch(
        t_path, facecolor=color.value, lw=line_weight, zorder=zorder)
    ax.add_patch(patch)
    return


def draw_square(ax, x, y, color, scale=0.1, line_weight=1.0, zorder=2):
    verts = np.array([
        (0.5, 0.5),
        (0.5, -0.5),
        (-0.5, -0.5),
        (-0.5, 0.5),
        (0.5, 0.5),
        (0., 0.),
    ]) * 2
    draw_shape([verts], ax, x, y, [color], scale=scale, line_weight=line_weight, zorder=zorder)
    return

def draw_bisected_square(ax, x, y, color, scale=0.1, line_weight=1.0, zorder=2):
    top_verts = np.array([
        (0.5, 0.5),
        (0.5, -0.5),
        (-0.5, 0.5),
        (0.5, 0.5),
        (0., 0.),
    ]) * 2
    
    bottom_verts = np.array([
        (0.5, -0.5),
        (-0.5, -0.5),
        (-0.5, 0.5),
        (0.5, -0.5),
        (0., 0.),
    ]) * 2
    draw_shape([top_verts, bottom_verts], ax, x, y, [color, GlycanColor.WHITE], scale=scale, line_weight=line_weight, zorder=zorder)
    return


def draw_bisected_diamond(ax, x, y, color, scale=0.1, line_weight=1.0, zorder=2):
    outside_verts = np.array([
        (1, 0),
        (0, 1),
        (-1, 0),
        (0, -1),
        (1, 0),
        (0, 0)
    ]) * sqrt(2)

    mid_vert = np.array([
        (1,0),
        (-1,0),
        (0, 0)
    ]) * sqrt(2)
    
    top_verts = np.array([
        (1, 0),
        (0, 1),
        (-1, 0),
        (1, 0),
        (0., 0.),
    ]) * sqrt(2)

    bottom_verts = np.array([
        (1, 0),
        (0, -1),
        (-1, 0),
        (1, 0),
        (0., 0.),
    ]) * sqrt(2)
    draw_shape([outside_verts, mid_vert], ax, x, y, [GlycanColor.CLEAR, GlycanColor.CLEAR], scale=scale, line_weight=line_weight, zorder=zorder+1)
    draw_shape([top_verts], ax, x, y, [color], scale=scale, line_weight=line_weight, zorder=zorder, edgecolor=GlycanColor.CLEAR)
    draw_shape([bottom_verts], ax, x, y, [GlycanColor.WHITE], scale=scale, line_weight=line_weight, zorder=zorder, edgecolor=GlycanColor.CLEAR)
    return

def draw_triangle(ax, x, y, color, scale=0.1, line_weight=1.0, zorder=2):
    verts = np.array([
        (0, 0.5),
        (0.625, -0.5),
        (-0.625, -0.5),
        (0, 0.5),
        (0., 0.),
    ]) * 1.9
    draw_shape([verts], ax, x, y, [color], scale=scale, line_weight=line_weight, zorder=zorder)
    return

def draw_bisected_triangle(ax, x, y, color, scale=0.1, line_weight=1.0, zorder=2):
    top_verts = np.array([
        (0, 0.5),
        (0.625, -0.5),
        (0, -0.5),
        (0, 0.5),
        (0., 0.),
    ]) * 1.9
    
    bottom_verts = np.array([
        (0, 0.5),
        (0, -0.5),
        (-0.625, -0.5),
        (0, 0.5),
        (0., 0.),
    ]) * 1.9
    draw_shape([top_verts, bottom_verts], ax, x, y, [color, GlycanColor.WHITE], scale=scale, line_weight=line_weight, zorder=zorder)
    return


def draw_diamond(ax, x, y, color, scale=0.1, line_weight=1.0, zorder=2):
    verts = np.array([
        (0, 1),
        (1, 0),
        (0, -1),
        (-1, 0),
        (0, 1),
        (0., 0.),
    ]) * sqrt(2)
    draw_shape([verts], ax, x, y, [color], scale=scale, line_weight=line_weight, zorder=zorder)
    return

def draw_rectangle(ax, x, y, color, scale=0.1, line_weight=1.0, zorder=2):
    verts = np.array([
        (0.5, 0.25),
        (0.5, -0.25),
        (-0.5, -0.25),
        (-0.5, 0.25),
        (0.5, 0.25),
        (0., 0.),
    ]) * 2
    draw_shape([verts], ax, x, y, [color], scale=scale, line_weight=line_weight, zorder=zorder)
    return

def draw_star(ax, x, y, color, scale=0.1, line_weight=1.0, zorder=2):
    outer_radius = 0.5
    inner_radius = 0.25
    coords_outer = [(outer_radius*cos(-2*pi*i/5 + pi/2), outer_radius*sin(-2*pi*i/5 + pi/2)) for i in range(6)]
    coords_inner = [(inner_radius*cos(-2*pi*i/5 + 7*pi/10), inner_radius*sin(-2*pi*i/5 + 7*pi/10)) for i in range(6)]
    verts = np.array([item for sublist in zip(coords_inner, coords_outer) for item in sublist]) * 2.5
    draw_shape([verts], ax, x, y, [color], scale=scale, line_weight=line_weight, zorder=zorder)
    return

def draw_pentagon(ax, x, y, color, scale=0.1, line_weight=1.0, zorder=2):
    outer_radius = 0.5
    coords_outer = [(outer_radius*cos(-2*pi*i/5 + pi/2), outer_radius*sin(-2*pi*i/5 + pi/2)) for i in range(6)]
    verts = np.array(coords_outer)* 2.5
    draw_shape([verts], ax, x, y, [color], scale=scale, line_weight=line_weight, zorder=zorder)
    return


def draw_flat_diamond(ax, x, y, color, scale=0.1, line_weight=1.0, zorder=2):
    verts = np.array([
        (0, 0.7),
        (1.4, 0),
        (0, -0.7),
        (-1.4, 0),
        (0, 0.7),
        (0., 0.),
    ]) * sqrt(2)
    draw_shape([verts], ax, x, y, [color], scale=scale, line_weight=line_weight, zorder=zorder)
    return

def draw_flat_hexagon(ax, x, y, color, scale=0.1, line_weight=1.0, zorder=2):
    verts = np.array([
        (0.4, 0.5),
        (0.7, 0),
        (0.4, -0.5),
        (-0.4, -0.5),
        (-0.7, 0),
        (-0.4, 0.5),
        (0.4, 0.5),
        (0., 0.),
    ]) * 2
    draw_shape([verts], ax, x, y, [color], scale=scale, line_weight=line_weight, zorder=zorder)
    return

def draw_cross(ax, x, y, color, scale=0.1, line_weight=2.0, zorder=2):
    top_verts = np.array([
        (0.5, 0.5),
        (-0.5, -0.5),
        (0., 0.),
    ]) * 1.5
    
    bottom_verts = np.array([
        (-0.5, 0.5),
        (0.5, -0.5),
        (0., 0.),
    ]) * 1.5
    draw_shape([top_verts, bottom_verts], ax, x, y, [GlycanColor.WHITE, GlycanColor.WHITE], edgecolor=color, 
               scale=scale, line_weight=line_weight, zorder=zorder)
    return

def _draw_modification_text(ax, x, y, color, scale=0.1, line_weight=0, zorder=2, text=''):
    verts = np.array([
        (0.5, 0.5),
        (0.5, -0.5),
        (-0.5, -0.5),
        (-0.5, 0.5),
        (0.5, 0.5),
        (0., 0.),
    ]) * 2
    draw_shape([verts], ax, x, y, [GlycanColor.WHITE], scale=scale, line_weight=0, zorder=zorder)
    ax.text(x, y, text, horizontalalignment='center',
            verticalalignment='center')
    return

def draw_modification_text(text):
    return functools.partial(_draw_modification_text, text=text)



def get_non_null_leaves(G):
    '''Get all non-null leaf nodes (i.e. nodes with no children except for a null node).'''
    leaves = [x for x in list(nx.dfs_preorder_nodes(G)) if not list(nx.neighbors(G, x))]
    new_leaves = []
    for leaf in leaves:
        label = nx.get_node_attributes(G, 'label')[leaf]
        if label == 'null' or label in MODIFICATION_CODES:
            new_leaf = next(G.predecessors(leaf))
            if len(non_null_children(G, new_leaf)) > 0:
                continue
            new_leaves.append(new_leaf)
        else:
            new_leaves.append(leaf)
    return new_leaves

def non_null_children(G, node):
    "Return the number of non-null children nodes."
    children = [x for x in G.successors(node) if nx.get_node_attributes(G, 'label')[x] != "null"
                and nx.get_node_attributes(G, 'label')[x] not in MODIFICATION_CODES]
    return children

    
def set_parent_node_positions(node_queue, G, node, y_positions, y_min_spacing):
    parent = next(G.predecessors(node), None)
    if parent in y_positions or parent is None:
        return
    children = non_null_children(G, parent)
    number_of_children = len(children)
    # If all children haven't been processed, then add parent to queue to be processed later,
    # and return to process other children
    if not set(children).issubset(set(y_positions.keys())):
        node_queue.append(node)
        return
    # If all children have been processed, y position is halfway between maximum and minimum child positions. 
    child_positions = [y_positions[child] for child in children]
    max_child = max(child_positions)
    min_child = min(child_positions)
    y_positions[parent] = (max_child + min_child) / 2
    set_parent_node_positions(node_queue, G, parent, y_positions, y_min_spacing)
    return

def set_y_positions(G):
    leaves = get_non_null_leaves(G)
    q = deque(leaves)
    # Just add one 'null node' spacer for each leaf. :)
    y_min = 0.0
    y_min_spacing = 1.0

    y_positions = {}

    i = 0
    while len(q) > 0:
        leaf = q.popleft()
        if leaf not in y_positions:
            y_positions[leaf] = y_min + (i) * y_min_spacing * 2 + y_min_spacing
            i += 1
        set_parent_node_positions(q, G, leaf, y_positions, y_min_spacing)
    for node in G.nodes():
        if nx.get_node_attributes(G, 'label')[node] == 'null':
            parent = next(G.predecessors(node), None)
            parent_position = y_positions[parent]
            child_position = parent_position - y_min_spacing * 0.6
            y_positions[node] = child_position
    for node in G.nodes():
        label = nx.get_node_attributes(G, 'label')[node]
        if label in MODIFICATION_CODES:
            parent = next(G.predecessors(node), None)
            parent_position = y_positions[parent]
            child_position = parent_position + y_min_spacing * 0.8
            y_positions[node] = child_position
    return y_positions

def set_x_positions(G):
    x_max = 0.0
    x_spacing = 1.0
    depth = nx.shortest_path_length(G, source=glycan_graph_methods.find_root_node(G))
    for node in G.nodes():
        label = nx.get_node_attributes(G, 'label')[node]
        if label == 'null' or label in MODIFICATION_CODES:
            depth[node] = depth[node] - 1
    x_positions = {x: x_max - (d + 1) * x_spacing for x, d in depth.items()}
    return x_positions

def glycan_tree_layout(G):
    x_pos = set_x_positions(G)
    y_pos = set_y_positions(G)
    return {i:(x,y_pos[i]) for i, x in x_pos.items()}


def get_edge_label(G, edge):
    label = nx.get_edge_attributes(G, 'label')[edge]
    return label


def draw_glycan_diagram(G, ax, draw_terminal_connection_labels=False):
    pos = glycan_tree_layout(G)
    #nx.draw_networkx_edges(G, pos, ax=ax)
    node_labels = nx.get_node_attributes(G, 'label')
    for edge in G.edges():
        x_text_offset = 0.01
        y_text_offset = 0.01
        edge_x, edge_y = (pos[edge[0]][0], pos[edge[1]][0]), (pos[edge[0]][1], pos[edge[1]][1])
        edge_label = get_edge_label(G, edge)
        ax.add_line(Line2D(edge_x, edge_y, zorder=1, color='black', linewidth=1))
        if edge_x[0] == edge_x[1]:
            alignment = 'left'
        elif (edge_y[0] - edge_y[1])/(edge_x[0] - edge_x[1]) > 0:
            x_text_offset = 0
            alignment = 'right'
        elif (edge_y[0] - edge_y[1])/(edge_x[0] - edge_x[1]) < 0:
            alignment = 'left'
        else:
            alignment = 'left'
        if edge_label[0] == 'a':
            ab_label = CONNECTION_CODES['alpha']
        elif edge_label[0] == 'b':
            ab_label = CONNECTION_CODES['beta']
        else:
            ab_label = edge_label[0]
        if node_labels[edge[0]] != 'null' and node_labels[edge[1]] != 'null' and \
            node_labels[edge[0]] not in MODIFICATION_CODES and \
            node_labels[edge[1]] not in MODIFICATION_CODES:
            ax.text(np.average(edge_x, weights=(0.3,0.7)) - x_text_offset,
                    np.average(edge_y, weights=(0.3,0.7)) + y_text_offset, ab_label, horizontalalignment=alignment)
            ax.text(np.average(edge_x, weights=(0.6,0.4)) - x_text_offset,
                    np.average(edge_y, weights=(0.6,0.4)) + y_text_offset, edge_label[2], horizontalalignment=alignment)
        if draw_terminal_connection_labels and node_labels[edge[1]] == 'null':
            ax.text(np.average(edge_x, weights=(0.5,0.5)),
                    np.average(edge_y, weights=(0.3,0.7)), edge_label[2], horizontalalignment=alignment)
    #First, process nodes with duplicate positions and combine into one (should only be modification codes that are overlayed)
    reverse_pos_dict = defaultdict(set)
    for node, xy in pos.items():
        reverse_pos_dict[xy].add(node)
    for xy, nodes in reverse_pos_dict.items():
        if len(nodes) > 1:
            names = []
            for node in nodes:
                glycan_name = node_labels[node]
                #Check that this is a modification code
                assert glycan_name in MODIFICATION_CODES
                glycan_edge = nx.get_edge_attributes(G, 'label')[[x for x in G.in_edges(node)][0]]
                glycan_name = str(glycan_edge[2]) + glycan_name
                names.append(glycan_name)
                del pos[node]
            combined_names = ', '.join(names)
            draw_modification_text(combined_names)(ax, xy[0], xy[1], GlycanColor.WHITE, scale=0.03)

    for node, xy in pos.items():
        glycan_name = node_labels[node]
        if glycan_name in MODIFICATION_CODES:
            glycan_edge = nx.get_edge_attributes(G, 'label')[[x for x in G.in_edges(node)][0]]
            glycan_name = str(glycan_edge[2]) + glycan_name
        draw_glycan(glycan_name, ax, xy[0], xy[1], scale=0.2)
    ax.axis('scaled')
    ax.axis('off')
    return ax




