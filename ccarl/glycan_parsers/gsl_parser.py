import re
from collections import defaultdict
from itertools import zip_longest

from pyparsing import (Forward, Group, LineEnd, OneOrMore,
                       Optional, Suppress, Word, nums, oneOf)

from .sugar_codes import MONOSACCHARIDE_CODES

SPACERS = {'(Me)AO',
           'AD',
           'AO',
           'C29',
           'C30',
           'Cer',
           'Cer28',
           'Cer32',
           'Cer34',
           'Cer36',
           'Cer42',
           'CerA',
           'CerB',
           'DH',
           'OX',
           'OY',
           'Succ-DH'
           }


MONOSACCHARIDE_CODES_TO_ADD = [
    'aGal',
    'GlcNS',
    'ΔUA',
    'Tyr',
    'Thr',
    'Ser',
    'S'  # Sulfur between sugars will be treated as an independent node in the graph.
]


MODIFICATIONS = [
    '9-OAc',
    '4-OAc',
    '3-deoxy',
    '4-deoxy',
    '6-deoxy',
    '7-deoxy',
    '8-deoxy',
    '9-deoxy',
    '4,8-deoxy',
    '4-OMe',
    '2-OMe',
    '9-OMe',
    '3-deoxy,3-carboxymethyl',
    '4,8-deoxy',
    'C7',
    'C8',
    'C8_diastereoisomer',  # Written as 'C8 diastereoisomer', but the space causes issues with parsing.
    '6-NAc',
    '6-NBz',
    '4,6-deoxy',
    '6-deoxy-6-carboxymethyl'
]


class GSLGlycanParser(object):
    '''A parser for Glycoscience Laboratory (Imperial College London) glycan strings.

    Provides access to a parser object that can be used to parse glycan strings.
    '''
    def __init__(self, monosaccharide_codes=None):
        '''Create a parser object which can then be used to parse multiple GSL glycan strings.

        Args:
            monosaccharide_codes (list, optional): A list of monosaccharide codes to initialise
                the parser with. If `None`, then uses the KEGG list of glycan codes, with the
                addition of `G-ol`, which appears in the CFG glycan arrays but is not defined
                by the KEGG codes. #TODO
        '''
        if monosaccharide_codes is None:
            self.monosaccharide_codes = get_default_monosaccharide_codes()
        self.parser = self._create_gsl_parser()

    def string_to_graph(self, glycan_string, parse_linker=False):
        '''Convert a GSL glycan string into a simple graph object.

        Args:
            glycan_string: A string containing a GSL glycan.
        Returns:
            defaultdict: Adjacency dictionary.
            dict: Node labels.
            dict: Vertex labels.
        '''
        glycan_string = glycosciences_to_cfg(glycan_string)
        parsed_result = self.parser.parseString(glycan_string)
        results = parse_gsl_structure_to_graph(parsed_result, parse_linker=parse_linker)
        return results

    def _create_gsl_parser(self):
        '''Create a parser object using the pyparsing module.

        GSL Glycans are encoded using an 'extended IUPAC condensed' format.

        This parsing module assumes the following EBNF grammar for the
        extended IUPAC condensed format:

        digit ::= "0" | "1" | ... | "9"
        anomeric_carbon_configuration ::= "a" | "b"
        link ::= [anomeric_carbon_configuration] [digit] '-' [digit]
        modifier ::= 'S' | 'P'
        branched_modification ::= '(' digit modifier ')'
        modification ::= (digit modifier) | branched_modification
        monosaccharide ::= "Gal" | "GalNAc" | "Fuc" | "GlcNAc" | ... | "Man" | "KDN"
        sugar_and_link ::= [modification]+ monosaccharide link
        spacer ::= 'AD' | "AO" | "C29" | "C30" | ... | "Cer"
        glycan ::= ( [('(' glycan ')')+] sugar_and_link+ )+ [spacer]
        '''
        digit = Word(nums, exact=1)
        anomeric_carbon_configuration = oneOf('a b')('Anomeric_carbon_configuration')
        link = Optional(anomeric_carbon_configuration) + Optional(digit)('c1') + Suppress('-') + \
            Optional(digit)('c2')
        modifier = oneOf('S P')("Modifier")
        branched_modification = Suppress('(') + digit + modifier + Suppress(')')
        alt_modifiers = Suppress('(') + oneOf(MODIFICATIONS)("Modifier") + Suppress(')')
        modification = Group((digit + modifier) | branched_modification | alt_modifiers)
        monosaccharide = oneOf(self.monosaccharide_codes)('Monosaccharide')

        sugar_and_link = Group(Optional(OneOrMore(modification))("Modification")) + \
            Group(monosaccharide + link)
        spacer = Group((oneOf(SPACERS) + LineEnd()))("Spacer")
        glycan = Forward()

        branch = Group(Suppress('(') + glycan + Suppress(')'))

        glycan << OneOrMore(
            Group(Optional(OneOrMore(branch))('Branch'))
            + OneOrMore(sugar_and_link)
            ) + Optional(spacer)
        return glycan


def get_default_monosaccharide_codes():
    '''Get a list of default monosaccharide codes from list of KEGG codes.'''
    return MONOSACCHARIDE_CODES + MONOSACCHARIDE_CODES_TO_ADD


def is_sugar(parser_object):
    '''Returns True if object contains a monosaccharide.'''
    return 'Monosaccharide' in parser_object


def is_branch(parser_object):
    '''Returns True if object contains a branch'''
    return 'Branch' in parser_object


def is_modification(parser_object):
    '''Returns True if object contains a modification.'''
    return 'Modification' in parser_object


def is_spacer(parser_object):
    '''Returns True if object is a spacer.'''
    return parser_object.getName() == 'Spacer'


def get_sugar(parser_object):
    '''Returns the monosaccharide code for a parsed sugar residue'''
    return parser_object['Monosaccharide']


def get_modification(parser_object):
    '''Returns the modification for a parsed sugar residue'''
    return parser_object['Modifier']


def get_sugar_linkage(parser_object):
    '''Returns the linkage for a parsed sugar residue'''
    linkage = (parser_object.get('Anomeric_carbon_configuration', ''),
               parser_object.get('c1', ''),
               parser_object.get('c2', ''))
    return linkage


def get_modification_linkage(parser_object):
    '''Returns the linkage for a parsed sugar residue'''
    linkage = ('', '', parser_object[0])
    # Some GSL modifications are not preceeded by linkage number.
    if linkage[2] in MODIFICATIONS:
        return ('', '', '')
    return linkage


def get_spacer(parser_object):
    '''Returns the spacer string for a parsed spacer.'''
    spacer = ''.join(parser_object)
    return spacer


def parse_gsl_structure_to_graph(parsed_result, parse_linker=False):
    '''Convert parsed GSL glycan into a simple graph object.

    Returns:
        defaultdict: Adjacency dictionary.
        dict: Node labels.
        dict: Vertex labels.
    '''
    index = -1
    output_graph = defaultdict(set)
    node_labels = {}
    _vertex_labels = {}
    _parse_gsl_structure_to_graph(parsed_result, index, output_graph,
                                  node_labels, _vertex_labels, parse_linker=parse_linker)
    vertex_labels = {(parent, child): _vertex_labels[parent] for parent, children in
                     output_graph.items() for child in children if parent < child}
    return output_graph, node_labels, vertex_labels


def _parse_gsl_structure_to_graph(parsed_result, index, output_graph,
                                  node_labels, vertex_labels, parse_linker=False):
    '''This function does most of the heavy lifting with regards to defining
    graph nodes and vertices from the parsed GSL string.

    Returns the index of the last sugar/residue processed.
    '''
    to_link_to_next_parent = []
    for item in parsed_result:
        if is_sugar(item):
            index += 1
            item_index = index
            # link to previous item(s) in structure
            for child_index in to_link_to_next_parent:
                output_graph[child_index].add(item_index)
                output_graph[item_index].add(child_index)
            node_labels[index] = get_sugar(item)
            vertex_labels[index] = get_sugar_linkage(item)
            to_link_to_next_parent = [item_index]
        if is_modification(item):
            for modification in item:
                index += 1
                item_index = index
                node_labels[index] = get_modification(modification)
                vertex_labels[index] = get_modification_linkage(modification)
                to_link_to_next_parent.append(item_index)
        if is_branch(item):
            for branch in item:
                index = _parse_gsl_structure_to_graph(branch, index, output_graph,
                                                      node_labels, vertex_labels)
                to_link_to_next_parent.append(index)
        # Treat a spacer the same as a sugar residue, except we don't
        # have a linkage to the next residue.
        if parse_linker and is_spacer(item):
            index += 1
            item_index = index
            # link to previous item(s) in structure
            for child_index in to_link_to_next_parent:
                output_graph[child_index].add(item_index)
                output_graph[item_index].add(child_index)
            node_labels[index] = get_spacer(item)
            to_link_to_next_parent = [item_index]
    return index


def _find_branch_end(glycan_lines, branch_point, direction_up=True):
    '''Find the end point of a branch (indicated by '|') in a glycosciences glycan string'''
    if direction_up:
        direction = 1
    else:
        direction = -1
    row = branch_point[0] - direction
    pos = branch_point[1]
    branch_end = glycan_lines[row][pos]
    if branch_end == "|":
        row, pos = _find_branch_end(glycan_lines, (row, pos), direction_up)
    return row, pos


def _identify_branch_insert_points(glycan_string):
    '''Identify branching insert points for a glycan string.

    Return a tuple of ((root line, root position), (branch line, branch position)).
    '''
    branches = []
    ends = []
    lines = glycan_string.split('\n')
    for i, line in enumerate(lines):
        for j, char in enumerate(line):
            if char == '|':
                branches.append((i, j))
    # Need to identify if branches stretch across multiple lines
    branches = sorted(branches)
    for i, branch in enumerate(branches):
        upper = _find_branch_end(lines, branch)
        lower = _find_branch_end(lines, branch, direction_up=False)
        # Identify if an end is a branch or main trunk
        # If upper is a branch, not root, then lines[upper[0]][upper[1] + 1] should be empty or newline.
        # Need to first check that upper[1] + 1 doesn't exceed line length.
        if len(lines[upper[0]]) > (upper[1] + 1) and lines[upper[0]][upper[1] + 1] != ' ':
            ends.append((upper, lower))
        else:
            ends.append((lower, upper))
    return set(ends)


def _identify_branch_end_points(glycan_string):
    '''Identify start and end points for individual branches'''
    branch_end_points = []
    branch_re = re.compile(r'[^\s|]+')
    for i, line in enumerate(glycan_string.split("\n")):
        for match in branch_re.finditer(line):
            branch_end_points.append((i, (match.start(), match.end())))
    return branch_end_points


def _get_branch_index(line_idx, pos_idx, branch_bounds):
    '''Get branch index and insert offset'''
    for i, x in enumerate(branch_bounds):
        if line_idx == x[0] and pos_idx in range(*x[1]):
            return i, pos_idx - x[1][0]
    raise IndexError(f"No glycan branch at line {line_idx} position {pos_idx}")


class _GlycanBranch:
    '''Represents a single glycan branch and its relationships to other branches'''
    def __init__(self, branch_str):
        self.children = []
        self.parent = None
        self.branch = branch_str
        self.parent_insert_index = None

    def get_condensed_representation(self):
        '''Return a one-line glycan representation similar to CFG format.'''
        children_reps = [f"{x.get_condensed_representation()}" for x in self.children]
        # Don't need to add double brackets if not needed. In fact, this breaks the parser for some modifications.
        children_reps = [x if (x[0] == '(' and x[-1] == ')') else f"({x})" for x in children_reps]
        indices = [x.parent_insert_index for x in self.children]
        children_reps = [x for _, x in sorted(zip(indices, children_reps), key=lambda pair: pair[0])]
        indices = sorted(indices)
        indices.insert(0, 0)
        split_rep = [self.branch[i:j] for i, j in zip_longest(indices, indices[1:])]
        rep = ''.join([x for sublist in zip_longest(split_rep, children_reps, fillvalue='') for x in sublist])
        return rep


def _prepare_for_parsing(glycan_str):
    '''Prepare a glycan string for parsing.

    Makes a few changes to ensure some sort of consistency with CFG-style format.

    In particular:

    - Sulfur between groups will be treated as a separate node. Replaces `-(S)-` with `-S-`.
    - All sulfation as a modification will be converted from `SU-2` to `(2S)` format.
    - All phosphate groups converted from P-6 to (6P) format
    - Fix minor issue with Rhα - convert to Rha.
    - All NeuAc and NeuGc converted to Neu5Ac and Neu5Gc to be consistent with CFG.
    - Replace all 'α' and 'ß' with 'a' and 'b'

    Args:
        glycan_str (str): A glycan string in one-line format
    Returns:
        str: A glycan string with minor modifications.
    '''
    glycan_str = glycan_str.replace('-(S)-', '-S-')
    for i in range(1, 7):
        glycan_str = glycan_str.replace(f'SU-{i}', f'{i}S')
    for i in range(1, 7):
        glycan_str = glycan_str.replace(f'P-{i}', f'{i}P')
    glycan_str = glycan_str.replace('SU', 'S')
    glycan_str = glycan_str.replace('Rhα', 'Rha')
    glycan_str = glycan_str.replace('NeuAc', 'Neu5Ac')
    glycan_str = glycan_str.replace('NeuGc', 'Neu5Gc')
    glycan_str = glycan_str.replace('α', 'a')
    glycan_str = glycan_str.replace('ß', 'b')
    return glycan_str


def preprocess_gsl(glycan_str):
    '''Small bit of preprocessing of the GSL string to prevent parsing errors.'''
    glycan_str = glycan_str.replace('\r', '\n')  # Extracted strings from GSL website have \r not \n as line ending.
    glycan_str = glycan_str.replace('C8 diastereoisomer', 'C8_diastereoisomer')  # Space causes issues in parsing.
    return glycan_str


def glycosciences_to_cfg(glycan_str):
    '''Return a condensed representation of a glycan string suitable for parsing.

    Input:
        glycan_str (str): A glycan in the representation used by Imperial College Glycosciences Lab
    Returns:
        str: A condensed representation similar to CFG format.
    '''
    glycan_str = preprocess_gsl(glycan_str)
    end_points = _identify_branch_end_points(glycan_str)
    branch_insert_points = _identify_branch_insert_points(glycan_str)
    branch_strings = [glycan_str.split("\n")[x[0]][x[1][0]:x[1][1]] for x in end_points]
    branches = [_GlycanBranch(x) for x in branch_strings]
    for root_idx, branch_idx in branch_insert_points:
        parent, insert_offset = _get_branch_index(root_idx[0], root_idx[1], end_points)
        child = _get_branch_index(branch_idx[0], branch_idx[1], end_points)[0]
        branches[parent].children.append(branches[child])
        branches[child].parent = branches[parent]
        branches[child].parent_insert_index = insert_offset
    possible_parents = [x for x in branches if x.parent is None]
    assert len(possible_parents) < 2, "Multiple parent/root branches found, error parsing"
    assert len(possible_parents) == 1, "Could not find a parent/root branch, error parsing"
    condensed_rep = possible_parents[0].get_condensed_representation()
    return _prepare_for_parsing(condensed_rep)
