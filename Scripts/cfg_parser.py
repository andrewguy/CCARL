from collections import defaultdict
import pandas as pd
import pyparsing as pp
from pyparsing import Word, Optional, nums, oneOf, Suppress, Literal, OneOrMore, Or
from pyparsing import Group, Forward, And, LineEnd

'''
This module contains a number of functions to convert a CFG glycan string into a more friendly
format.

The main class within this module is the `CFGGlycanParser` class:

>>> glycan_string = '(6S)Galb1-4(Fuca1-3)(Fuca1-5)Ara-olb1-6Galb1-4Glc-Sp21'
>>> parser = CFGGlycanParser()
>>> results = parser.cfg_string_to_graph(glycan_string)

The output from this function can be used to build graphs using one of several Python graph
libraries, such as `networkx`.
'''

# Spacer anomer types taken from Grant et al., 2014 (Glycobiology)
SPACER_ANOMERS = {
    'Sp0': 'a/b',
    'Sp8': 'a/b',
    'Sp9': 'a/b',
    'Sp10': 'a/b',
    'Sp11': 'a',
    'Sp12': 'a/b',
    'Sp13': 'b',
    'Sp14': 'a/b',
    'Sp15': 'a',
    'Sp16': 'a',
    'Sp17': 'b',
    'Sp18': 'b',
    'Sp19': 'a/b',
    'Sp20': 'b',
    'Sp21': 'a/b',
    'Sp22': 'b',
    'Sp23': 'b',
    'Sp24': 'a/b',
    'Sp25': 'a/b'
}

class CFGGlycanParser(object):
    '''A parser for CFG glycan strings.

    Provides access to a parser object that can be used to parse CFG glycan strings.
    '''
    def __init__(self, monosaccharide_codes=None):
        '''Create a parser object which can then be used to parse multiple CFG glycan strings.

        Args:
            monosaccharide_codes (list, optional): A list of monosaccharide codes to initialise
                the parser with. If `None`, then uses the KEGG list of glycan codes, with the
                addition of `G-ol`, which appears in the CFG glycan arrays but is not defined
                by the KEGG codes.
        '''
        if monosaccharide_codes is None:
            self.monosaccharide_codes = get_default_monosaccharide_codes()
        self.parser = self._create_cfg_parser()

    def cfg_string_to_graph(self, glycan_string, parse_linker=False):
        '''Convert a CFG glycan string into a simple graph object.

        Args:
            glycan_string: A string containing a CFG glycan in extended IUPAC condensed format.
        Returns:
            defaultdict: Adjacency dictionary.
            dict: Node labels.
            dict: Vertex labels.
        '''
        parsed_result = self.parser.parseString(glycan_string)
        results = parse_cfg_structure_to_graph(parsed_result, parse_linker=parse_linker)
        return results

    def _create_cfg_parser(self):
        '''Create a parser object using the pyparsing module.

        CFG Glycans are encoded using an 'extended IUPAC condensed' format.

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
        spacer ::= 'Sp' digit+
        glycan ::= ( [('(' glycan ')')+] sugar_and_link+ )+ [spacer]
        '''
        digit = Word(nums, exact=1)
        anomeric_carbon_configuration = oneOf('a b')('Anomeric_carbon_configuration')
        link = Optional(anomeric_carbon_configuration) + Optional(digit)('c1') + Suppress('-') + \
               Optional(digit)('c2')
        modifier = oneOf('S P')("Modifier")
        branched_modification = Suppress('(') + digit + modifier + Suppress(')')
        modification = Group((digit + modifier) | branched_modification)
        monosaccharide = oneOf(self.monosaccharide_codes)('Monosaccharide')

        sugar_and_link = Group(Optional(OneOrMore(modification))("Modification")) + \
                         Group(monosaccharide + link)
        spacer = Group((Literal('Sp') + OneOrMore(digit) + LineEnd()) | Word('-MDPLys'))("Spacer")
        glycan = Forward()

        branch = Group(Suppress('(') + glycan + Suppress(')'))

        glycan << OneOrMore(
            Group(Optional(OneOrMore(branch))('Branch'))
            + OneOrMore(sugar_and_link)
            ) + Optional(spacer)
        return glycan


def get_default_monosaccharide_codes():
    '''Get a list of default monosaccharide codes from list of KEGG codes.'''
    kegg_data = pd.read_csv('../Data/KEGG_Glycan_Codes.csv', header=0)
    monosaccharide_codes = list(kegg_data['Code'])
    # Manually add G-ol - probably Glucose alcohol, which appears in CFG arrays,
    # but not in KEGG list
    monosaccharide_codes.append('G-ol')
    monosaccharide_codes.append('KDN')
    monosaccharide_codes.append('Neu5,9Ac2')
    monosaccharide_codes.append('GlcN(Gc)b')
    monosaccharide_codes.append('GlcNA')

    return monosaccharide_codes


def is_sugar(parser_object):
    '''Returns True if object contains a monosaccharide.'''
    if 'Monosaccharide' in parser_object:
        return True
    else:
        return False


def is_branch(parser_object):
    '''Returns True if object contains a branch'''
    if 'Branch' in parser_object:
        return True
    else:
        return False


def is_modification(parser_object):
    '''Returns True if object contains a modification.'''
    if 'Modification' in parser_object:
        return True
    else:
        return False

def is_spacer(parser_object):
    '''Returns True if object is a spacer.'''
    if parser_object.getName() == 'Spacer':
        return True
    else:
        return False


def get_sugar(parser_object):
    '''Returns the monosaccharide code for a parsed sugar residue'''
    sugar_code = parser_object['Monosaccharide']
    return sugar_code


def get_modification(parser_object):
    '''Returns the modification for a parsed sugar residue'''
    modifier = parser_object['Modifier']
    return modifier


def get_sugar_linkage(parser_object):
    '''Returns the linkage for a parsed sugar residue'''
    linkage = (parser_object.get('Anomeric_carbon_configuration', ''),
               parser_object.get('c1', ''),
               parser_object.get('c2', ''))
    return linkage


def get_modification_linkage(parser_object):
    '''Returns the linkage for a parsed sugar residue'''
    linkage = ('', '', parser_object[0])
    return linkage


def get_spacer(parser_object):
    '''Returns the spacer string for a parsed spacer.'''
    spacer = ''.join(parser_object)
    return spacer


def parse_cfg_structure_to_graph(parsed_result, parse_linker=False):
    '''Convert parsed CFG glycan into a simple graph object.

    Returns:
        defaultdict: Adjacency dictionary.
        dict: Node labels.
        dict: Vertex labels.
    '''
    index = -1
    output_graph = defaultdict(set)
    node_labels = {}
    _vertex_labels = {}
    _parse_cfg_structure_to_graph(parsed_result, index, output_graph,
                                  node_labels, _vertex_labels, parse_linker=parse_linker)
    vertex_labels = {(parent, child): _vertex_labels[parent] for parent, children in output_graph.items()
                     for child in children if parent < child}
    return output_graph, node_labels, vertex_labels


def _parse_cfg_structure_to_graph(parsed_result, index, output_graph,
                                  node_labels, vertex_labels, parse_linker=False):
    '''This function does most of the heavy lifting with regards to defining
    graph nodes and vertices from the parsed CFG string.

    Returns the index of the last sugar/residue processed.
    '''
    to_link_to_next_parent = []
    for item in parsed_result:
        if is_sugar(item):
            index += 1
            item_index = index
            #link to previous item(s) in structure
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
                index = _parse_cfg_structure_to_graph(branch, index, output_graph,
                                                      node_labels, vertex_labels)
                to_link_to_next_parent.append(index)
        # Treat a spacer the same as a sugar residue, except we don't
        # have a linkage to the next residue.
        if parse_linker and is_spacer(item):
            index += 1
            item_index = index
            #link to previous item(s) in structure
            for child_index in to_link_to_next_parent:
                output_graph[child_index].add(item_index)
                output_graph[item_index].add(child_index)
            node_labels[index] = get_spacer(item)
            to_link_to_next_parent = [item_index]
    return index
