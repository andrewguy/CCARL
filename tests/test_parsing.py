from ccarl import ccarl
from ccarl.glycan_parsers import cfg_parser, gsl_parser
from collections import Counter


class TestCFGParsing:
    @classmethod
    def setup_class(cls):
        cls.test_glycan_string_1 = '(6S)Galb1-4(Fuca1-3)(Fuca1-5)Ara-olb1-6Galb1-4Glc-Sp21'
        cls.test_gsl_string_1 = '''
       Fuca-1
            |
Galα-3Galß-4Glc-AO
      |     |
      |     |
 Fucα-2     |
       Fucα-3
        '''
        cls.test_gsl_string_2 = '''
Galα-3Galß-4Glc-AO
      |     |
 Fucα-2     |
       Fucα-3'''
        cls.test_gsl_sulf_string_1 = 'NeuAcα-(S)-6Galß-4Glcß-(S)-Cer36'
        cls.test_gsl_sulf_string_2 = '(6S)NeuAcα-(S)-6Galß-4Glcß-(S)-Cer36'
        cls.test_gsl_sulf_string_3 = 'SU-3GlcAß-3Galß-4Glcß-C30'
        cls.test_gsl_3 = 'ΔUA-4GlcNS-AO'
        cls.test_gsl_comp_mod = '(3-deoxy,3-carboxymethyl)Galß-4Glcß-C30'
        cls.test_gsl_agal = 'aGalß-4Glcß-C30'
        cls.cfg_parser = cfg_parser.CFGGlycanParser()
        cls.gsl_parser = gsl_parser.GSLGlycanParser()
        cls.test_small_branched = '''
                                          Galß-4GlcNAcß-6
                                                        |
                                                        Galß-4GlcNAcß-3Galß-4Glcß-Cer
                                                        |
                                          Galß-4GlcNAcß-3'''
        cls.test_large_branched = '''
                                          Galß-4GlcNAcß-6
                                                        |
                            Galß-4GlcNAcß-6             Galß-4GlcNAcß-3Galß-4Glcß-Cer
                                          |             |
              Galß-4GlcNAcß-6             Galß-4GlcNAcß-3
                            |             |
Galß-4GlcNAcß-6             Galß-4GlcNAcß-3
              |             |
              Galß-4GlcNAcß-3
              |
Galß-4GlcNAcß-3'''
        cls.test_med_branched = '''
                                          Galß-4GlcNAcß-6
                                  GlcNAcß-3             |
                                          |             Galß-4GlcNAcß-3Galß-4Glcß-Cer
                                          |             |
                                          Galß-4GlcNAcß-3
                                          |
                            Galß-4GlcNAcß-3'''
        cls.test_c8_mod = '''(C8 diastereoisomer)NeuAcα-3Galß-4Glcß-Cer36'''
        cls.test_complex_su_mod = '''
SU-2
   |
   ΔUA-4GlcNSα-4IdoAα-4GlcNSα-4IdoAα-4GlcNSα-4IdoAα-4GlcNSα-4IdoAα-4GlcNSα-4IdoAα-4GlcNSα-4IdoAα-4GlcNS-AO
        |       |      |       |      |       |      |       |      |       |      |       |      |
     SU-6    SU-2   SU-6    SU-2   SU-6    SU-2   SU-6    SU-2   SU-6    SU-2   SU-6    SU-2   SU-6'''

    def test_parsing_glycan_string(self):
        parsed_graph, node_labels, vertex_labels = self.cfg_parser.string_to_graph(self.test_glycan_string_1)
        assert isinstance(node_labels, dict)
        assert isinstance(vertex_labels, dict)
        assert len(node_labels) == 7 # Not parsing linker
        assert Counter(node_labels.values()) == Counter(['S', 'Gal', 'Fuc', 'Fuc', 'Ara-ol', 'Gal', 'Glc'])
        assert Counter(vertex_labels.values()) == Counter([('', '', '6'), ('b', '1', '4'), ('a', '1', '3'),
                                                           ('a', '1', '5'), ('b', '1', '6'), ('b', '1', '4')])

    def test_parsing_glycan_string_with_linker(self):
        parsed_graph, node_labels, vertex_labels = self.cfg_parser.string_to_graph(self.test_glycan_string_1, parse_linker=True)
        assert isinstance(node_labels, dict)
        assert isinstance(vertex_labels, dict)
        assert len(node_labels) == 8  # Parsing linker
        assert Counter(node_labels.values()) == Counter(['S', 'Gal', 'Fuc', 'Fuc', 'Ara-ol', 'Gal', 'Glc', 'Sp21'])
        assert Counter(vertex_labels.values()) == Counter([('', '', '6'), ('b', '1', '4'), ('a', '1', '3'),
                                                           ('a', '1', '5'), ('b', '1', '6'), ('b', '1', '4'),
                                                           ('', '', '')])

    def test_parsing_gsl_string(self):
        parsed_graph, node_labels, vertex_labels = self.gsl_parser.string_to_graph(self.test_gsl_string_1)
        assert isinstance(node_labels, dict)
        assert isinstance(vertex_labels, dict)
        assert len(node_labels) == 6  # Parsing linker
        assert Counter(node_labels.values()) == Counter(['Fuc', 'Gal', 'Gal', 'Glc', 'Fuc', 'Fuc'])
        assert Counter(vertex_labels.values()) == Counter([('a', '', '1'), ('a', '', '3'), ('b', '', '4'),
                                                           ('a', '', '2'), ('a', '', '3')])

    def test_parsing_gsl_string_with_linker(self):
        parsed_graph, node_labels, vertex_labels = self.gsl_parser.string_to_graph(self.test_gsl_string_1, parse_linker=True)
        assert isinstance(node_labels, dict)
        assert isinstance(vertex_labels, dict)
        assert len(node_labels) == 7  # Parsing linker
        assert Counter(node_labels.values()) == Counter(['Fuc', 'Gal', 'Gal', 'Glc', 'Fuc', 'Fuc', 'AO'])
        assert Counter(vertex_labels.values()) == Counter([('a', '', '1'), ('a', '', '3'), ('b', '', '4'),
                                                           ('a', '', '2'), ('a', '', '3'), ('', '', '')])

    def test_parsing_gsl_string_2_with_linker(self):
        parsed_graph, node_labels, vertex_labels = self.gsl_parser.string_to_graph(self.test_gsl_string_2, parse_linker=True)
        assert isinstance(node_labels, dict)
        assert isinstance(vertex_labels, dict)
        assert len(node_labels) == 6  # Parsing linker
        assert Counter(node_labels.values()) == Counter(['Gal', 'Gal', 'Glc', 'Fuc', 'Fuc', 'AO'])
        assert Counter(vertex_labels.values()) == Counter([('a', '', '3'), ('b', '', '4'),
                                                           ('a', '', '2'), ('a', '', '3'), ('', '', '')])

    def test_parsing_gsl_string_1_with_S(self):
        parsed_graph, node_labels, vertex_labels = self.gsl_parser.string_to_graph(self.test_gsl_sulf_string_1, parse_linker=True)
        print(node_labels)
        assert isinstance(node_labels, dict)
        assert isinstance(vertex_labels, dict)
        assert len(node_labels) == 6  # Parsing linker
        assert Counter(node_labels.values()) == Counter(['Neu5Ac', 'S', 'S', 'Gal', 'Glc', 'Cer36'])
        assert Counter(vertex_labels.values()) == Counter([('a', '', ''), ('', '', '6'),
                                                           ('b', '', '4'), ('b', '', ''), ('', '', '')])

    def test_parsing_gsl_string_2_with_S(self):
        parsed_graph, node_labels, vertex_labels = self.gsl_parser.string_to_graph(self.test_gsl_sulf_string_2, parse_linker=True)
        assert isinstance(node_labels, dict)
        assert isinstance(vertex_labels, dict)
        assert len(node_labels) == 7  # Parsing linker
        assert Counter(node_labels.values()) == Counter(['S', 'Neu5Ac', 'S', 'S', 'Gal', 'Glc', 'Cer36'])
        assert Counter(vertex_labels.values()) == Counter([('', '', '6'), ('a', '', ''), ('', '', '6'),
                                                           ('b', '', '4'), ('b', '', ''), ('', '', '')])

    def test_parsing_gsl_string_with_UA(self):
        parsed_graph, node_labels, vertex_labels = self.gsl_parser.string_to_graph(self.test_gsl_3, parse_linker=True)
        assert isinstance(node_labels, dict)
        assert isinstance(vertex_labels, dict)
        assert len(node_labels) == 3  # Parsing linker
        assert Counter(node_labels.values()) == Counter(['ΔUA', 'GlcNS', 'AO'])
        assert Counter(vertex_labels.values()) == Counter([('', '', ''), ('', '', '4')])

    def test_parsing_gsl_string_with_complex_modification(self):
        parsed_graph, node_labels, vertex_labels = self.gsl_parser.string_to_graph(self.test_gsl_comp_mod, parse_linker=True)
        assert isinstance(node_labels, dict)
        assert isinstance(vertex_labels, dict)
        assert len(node_labels) == 4  # Parsing linker
        assert Counter(node_labels.values()) == Counter(['3-deoxy,3-carboxymethyl', 'Gal', 'Glc', 'C30'])
        assert Counter(vertex_labels.values()) == Counter([('', '', ''), ('b', '', '4'), ('b', '', '')])

    def test_parsing_gsl_string_with_SU(self):
        parsed_graph, node_labels, vertex_labels = self.gsl_parser.string_to_graph(self.test_gsl_sulf_string_3, parse_linker=True)
        assert isinstance(node_labels, dict)
        assert isinstance(vertex_labels, dict)
        assert len(node_labels) == 5  # Parsing linker
        assert Counter(node_labels.values()) == Counter(['S', 'GlcA', 'Gal', 'Glc', 'C30'])
        assert Counter(vertex_labels.values()) == Counter([('', '', '3'), ('b', '', '3'), ('b', '', '4'), ('b', '', '')])

    def test_parsing_gsl_string_with_agal(self):
        parsed_graph, node_labels, vertex_labels = self.gsl_parser.string_to_graph(self.test_gsl_agal, parse_linker=True)
        assert isinstance(node_labels, dict)
        assert isinstance(vertex_labels, dict)
        assert len(node_labels) == 3 # Parsing linker
        assert Counter(node_labels.values()) == Counter(['aGal', 'Glc', 'C30'])
        assert Counter(vertex_labels.values()) == Counter([('b', '', '4'), ('b', '', '')])

    def test_parsing_gsl_small_branch(self):
        parsed_graph, node_labels, vertex_labels = self.gsl_parser.string_to_graph(self.test_small_branched, parse_linker=True)
        assert isinstance(node_labels, dict)
        assert isinstance(vertex_labels, dict)
        assert len(node_labels) == 9 # Parsing linker
        assert Counter(node_labels.values()) == Counter(['Gal', 'GlcNAc', 'Gal', 'GlcNAc', 'Gal', 'GlcNAc', 'Gal', 'Glc', 'Cer'])
        #assert Counter(vertex_labels.values()) == Counter([('b', '', '4'), ('b', '', '')])

    def test_parsing_gsl_med_branch(self):
        parsed_graph, node_labels, vertex_labels = self.gsl_parser.string_to_graph(self.test_med_branched, parse_linker=True)
        assert isinstance(node_labels, dict)
        assert isinstance(vertex_labels, dict)
        assert len(node_labels) == 12 # Parsing linker
        #assert Counter(node_labels.values()) == Counter(['aGal', 'Glc', 'C30'])
        #assert Counter(vertex_labels.values()) == Counter([('b', '', '4'), ('b', '', '')])


    def test_parsing_gsl_large_branch(self):
        parsed_graph, node_labels, vertex_labels = self.gsl_parser.string_to_graph(self.test_large_branched, parse_linker=True)
        assert isinstance(node_labels, dict)
        assert isinstance(vertex_labels, dict)
        assert len(node_labels) == 21 # Parsing linker
        #assert Counter(node_labels.values()) == Counter(['aGal', 'Glc', 'C30'])
        #assert Counter(vertex_labels.values()) == Counter([('b', '', '4'), ('b', '', '')])

    def test_parsing_gsl_c8(self):
        parsed_graph, node_labels, vertex_labels = self.gsl_parser.string_to_graph(self.test_c8_mod, parse_linker=True)
        assert isinstance(node_labels, dict)
        assert isinstance(vertex_labels, dict)
        assert len(node_labels) == 5 # Parsing linker
        assert Counter(node_labels.values()) == Counter(['C8_diastereoisomer', 'Neu5Ac', 'Gal', 'Glc', 'Cer36'])
        assert Counter(vertex_labels.values()) == Counter([('', '', ''), ('a', '', '3'), ('b', '', '4'), ('b', '', '')])

    def test_parsing_gsl_su_mod(self):
        parsed_graph, node_labels, vertex_labels = self.gsl_parser.string_to_graph(self.test_complex_su_mod, parse_linker=True)
        assert isinstance(node_labels, dict)
        assert isinstance(vertex_labels, dict)
        assert len(node_labels) == 29 # Parsing linker
        #assert Counter(node_labels.values()) == Counter(['C8_diastereoisomer', 'Neu5Ac', 'Gal', 'Glc', 'Cer36'])
        #assert Counter(vertex_labels.values()) == Counter([('', '', ''), ('a', '', '3'), ('b', '', '4'), ('b', '', '')])

