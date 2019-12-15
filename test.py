from hmm import State, StateType
from hmm import START_CODON, INTERGENE, STOP_CODON_A, STOP_CODON_B
from hmm import get_central_probabilities, get_states
from hmm import get_initial_probabilities, get_emission_probabilities
from process import output_genes
from process import get_complement_sequence
from process import read_genes
from process import get_metrics
from rnn_utils import tag_sequence
from rnn_utils import split_labels
from rnn_utils import prepare_sequence
from rnn_utils import output_rnn_genes
import unittest
from math import log

class TestHmm(unittest.TestCase):
    def setUp(self):
        self.codon_frequency = dict()
        atg_codon = {"Aa": "Met", "Usage": 1, "Random": 1}
        self.codon_frequency["ATG"] = atg_codon #start codon
        tag_codon = {"Aa": "*", "Usage": 0, "Random": 0}
        self.codon_frequency["TAG"] = tag_codon #stop codon

        self.start_prob, self.inter_prob, self.stop_prob = get_emission_probabilities()
        self.init_prob = get_initial_probabilities(self.codon_frequency,1)


    def test_get_central_probabilities(self):
        central_prob = get_central_probabilities(self.codon_frequency, 1)
        expected_cp = dict()
        expected_cp["ATG0"] = log(1)
        expected_cp[STOP_CODON_A + "0"] = log(0.1)
        expected_cp[STOP_CODON_B + "0"] = log(0.1)
        
        self.assertEqual(central_prob, expected_cp)


    def test_get_states(self):
        states = get_states(self.codon_frequency, self.start_prob, self.inter_prob, self.stop_prob)
        self.assertEqual(13, len(states))
        self.assertIn("ATG0", states)
        self.assertIn("ATG1", states)
        self.assertIn("ATG2", states)
        self.assertIn("SPA0", states)
        self.assertIn("SPA1", states)
        self.assertIn("SPA2", states)
        self.assertIn("SPB0", states)
        self.assertIn("SPB1", states)
        self.assertIn("SPB2", states)
        self.assertIn("STT0", states)
        self.assertIn("STT1", states)
        self.assertIn("STT2", states)
        self.assertIn("INT", states)


    def test_get_genes(self):
        states = get_states(self.codon_frequency, self.start_prob, self.inter_prob, self.stop_prob)
        state_sequence = ["INT", "INT", "INT", "STT0", "STT1", "STT2", "ATG0", "ATG1", "ATG2", "SPA0", "SPA1", "SPA2", "INT", "INT"]
        genes = output_genes(state_sequence, states)

        expected_genes = [(6,8)]
        self.assertEqual(genes, expected_genes)

    
    def test_get_complement_sequence(self):
        sequence = "ATGCA"
        complement = get_complement_sequence(sequence)
        self.assertEqual("TACGT", complement)
        self.assertEqual(sequence, get_complement_sequence(complement))


    def test_read_genes(self):
        genes = read_genes("../test.gene", complement=False)
        expected = [(2,5),(8,10)]
        self.assertEqual(expected, genes)

        genes = read_genes("../test.gene", complement=True)
        expected = [(8,10)]
        self.assertEqual(expected, genes)


    def test_get_metrics(self):
        expected = [(1,2), (4,5), (7,9), (10,12)]
        genes_tested = [(1,2), (5,5), (6,7), (10,11)]

        expected_metrics = (1, 3, 3)
        self.assertEqual(expected_metrics, get_metrics(expected, genes_tested))


    def test_tag_sequence(self):
        sequence = "ABCDEF"
        genes_list = [(1,2), (4,4)]
        labels = tag_sequence(sequence, genes_list)

        expected = [("A", 0), ("B", 1),("C", 1), ("D", 0), ("E", 1), ("F", 0)]
        self.assertEqual(expected, labels)

    def test_split_labels(self):
        labels = [("A", 0), ("B", 1),("C", 1), ("D", 0), ("E", 1), ("F", 0)]
        genes_list = [(1,2), (4,4)]
        split = split_labels(labels, genes_list)
        expected = [[('A', 0), ('B', 1), ('C', 1)], [('D', 0), ('E', 1), ('F', 0)]]
        self.assertEqual(expected, split)

    def test_prepare_sequence(self):
        seq_labels = [("A", 0), ("T", 1),("C", 1), ("T", 0), ("G", 1), ("G", 0)]
        mapping = {'A': [0,0], 'T': [0,1], 'G': [1,0], 'C': [1,1]}
        seq, labels = prepare_sequence(seq_labels, mapping)
        expected_seq = [[0,0], [0,1], [1,1], [0,1], [1,0], [1,0]]
        expected_labels = [0, 1, 1, 0, 1, 0]
        self.assertEqual(expected_seq, seq)
        self.assertEqual(expected_labels, labels)

    def test_output_rnn_genes(self):
        rnn_output = [1,0,0,1,1,0,1]
        expected = [(0,0),(3,4),(6,6)]
        genes = output_rnn_genes(rnn_output)
        self.assertEqual(expected, genes)


if __name__ == "__main__":
    unittest.main()
