#!/usr/bin/env python3

from hmm import get_codon_frequencies
from hmm import get_initial_probabilities
from hmm import get_central_probabilities
from hmm import get_emission_probabilities
from hmm import get_states
from process import read_fasta
from process import read_genes
from process import output_genes
from process import get_complement_sequence
from process import get_metrics
from process import get_mean_genes
from viterbi import viterbi
from rnn import run_rnn
from rnn_utils import tag_sequence, split_labels, output_rnn_genes, prepare_sequence

import argparse


def main():
    parser = argparse.ArgumentParser(
        description='Parse a sequence into gene regions using Viterbi.')
    parser.add_argument('-f', action="store", dest="f", type=str, required=True)
    parser.add_argument('-g', action="store", dest="g", type=str, required=True)
    parser.add_argument('-e', action="store", dest="e", type=str, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    gene_file = args.g
    epochs = int(args.e)

    (codon_frequency, column_names, count) = get_codon_frequencies("codon_frequency.csv")
    central_prob = get_central_probabilities(codon_frequency, 100)
    init_prob = get_initial_probabilities(central_prob, 100)
    start_prob, inter_prob, stop_prob = get_emission_probabilities()
    states = get_states(codon_frequency, start_prob, inter_prob, stop_prob, central_prob)

    """
    for state_name in states:
        print(states[state_name])
    """

    # Normal strand
    obs_sequence = read_fasta(fasta_file)
    expected_genes = read_genes(gene_file, complement=False)
    train_seq_tags = tag_sequence(obs_sequence, expected_genes)
    all_samples = split_labels(train_seq_tags, expected_genes)
    train = all_samples[0:int(len(all_samples)/5*4)]
    test = all_samples[int(len(all_samples)/5*4):]
    tp = 0
    fp = 0
    fn = 0
    exact = 0
    all_genes = []
    all_expected = []

    """
    for short_test, i in zip(test, range(len(test))):
        print("Test", i)
        # Viterbi
        seq, labels = prepare_sequence(short_test, {'A':'A', 'T':'T', 'G':'G', 'C':'C'})
        state_sequence, p = viterbi(seq, init_prob, states)
        genes = output_genes(state_sequence, states)
        expected_genes = output_rnn_genes(labels)
        all_genes.append(genes)
        all_expected.append(expected_genes)
        print("Actual", genes, "Expected", expected_genes)
        tp_i, fp_i, fn_i, e_i = get_metrics(expected_genes, genes)
        tp += tp_i
        fp += fp_i
        fn += fn_i
        exact += e_i
        print("Viterbi: TP:", tp, "FP:", fp, "FN:", fn, "Exact:", exact) 
    
    with open('output.csv', 'w') as f:
        for g, e in zip(all_genes, all_expected):
            f.write(str(g)+"|"+str(e)+"\n")
    """

    """
    print(len(expected_genes))
    mean = get_mean_genes(expected_genes)
    print("Average gene length", mean)
    print("Total gene length", len(expected_genes) * mean)
    print("Sequence length", len(obs_sequence))
    print("Intergene length", len(obs_sequence) - (len(expected_genes) * mean))
    print("num of intergenes", len(expected_genes) + 1)
    print("Mean intergene length", (len(obs_sequence) - (len(expected_genes) * mean))/(len(expected_genes)+1))
    """

    # RNN
    test_predicted = run_rnn(train, test, epochs)
    genes = output_rnn_genes(test_predicted)
    tp, fp, fn, exact = get_metrics(expected_genes, genes)
    print(genes)
    print("RNN TP:", tp, "FP:", fp, "FN:", fn, "Exact:", exact)

    with open('rnn-output.csv', 'w') as f:
        for g, e in zip(all_genes, all_expected):
            f.write(str(g)+"|"+str(e)+"\n")

    """
    # Complement strand
    com_sequence = get_complement_sequence(obs_sequence)
    expected_genes = read_genes(gene_file, complement=True)

    # Viterbi complement strand
    state_sequence, p = viterbi(com_sequence, init_prob, states)
    genes = output_genes(state_sequence, states)
    tp, fp, fn = get_metrics(expected_genes, genes)
    print("Viterbi TP:", tp, "FP:", fp, "FN:", fn)

    # RNN
    train_seq_tags = tag_sequence(com_sequence, expected_genes)
    all_samples = split_labels(train_seq_tags, expected_genes)
    train = all_samples[0:int(len(all_samples)/5*4)]
    test = all_samples[int(len(all_samples)/5*4):]
    test_predicted = run_rnn(train, test, epochs)
    genes = output_rnn_genes(test_predicted)
    tp, fp, fn = get_metrics(expected_genes, genes)
    print("RNN TP:", tp, "FP:", fp, "FN:", fn)
    """

if __name__ == "__main__":
    main()
