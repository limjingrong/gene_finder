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
from viterbi import viterbi
from rnn import run_rnn
from rnn_utils import tag_sequence, split_labels, output_rnn_genes

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
    init_prob = get_initial_probabilities(codon_frequency, 100)
    central_prob = get_central_probabilities(codon_frequency, 100)
    start_prob, inter_prob, stop_prob = get_emission_probabilities()
    states = get_states(codon_frequency, start_prob, inter_prob, stop_prob)

    # Normal strand
    obs_sequence = read_fasta(fasta_file)
    expected_genes = read_genes(gene_file, complement=False)

    # Viterbi
    state_sequence, p = viterbi(obs_sequence, init_prob, states)
    genes = output_genes(state_sequence, states)
    tp, fp, fn = get_metrics(expected_genes, genes)
    print("Viterbi: TP:", tp, "FP:", fp, "FN:", fn)

    # RNN
    train_seq_tags = tag_sequence(obs_sequence, expected_genes)
    all_samples = split_labels(train_seq_tags, expected_genes)
    train = all_samples[0:int(len(all_samples)/5*4)]
    test = all_samples[int(len(all_samples)/5*4):]
    test_predicted = run_rnn(train, test, epochs)
    genes = output_rnn_genes(test_predicted)
    tp, fp, fn = get_metrics(expected_genes, genes)
    print("RNN TP:", tp, "FP:", fp, "FN:", fn)

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


if __name__ == "__main__":
    main()
