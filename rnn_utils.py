import random

"""
Tag a FASTA sequence with gene labels.

Arguments:
    sequence: the DNA sequence
    genes_list: A list of tuples (start_inc, end_inc)

Returns:
    List of tuples (DNA_base, 0 or 1)
"""
def tag_sequence(sequence, genes_list):
    labels = []
    genes_index = 0

    prev_end = -1
    for (start_inc, end_inc) in genes_list:
        no_genes_length = start_inc - (prev_end+1)
        assert no_genes_length == len(sequence[prev_end+1:start_inc]), str(no_genes_length) + " " +str(start_inc) + " " + str(prev_end) + " " + str(len(sequence[prev_end+1:start_inc]))
        no_genes = zip(sequence[prev_end+1:start_inc], [0] * no_genes_length)
        no_genes = [(a,b) for (a,b) in no_genes]
        labels.extend(no_genes)

        genes_length = end_inc - start_inc + 1
        genes = zip(sequence[start_inc:end_inc+1], [1] * genes_length)
        genes = [(a,b) for (a,b) in genes]
        labels.extend(genes)

        prev_end = end_inc

    if prev_end < len(sequence) - 1:
        no_genes_length = len(sequence) - (prev_end+1)
        no_genes = zip(sequence[prev_end+1:], [0] * no_genes_length)
        no_genes = [(a,b) for (a,b) in no_genes]
        labels.extend(no_genes)

    assert(len(labels)==len(sequence))

    return labels


"""
Splits a list of labels into samples where each sample contains 1 gene.

Arguments:
    labels: A list of tuples (DNA_base, 0 or 1)
    genes_list: A list of tuples (start_inc, end_inc)

Returns:
    List of list of tuples (DNA_base, 0 or 1)
"""
def split_labels(labels, genes_list):
    if len(genes_list) <= 1:
        return [labels]

    split = []
    start = 0
    gene_end = genes_list[0][1]
    for i in range(1, len(genes_list)):
        next_gene_start = genes_list[i][0]
        # can only start at an intergene
        new_start = random.sample(range(gene_end+1, next_gene_start), 1)[0]
        split.append(labels[start:new_start])
        start = new_start
        gene_end = genes_list[i][1]

    split.append(labels[start:])

    assert(len(split)==len(genes_list))
    return split

"""
Convert sequence with labels to its idx

Arguments:
    labels: A list of tuples (DNA_base, 0 or 1)
    tag_mapping: A dictionary mapping bases to numbers

Returns:
    List of numbers corresponding to each base
    List of numbers corresponding to each label
"""
def prepare_sequence(seq_labels, tag_mapping):
    sequence = []
    labels = []
    for base, label in seq_labels:
        sequence.append(tag_mapping[base])
        labels.append(label)

    return sequence, labels

"""
Computes genes from an RNN output.

Arguments:
    output: List of 0s and 1s.

Returns:
    Tuples of (start_inc, end_inc) for each gene.
"""
def output_rnn_genes(output):
    genes = []
    gene_start = None
    for i in range(len(output)):
        if output[i] == 0 and gene_start is None:
            continue
        elif output[i] == 0 and gene_start is not None:
            gene_end = i-1
            genes.append((gene_start, gene_end))
            gene_start = None
        elif output[i] == 1 and gene_start is None:
            gene_start = i
        elif output[i] == 1 and gene_start is not None:
            continue

    if gene_start is not None:
        genes.append((gene_start, len(output)-1))

    return genes
