from hmm import StateType
from hmm import INTERGENE, STOP_CODON_A, STOP_CODON_B, START_CODON


'''Reads the fasta file and outputs the sequence to analyze.
Arguments:
	filename: name of the fasta file
Returns:
	s: string with relevant sequence
'''
def read_fasta(filename):
    with open(filename, "r") as f:
        s = ""
        for l in f.readlines()[1:]:
            s += l.strip()
    return s


"""
Calculate the average gene length.

Argument:
    genes: List of tuples (start_inc, end_inc) of gene locations.

Returns:
    The average gene length
"""
def get_mean_genes(genes):
    length = 0
    for (start_inc, end_inc) in genes:
        length += end_inc - start_inc + 1

    return length / len(genes)

"""
Reads the genes file and outputs the gene locations to analyze.

Arguments:
    filename: name of the gene info file

Returns:
    List of tuples (start_inc, end_inc) of gene locations.
"""
def read_genes(filename, complement=False):
    genes = []
    prev_end = -1
    with open(filename, "r") as f:
        for l in f.readlines():
            if " gene " in l and (complement != ("complement" not in l)):
                l = l.replace('complement(', '').replace('join(', '').replace(')', '')
                location = l.split()[1]
                location = location.split(',')
                for pair in location:
                    pair = pair.split("..")
                    if (int(pair[0])-2) <= prev_end:
                        continue

                    genes.append((int(pair[0])-1, int(pair[1])-1))
                    prev_end = int(pair[1])-1

    return genes

"""
Gets the complement of a sequence from the original.

Arguments:
    obs_sequence: sequence

Returns:
    String of complement sequence
"""
def get_complement_sequence(obs_sequence):
    complement = []
    for obs_i in obs_sequence:
        if obs_i=='A':
            complement.append('T')
        elif obs_i=='T':
            complement.append('A')
        elif obs_i=='G':
            complement.append('C')
        elif obs_i=='C':
            complement.append('G')

    return ''.join(complement)


"""
Given a list of state names and states, outputs gene locations.

Arguments:
    state_sequence: list of state names
    states: mapping of state names to state objects

Returns:
    List of tuples (start_inc, end_inc) of gene locations.
"""
def output_genes(state_sequence, states):
    genes = []
    gene_start = None

    for i, state_name in zip(range(len(state_sequence)), state_sequence):
        state = states[state_name]

        if gene_start is None:
            if state.name == START_CODON + "0":
                gene_start = i
        elif state.name == STOP_CODON_A + "2" or state.name == STOP_CODON_B + "2":
            gene_end = i
            genes.append((gene_start, gene_end))
            gene_start = None

    if gene_start:
        genes.append((gene_start, len(state_sequence)-1))
    return genes


"""
Compare a gene output with the expected.
Precondition: expected and genes_tested are sorted by start_inc

Arguments:
    expected: List of tuples of genes (start_inc, end_inc)
    genes_tested: List of tuples of genes (start_inc, end_inc)

Returns a tuple of:
    tp, fp, fn
"""
def get_metrics(expected, genes_tested):
    true_positive = 0
    false_positive = 0
    false_negative = 0
    exact = 0

    expected_index = 0
    genes_index = 0

    # definition of true positive: any overlapping guesses
    # definition of false positive: predicts a gene that does not overlap
    # definition of false negative: real gene does not overlap with predictions

    while expected_index < len(expected) and genes_index < len(genes_tested):
        (expected_start, expected_end) = expected[expected_index]
        (gene_start, gene_end) = genes_tested[genes_index]
        
        if overlaps(expected[expected_index], genes_tested[genes_index]):
            expected_index += 1
            genes_index += 1
            true_positive += 1
            if expected_start == gene_start or expected_end == gene_end:
                exact += 1
            continue

        if expected_start < gene_start:
            false_negative += 1
            expected_index += 1
        elif expected_start > gene_start:
            false_positive += 1
            genes_index += 1

    false_positive += len(genes_tested) - genes_index
    false_negative += len(expected) - expected_index

    return (true_positive, false_positive, false_negative, exact)


"""
Determines if two genes overlap.

Arguments:
    gene_1: (start_inc, end_inc) of gene 1
    gene_2: (start_inc, end_inc) of gene 2

Returns:
    True or False
"""
def overlaps(gene_1, gene_2):
    # start of gene 1 falls between gene 2
    if gene_1[0] >= gene_2[0] and gene_1[0] <= gene_2[1]:
        return True

    # start of gene 2 falls between gene 1
    if gene_2[0] >= gene_1[0] and gene_2[0] <= gene_1[1]:
        return True

    return False
