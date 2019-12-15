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
Reads the genes file and outputs the gene locations to analyze.

Arguments:
    filename: name of the gene info file

Returns:
    List of tuples (start_inc, end_inc) of gene locations.
"""
def read_genes(filename, complement=False):
    genes = []
    with open(filename, "r") as f:
        for l in f.readlines():
            if not complement:
                if " gene " in l and "complement" not in l:
                    location = l.split()[1].split("..")
                    genes.append((int(location[0]),int(location[1])))
            else:
                if " gene " in l and "complement" in l:
                    location = l.split()[1]
                    location = location.replace('complement(', '').replace(')', '')
                    location = location.split("..")
                    genes.append((int(location[0]),int(location[1])))
    
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
            if state.name == START_CODON + "2":
                gene_start = i+1
        elif state.name == STOP_CODON_A + "0" or state.name == STOP_CODON_B + "0":
            gene_end = i-1
            genes.append((gene_start, gene_end))
            gene_start = None

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

    expected_index = 0
    genes_index = 0

    while expected_index < len(expected) and genes_index < len(genes_tested):
        (expected_start, expected_end) = expected[expected_index]
        (gene_start, gene_end) = genes_tested[genes_index]
        
        if expected_start == gene_start:
            if expected_end == gene_end:
                true_positive += 1
            else:
                false_positive += 1
                false_negative += 1
            expected_index += 1
            genes_index += 1
        elif expected_start < gene_start:
            false_negative += 1
            expected_index += 1
        elif expected_start > gene_start:
            false_positive += 1
            genes_index += 1

    false_positive += len(genes_tested) - genes_index
    false_negative += len(expected) - expected_index

    return (true_positive, false_positive, false_negative)

