#!/usr/bin/env python3

import csv
from enum import Enum
from math import log, exp

INTERGENE = "INT"
STOP_CODON_A = "SPA"
STOP_CODON_B = "SPB"
START_CODON = "STT"


"""
Reads absolute counts of a codon from a given file.

Arguments:
    file_name: path to a csv file containing frequencies.

Returns:
    Dictionary mapping codon name to its values
    List of column names
    Total number of codons
"""
def get_codon_frequencies(file_name):
    frequencies = dict()
    with open(file_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        column_names = []
        for row in csv_reader:
            if line_count==0:
                column_names = row
            else:
                codon_name = row[0]
                codon_information = dict()
                codon_information[column_names[1]] = row[1]
                codon_information[column_names[2]] = float(row[2])
                codon_information[column_names[3]] = float(row[3])
                frequencies[codon_name] = codon_information

            line_count += 1

    return frequencies, column_names, line_count-1

"""
Defines an enum to classify types of states.
"""
class StateType(Enum):
    DEFAULT = 0
    AMINO_ACID = 1
    STOP_CODON = 2
    START_CODON = 3
    INTERGENE = 4

"""
Defines a state containing transition and emission probabilities.
"""
class State:

    """
    Constructor.
    Arguments:
        meta_name: "Group" name, e.g codon or type of stop state
        state_type: StateType enum
        state_num: 0-2 if not an intergene state

    Properties (examples)
        meta_name: "ATT"
        name: "ATT0"
        state_type: StateType.AMINO_ACID
        self.state_num: 0
        self.transitions: {"ATT1": log(1)}
        self.emissions: {"A": log(1)}
    """
    def __init__(self, meta_name, state_type, state_num):
        self.meta_name = meta_name
        self.name = None
        self.state_type = state_type
        self.state_num = int(state_num)
        self.transitions = dict()
        self.emissions = dict()

    """
    Initializes transition probabilities.

    Arguments:
        central_prob: Dict of transitions through central state
    """
    def init_transitions(self, central_prob):
        # leads into central state
        if self.state_type == StateType.AMINO_ACID or \
                self.state_type == StateType.START_CODON:
            self.name = self.meta_name + str(self.state_num)
            if self.state_num != 2:
                next_state = self.meta_name + str(self.state_num+1)
                self.transitions[next_state] = log(1.0)
            else:
                self.transitions = central_prob
        elif self.state_type == StateType.STOP_CODON:
            self.name = self.meta_name + str(self.state_num)
            if self.state_num != 2:
                next_state = self.meta_name + str(self.state_num+1)
                self.transitions[next_state] = log(1.0)
            else:
                self.transitions[INTERGENE] = log(1.0)
        elif self.state_type == StateType.INTERGENE:
            self.name = self.meta_name
            self.transitions[INTERGENE] = log(0.5) #some prob of staying
            self.transitions[START_CODON] = log(0.5)

    """
    Initializes emission probabilities

    Arguments:
        start_prob: Emission probabilities of the start states
        inter_prob: Emission probabilities of the intergene state
        stop_prob: Emission probabilities of the stop states
    """
    def init_emissions(self, start_prob, inter_prob, stop_prob):
        if self.state_type == StateType.AMINO_ACID:
            base = self.meta_name[self.state_num]
            self.emissions[base] = log(1.0)
        elif self.state_type == StateType.STOP_CODON:
            self.emissions = stop_prob[self.meta_name][self.state_num]
        elif self.state_type == StateType.INTERGENE:
            self.emissions = inter_prob
        elif self.state_type == StateType.START_CODON:
            self.emissions = start_prob[self.state_num]


def get_initial_probabilities(codon_frequency, total_freq):
    init_prob = dict()

    for codon in codon_frequency:
        if codon_frequency[codon]["Aa"] != "*":
            init_prob[codon + "0"] = log(codon_frequency[codon]["Usage"]) - log(total_freq)

    init_prob[START_CODON + "0"] = log(0.1) # to be replaced
    init_prob[INTERGENE] = log(0.1)

    # assume you cannot start in the middle of a codon or at a stop codon

    return init_prob

"""
Returns emission probabilities for a few types of states.

Returns:
    List emission probabilities (as dicts) for the start states
    Dict emission probabilities (as dicts) for the stop states
    Emission probabilities (as dicts) for the intergene state
"""
def get_emission_probabilities():
    start_prob = [{'A': log(0.25), 'C': log(0.25), 'G': log(0.25), 'T': log(0.25)},
            {'A': log(0.25), 'C': log(0.25), 'G': log(0.25), 'T': log(0.25)},
            {'A': log(0.25), 'C': log(0.25), 'G': log(0.25), 'T': log(0.25)}]
    stop_prob = {
            STOP_CODON_A: [{'A': log(0.25), 'C': log(0.25), 'G': log(0.25), 'T': log(0.25)},
                {'A': log(0.25), 'C': log(0.25), 'G': log(0.25), 'T': log(0.25)},
                {'A': log(0.25), 'C': log(0.25), 'G': log(0.25), 'T': log(0.25)}],
            STOP_CODON_B: [{'A': log(0.25), 'C': log(0.25), 'G': log(0.25), 'T': log(0.25)},
                {'A': log(0.25), 'C': log(0.25), 'G': log(0.25), 'T': log(0.25)},
                {'A': log(0.25), 'C': log(0.25), 'G': log(0.25), 'T': log(0.25)}]
            }
    inter_prob = {'A': log(0.25), 'C': log(0.25), 'G': log(0.25), 'T': log(0.25)}
    return start_prob, inter_prob, stop_prob


"""
Returns transition probabilities for any state going through the central state.

Arguments:
    codon_frequency: A mapping of codon name to dictionary of values
    total_freq: An integer value of the total observation counts

Returns:
    Dictionary of end_state : log probability
"""
def get_central_probabilities(codon_frequency, total_freq):
    central_prob = dict()

    for codon in codon_frequency:
        if codon_frequency[codon]["Aa"] != "*":
            central_prob[codon + "0"] = log(codon_frequency[codon]["Usage"]) - log(total_freq)

    central_prob[STOP_CODON_A + "0"] = log(0.1) # to be replaced
    central_prob[STOP_CODON_B + "0"] = log(0.1) # to be replaced

    return central_prob


"""
Returns a dictionary mapping state name to State object

Arguments:
    codon_frequency: A mapping of codon name to dictionary of values
    start_prob: Emissions of the start codon states
    inter_prob: Emissions of the intergene state
    Stop_prob: Dict mapping stop codon type to their emissions

Returns:
    Dictionary mapping state name to State object (initialized)
"""
def get_states(codon_frequency, start_prob, inter_prob, stop_prob):
    states = dict()
    for codon in codon_frequency:
        state_type = None
        # what should we do about start codons?
        if codon_frequency[codon]["Aa"] != "*":
            state_type = StateType.AMINO_ACID
            for i in range(3):
                state = State(codon, state_type, i)
                state.init_transitions(codon_frequency)
                state.init_emissions(start_prob, inter_prob, stop_prob)
                states[codon + str(i)] = state

    for i in range(3):
        state = State(STOP_CODON_A, StateType.STOP_CODON, i)
        state.init_transitions(codon_frequency)
        state.init_emissions(start_prob, inter_prob, stop_prob)
        states[STOP_CODON_A + str(i)] = state

        state = State(STOP_CODON_B, StateType.STOP_CODON, i)
        state.init_transitions(codon_frequency)
        state.init_emissions(start_prob, inter_prob, stop_prob)
        states[STOP_CODON_B + str(i)] = state
 
        state = State(START_CODON, StateType.START_CODON, i)
        state.init_transitions(codon_frequency)
        state.init_emissions(start_prob, inter_prob, stop_prob)
        states[START_CODON + str(i)] = state

    state = State(INTERGENE, StateType.INTERGENE, -1)
    state.init_transitions(codon_frequency)
    state.init_emissions(start_prob, inter_prob, stop_prob)
    states[INTERGENE] = state

    return states
