#!/usr/bin/env python3

import numpy as np


''' Outputs the Viterbi decoding of a given observation.
Arguments:
	obs: observed sequence of emitted states (list of emissions)
	init_probs: initial log-probabilities for each hidden state (dictionary)
Returns:
	l: list of most likely hidden states at each position (list of hidden
           states)
	p: log-probability of the returned hidden state sequence
'''
def viterbi(obs, init_probs, states):
    ''' Complete this function. '''
    K = len(init_probs) # states
    N = len(obs) # emissions
    V = dict() # viterbi matrix, V[state] returns a list
    B = dict() # backpointers, l[state] returns a list

    # V and B may contain "None" values for impossible paths

    for state_name in states:
        state = states[state_name]
        V[state_name] = []
        B[state_name] = []
        B[state_name].append("")
        value = None
        if state_name in init_probs and obs[0] in state.emissions:
            # if possible to start at this state, and emit the first base
            value = init_probs[state_name] + state.emissions[obs[0]]
        V[state_name].append(value)

    for i in range(1, len(obs)):
        for state_name in states:
            state = states[state_name]
            (prev_max, max_state_name) = get_max_prev_state(V, state_name, i, states)
            value = None
            if prev_max is not None and obs[i] in state.emissions:
                # if possible to come to this state, and emit this base
                value = state.emissions[obs[i]] + prev_max

            V[state_name].append(value)
            B[state_name].append(max_state_name)

    final_state_name = None
    max_final = None
    for state_name in states:
        value = V[state_name][len(obs)-1]
        if value is None:
            # not possible to end at this state
            continue

        if max_final is None or value > max_final:
            max_final = value
            final_state_name = state_name

    print("If this assert fails, no viable path.")
    assert(final_state_name is not None)
    print("Assert passed.")

    l = []
    state_name = final_state_name
    for i in range(len(obs)-1, -1, -1):
        assert(state_name is not None)
        l.append(state_name)
        state_name = B[state_name][i]

    l.reverse()
    print("Final backpointer trace: ", l)
    return (l, max_final)


def get_max_prev_state(V, state_name, i, states):
    prev_max = None
    max_state_name = None

    for prev_state_name in states:
        prev_V = V[prev_state_name][i-1]
        prev_state = states[prev_state_name]
        if prev_V is None or state_name not in prev_state.transitions:
            continue

        # if possible to come from prev_state_name and transit here
        trans_value = prev_state.transitions[state_name]
        value = prev_V + trans_value

        if prev_max is None or value > prev_max:
            prev_max = value
            max_state_name = prev_state_name

    # print(prev_max, max_state_name)
    return (prev_max, max_state_name)
