import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from rnn_utils import prepare_sequence
import numpy as np

torch.manual_seed(1)

# From PyTorch tutorials
class LSTMTagger(nn.Module):

    def __init__(self, embedding_dim, hidden_dim, tagset_size):
        super(LSTMTagger, self).__init__()
        self.hidden_dim = hidden_dim

        # The LSTM takes word embeddings as inputs, and outputs hidden states
        # with dimensionality hidden_dim.
        self.lstm = nn.LSTM(input_size=embedding_dim, hidden_size=hidden_dim)

        # The linear layer that maps from hidden state space to tag space
        self.hidden2tag = nn.Linear(hidden_dim, tagset_size)
        self.out = nn.Sigmoid()

    def forward(self, sentence):
        lstm_out, lstm_hidden = self.lstm(sentence)
        tag_space = self.hidden2tag(lstm_hidden[0].view(len(sentence[0]), -1))
        tag_scores = self.out(tag_space)
        return tag_scores

"""
Trains LSTM on the labels.

Arguments:
    training_labels: List of list of tuples (DNA_base, 0 or 1)
    test_labels: List of list of tuples (DNA_base, 0 or 1)
    epochs: number of epochs

Returns:

"""
def run_rnn(training_labels, test_labels, epochs):
    tag_embedding = {'A': [0,0,0,1], 'T': [0,0,1,0], 'G': [0,1,0,0], 'C': [1,0,0,0]}
    model = LSTMTagger(4, 10, 1)
    loss_function = nn.BCELoss()
    optimizer = optim.SGD(model.parameters(), lr=0.1)

    torch.set_printoptions(threshold=5000)

    #training_labels = [[('A',0), ('T',1), ('C',0)]]

    for e in range(epochs):
        print("Epochs", e)
        for training_label, i in zip(training_labels, range(len(training_labels))):
            if i%100==0:
                print("Training label", i)
            model.zero_grad()
            optimizer.zero_grad()
            inputs, tags = prepare_sequence(training_label, tag_embedding)
            tag_scores = model(torch.tensor(np.array([inputs], dtype=float), requires_grad=True).float()) #batch of 1
            loss = loss_function(tag_scores, torch.FloatTensor(tags).reshape((-1,1)))
            loss.backward()
            optimizer.step()
            
    with torch.no_grad():
        predicted_gene = []
        for test_label in test_labels:
            inputs, tags = prepare_sequence(test_label, tag_embedding)
            tag_scores = model(torch.tensor(np.array([inputs], dtype=float), requires_grad=True).float())#batch of 1
            #tag_scores = [0 if i < 0.21 else 1 for i in tag_scores]
            predicted_gene.extend(tag_scores.reshape(1,-1).tolist()[0])

    predicted_gene_average = sum(predicted_gene) / len(predicted_gene)
    predicted_gene = [0 if i < predicted_gene_average else 1 for i in predicted_gene]

    print(predicted_gene)
    return predicted_gene

"""
split = [[('A', 0), ('T', 1), ('C', 1)], [('G', 0), ('T', 1), ('A', 0)]]
run_rnn(split, split, 1000)
"""
