import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from rnn_utils import prepare_sequence

torch.manual_seed(1)

# From PyTorch tutorials
class LSTMTagger(nn.Module):

    def __init__(self, embedding_dim, hidden_dim, tagset_size):
        super(LSTMTagger, self).__init__()
        self.hidden_dim = hidden_dim

        # The LSTM takes word embeddings as inputs, and outputs hidden states
        # with dimensionality hidden_dim.
        self.lstm = nn.LSTM(embedding_dim, hidden_dim)

        # The linear layer that maps from hidden state space to tag space
        self.hidden2tag = nn.Linear(hidden_dim, tagset_size)

    def forward(self, sentence):
        lstm_out, _ = self.lstm(sentence)
        tag_space = self.hidden2tag(lstm_out.view(len(sentence[0]), -1))
        tag_scores = F.log_softmax(tag_space, dim=1)
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
    tag_embedding = {'A': [0,0], 'T': [1,1], 'G': [1,0], 'C': [1,1]}
    model = LSTMTagger(2, 5, 2)
    loss_function = nn.NLLLoss()
    optimizer = optim.SGD(model.parameters(), lr=0.1)

    for _ in range(epochs):
        for training_label in training_labels:
            inputs, tags = prepare_sequence(training_label, tag_embedding)
            tag_scores = model(torch.FloatTensor([inputs])) #batch of 1
            loss = loss_function(tag_scores, torch.tensor(tags))
            loss.backward()
            optimizer.step()

    with torch.no_grad():
        predicted_gene = []
        for test_label in test_labels:
            inputs, tags = prepare_sequence(test_label, tag_embedding)
            tag_scores = model(torch.FloatTensor([inputs])) #batch of 1
            values, predicted = torch.max(tag_scores, 1)
            predicted_gene.extend(predicted.tolist())

    print(predicted_gene)
    return predicted_gene

split = [[('A', 0), ('T', 1), ('C', 1)], [('G', 0), ('T', 1), ('A', 0)]]
run_rnn(split, split, 10)
