import matplotlib
import matplotlib.pyplot as plt
from Bio import Phylo


def show_tree(tree, label=False):
    matplotlib.rc('font', size=10)
    # set the size of the figure
    fig = plt.figure(figsize=(25, 25), dpi=100)

    axes = fig.add_subplot(1, 1, 1)
    if label:
        Phylo.draw(tree,branch_labels=lambda c: c.branch_length, axes=axes)
    else:
        Phylo.draw(tree, axes=axes)


def store_res(seq, filename):

    seq = ''.join(seq)

    f = open("./results/node_" + filename + ".fas", 'w')
    f.write(">" + filename + "\n")
    f.write(seq)
    f.close()