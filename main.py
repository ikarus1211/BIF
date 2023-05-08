from Bio import Phylo
import matplotlib
import matplotlib.pyplot as plt
import pylab
import csv


PATH_FASTA = "./msa.fasta"
PATH_TREE  = "./tree.tre"
PATH_ANCES = "./ancestrals.csv"
SEQ_LEN = 96



def show_tree(tree, label=False):
    matplotlib.rc('font', size=10)
    # set the size of the figure
    fig = plt.figure(figsize=(25, 25), dpi=100)
    # alternatively
    # fig.set_size_inches(10, 20)
    axes = fig.add_subplot(1, 1, 1)
    if label:
        Phylo.draw(tree,branch_labels=lambda c: c.branch_length, axes=axes)
    else:
        Phylo.draw(tree, axes=axes)


def load_tree(path):
    f = open(path, 'r')
    tree = Phylo.read(path, "newick")
    tree.ladderize()

    return tree


def load_fasta(path):

    f = open(path, 'r')
    reads = f.readlines()

    msa = []
    name = ""
    dat = []

    for line in reads:
        line = line.strip('\n')
        # if first read
        if line[0] == '>' and name == "":
            name = line[1:]
        # If new read store the old and start new
        elif line[0] == '>' and name != "":
            msa.append({"name": name, "reads": dat[0]+dat[1]})
            dat = []
            name = line[1:]
        else:
            dat.append(line)
    # Append last read
    msa.append({"name": name, "reads": dat[0]+dat[1]})
    return msa


def load_ances(path):
    msa = {}
    with open(path, 'r') as csvfile:
        reads = list(csv.reader(csvfile, delimiter=','))
        global AMIN_POS
        AMIN_POS = reads[0][2:]
        for read in reads[1:]:

            if read[0] in msa:

                msa[read[0]].append(read[2:])
            else:

                msa[read[0]] = [read[2:]]
    return msa

def find_max_prob(prob):
    maximum = 0.0
    max_indx = 0
    for i,p in enumerate(prob):
        if p == '-':
            continue
        else:
            if float(p) > maximum:
                maximum = float(p)
                max_indx = i
    return maximum, max_indx

def get_resulting_seq(current, anc):
    seq = []
    probs = anc[str(current)]
    for prob in probs:
        maxprob, idx = find_max_prob(prob)
        seq.append(AMIN_POS[idx])
    print("done")



def travel_tree(tree, msa, anc, seq):

    if len(tree.clades) <= 0:
        meta.append({"name": tree.name, "conf": tree.confidence, "bl": tree.branch_length})
        return meta, []
    # Left travel
    ls = travel_tree(tree.clades[0], msa, anc, seq)
    # Right travel
    rs = travel_tree(tree.clades[1], msa, anc, seq)

    res_seq = get_resulting_seq(tree.confidence, anc)
    return

def fil_tree(tree, msa, anc):

    r_node = tree.root


    # Left travel
    travel_tree(tree.clade.clades[0], msa, anc, [])
    # Right travel
    travel_tree(tree.clade.clades[1], msa, anc, [])
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    t = load_tree(PATH_TREE)
    m = load_fasta(PATH_FASTA)
    a = load_ances(PATH_ANCES)
    fil_tree(t,m,a)