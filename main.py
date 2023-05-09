from Bio import Phylo
import csv

import utils

PATH_FASTA = "./msa.fasta"
PATH_TREE  = "./tree.tre"
PATH_ANCES = "./ancestrals.csv"
SEQ_LEN = 96




def load_tree(path):
    f = open(path, 'r')
    tree = Phylo.read(path, "newick")
    tree.ladderize()
    return tree


def load_fasta(path):

    f = open(path, 'r')
    reads = f.readlines()

    msa = {}
    name = ""
    dat = []

    for line in reads:
        line = line.strip('\n')
        # if first read
        if line[0] == '>' and name == "":
            name = line[1:]
        # If new read store the old and start new
        elif line[0] == '>' and name != "":
            msa[name] = dat[0] + dat[1]
            dat = []
            name = line[1:]
        else:
            dat.append(line)
    # Append last read
    msa[name] = dat[0]+dat[1]
    if len(dat[0]+dat[1]) < 96:
        print("Smaller")
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
    return seq



def travel_tree(tree, msa, anc, seq):

    met = []
    if len(tree.clades) <= 0:
        met.append({"name": tree.name, "conf": tree.confidence, "bl": tree.branch_length})
        return met
    # Left travel
    ls = travel_tree(tree.clades[0], msa, anc, seq)
    # Right travel
    rs = travel_tree(tree.clades[1], msa, anc, seq)

    lsrs = ls + rs
    res_seq = get_resulting_seq(tree.confidence, anc)


    for i in range(0, len(res_seq)):
        amin_sum, gap_sum = [0, 0]
        for rec in lsrs:
            amin = msa[rec["name"]][i]
            if amin == "-":
                gap_sum += rec["bl"]
            else:
                amin_sum += rec["bl"]
        if gap_sum > amin_sum:
            res_seq[i] = '-'

    utils.store_res(res_seq, str(tree.confidence))
    # Adjust weights by adding new branch len
    for rec in lsrs:
        rec["bl"] += tree.branch_length
    return lsrs

def fil_tree(tree, msa, anc):

    r_node = tree.root
    # Left travel
    ls = travel_tree(tree.clade.clades[0], msa, anc, [])
    # Right travel
    rs = travel_tree(tree.clade.clades[1], msa, anc, [])

    lsrs = ls + rs
    res_seq = get_resulting_seq(r_node.confidence, anc)


    for i in range(0, len(res_seq)):
        amin_sum, gap_sum = [0, 0]
        for rec in lsrs:
            amin = msa[rec["name"]][i]
            if amin == "-":
                gap_sum += rec["bl"]
            else:
                amin_sum += rec["bl"]
        if gap_sum > amin_sum:
            res_seq[i] = '-'

    utils.store_res(res_seq, str(r_node.confidence))


if __name__ == '__main__':
    t = load_tree(PATH_TREE)
    m = load_fasta(PATH_FASTA)
    a = load_ances(PATH_ANCES)
    fil_tree(t, m, a)