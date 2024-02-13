# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    seqs_dict = {'Gallus gallus': (gg_seq, gg_header),
                 'Mus musculus': (mm_seq, mm_header),
                 'Balaeniceps rex': (br_seq, br_header),
                 'Tursiops truncatus': (tt_seq, tt_header),
                 }

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    alignment_dict = {}
    for species, (seq, _) in seqs_dict.items():
        nw = NeedlemanWunsch('substitution_matrices/BLOSUM62.mat', -10, -1)
        alignment_dict[species] = nw.align(hs_seq, seq)
    print(list(sorted(alignment_dict.keys(), key=lambda x: alignment_dict[x][0])))

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    print('Scores:')
    for species, score in alignment_dict.items():
        print(species, score[0])
    

if __name__ == "__main__":
    main()
