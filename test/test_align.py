# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    nw = NeedlemanWunsch('substitution_matrices/BLOSUM62.mat', -10, -1)

    # test outputs
    assert (4.0, 'MYQR', 'M-QR') == nw.align(seq1, seq2)
    
    # test score matrix
    assert nw._F == np.darray([[0, -10, -11, -12,],[-10,   5,  -6,  -7,],[-11,  -6,   4,  -7,],[-12,  -7,  -1,   5,],[-13,  -8,  -6,   4,]])

    # test backtrace matrix
    assert nw._backmat == np.ndarray([[0, 2, 2, 2],[1, 0, 1, 1],[1, 2, 0, 1],[1, 2, 0, 0],[1, 2, 0, 0]])
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    nw = NeedlemanWunsch('substitution_matrices/BLOSUM62.mat', -10, -1)
    assert (17.0, 'MAVHQLIRRP', 'M---QLIRHP') == nw.align(seq3, seq4)




