import kron as kron_ops
import numpy as np
import expansion
import copy
import matplotlib.pyplot as plt


def expanded_cnot2():
    """ 2 vars only """
    return np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])


def swap():
    """ 2 vars only """
    return np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])



def expanded_cnot(start, end, length):
    """
    start, end, length use ordinary python count from 0
    Extremely slow but will be precached in timing test
    """
    # may need a mapping from possible inputs to index in mat, and back
    # this mapping can be calculated, not just a dict

    # go through the possible inputs
    # put a 1 on the possible outputs
    res = np.zeros((2 ** length, 2 ** length))
    for possible_interp in expansion.possible_interps(length):
        corresponding_interp = list(copy.copy(possible_interp))
        if possible_interp[start] == 1:
            corresponding_interp[end] = 1 - corresponding_interp[end]
        res[
            expansion.interp_to_idx(possible_interp),
            expansion.interp_to_idx(corresponding_interp),
        ] = 1
    return res


def expanded_toffoli(fst, snd, end, length):
    res = np.zeros((2 ** length, 2 ** length))
    for possible_interp in expansion.possible_interps(length):
        corresponding_interp = list(copy.copy(possible_interp))
        corresponding_interp[end] ^= (possible_interp[fst] * possible_interp[snd])
        res[
            expansion.interp_to_idx(possible_interp),
            expansion.interp_to_idx(corresponding_interp),
        ] = 1
    return res


def test_expanded_cnot():
    print(expanded_cnot(1, 2, 3))


def test_expanded_toffoli():
    plt.imshow(expanded_toffoli(0,1,2,3))
    plt.show()

if __name__ == "__main__":
    test_expanded_toffoli()
