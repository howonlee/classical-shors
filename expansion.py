import numpy as np
import itertools
import string


def expand_res_by_function(fn, num_vars):
    expanded_res = np.zeros(2 ** num_vars)
    for possible_interp in possible_interps(num_vars):
        assignments = ls_to_dict(possible_interp)
        res = fn(assignments)
        expanded_res[interp_to_idx(possible_interp)] = res
    return expanded_res


def ls_to_dict(ls):
    assert len(ls) < len(string.ascii_uppercase)
    return dict(zip(string.ascii_uppercase[: len(ls)], ls))


def interp_to_idx(interp):
    return sum([b << i for i, b in enumerate(reversed(interp))])

def idx_to_interp(idx, num_vars):
    return tuple(int(x) for x in format(idx, 'b').zfill(num_vars))


def possible_interps(num_vars):
    return itertools.product(*[[0, 1] for _ in range(num_vars)])

def test_interp_to_idx():
    for orig_idx in range(200):
        interp = idx_to_interp(orig_idx, 8)
        new_idx = interp_to_idx(interp)
        assert new_idx == orig_idx

if __name__ == "__main__":
    print(expand_res_by_function(lambda x: 1, 6))
