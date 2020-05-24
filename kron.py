import numpy as np
import scipy.linalg as sci_lin
import matplotlib.pyplot as plt
import random
import functools
import collections
import copy

"""
General utils for dealing with prepositional statements with linear algebra
Also gets used by the other files.
Also has implementation of the linear operator application to mix of
Kronecker product and Kronecker sum, which is how we get to deal with entangling in a fast manner
"""

"""
All this stuff is bras instead of kets,
because I don't think about the outer product operator making at all here yet,
and row vectors are easier to deal with in numpy.
The exposition is in kets, though. Just transpose stuff in your head
"""

"""
Operations are encoded with lists of symbols, * for Kronecker product, + for Kronecker sum.
This is a problem because Kronecker product and Kronecker sum are associative,
but not with respect to each other (think of ordinary product and sum).

But making an actual tree is straightforward, it's just that I'm not doing it.
Consider a list [* + + *] to encode ((((A * B) + C) + D) * E), all Kronecker operators of course
"""

IN_GF2 = False


def gf2(arr):
    """
    Coerces arr to say that we are in Galois field of 2 elements.
    Assumes integer valued inputs
    """
    if IN_GF2:
        return np.mod(arr, 2)
    else:
        return arr


def kron_and(pair):
    return np.kron(pair[0], pair[1])


def kron_xor(pair):
    return gf2(
        np.kron(pair[0], np.ones_like(pair[1]))
        + np.kron(np.ones_like(pair[0]), pair[1])
    )


def kron_or(pair):
    xor_res = kron_xor(pair)
    and_res = kron_and(pair)
    return gf2(xor_res + and_res)


def apply_op(op, pair):
    ops_to_funcs = {
        "*": kron_and,
        "+": kron_xor,
        "v": kron_or,
    }
    return ops_to_funcs[op](pair)


def expand(vecs, ops):
    vecs, ops = vecs[:], ops[:]
    curr_pair = (vecs[0], vecs[1])
    vecs.pop(0)
    while vecs:
        vecs.pop(0)
        curr_op = ops.pop(0)
        curr_res = apply_op(curr_op, curr_pair)
        if vecs:
            curr_pair = (curr_res, vecs[0])
    return curr_res


def sierpinski(order):
    """
    Extraneous but fun!
    Changes bases to Zhegalkin basis, iirc.
    Factorizing to linear operator composed of
    bunch of Kronecker products is trivial
    """
    root = np.array([[1, 0], [1, 1]])
    curr_res = root.copy()
    for curr_order in range(order - 1):
        curr_res = np.kron(curr_res, root)
    return curr_res


def rand_vec():
    res = random.choice(
        [np.array([1, -1]), np.array([1, 0]), np.array([0, 1]), np.array([1, 1])]
    )
    return res.copy()


def rand_op():
    return random.choice(["*", "+"])


def rand_proplin(length):
    assert length >= 2, "Prop lin length >= 2 for now"
    vecs = [rand_vec() for _ in range(length)]
    ops = [rand_op() for _ in range(length - 1)]
    return (vecs, ops)


def test_entangled_property():
    assert not IN_GF2, "Needs to be not in GF2, futz with the constant"
    hads = [np.array([[0.2, 0.8], [1.2, 2.2]]) for _ in range(3)]
    vecs = list(map(np.array, [[1, 0], [1, 0], [1, 1]]))
    ops = ["+", "+"]
    expand_res = expand(vecs, ops)
    mult = functools.reduce(np.kron, hads)
    print(mult)
    print(expand_res)
    print(mult.dot(expand_res))

    fst = [hads[0].dot(vecs[0]), hads[1].dot([1, 1]), hads[2].dot([1, 1])]
    snd = [hads[0].dot([1, 1]), hads[1].dot(vecs[1]), hads[2].dot([1, 1])]
    thd = [hads[0].dot([1, 1]), hads[1].dot([1, 1]), hads[2].dot(vecs[2])]
    allprod = ["*", "*"]
    print(expand(fst, allprod) + expand(snd, allprod) + expand(thd, allprod))


def apply_linear_op_product(linear_op_mats, vecs):
    return [linear_op_mats[i].dot(vec) for i, vec in enumerate(vecs)]


def apply_linear_op_entangled(linear_op_mats, vecs, ops):
    """
    Applies it our quick way
    Ignores normalization. I am pretending I am a theorist, beware!

    Can't return one proposition because entangled, returns a list of propositions.
    So understand that that list of propositions consists of an entangled state.
    But by construction it's a list of product propositions, at least.
    """
    num_summations = collections.Counter(ops)["+"]
    our_res = [[None for _ in range(len(vecs))] for _ in range(num_summations + 1)]
    """
    Reverse these because of the associative assumption
    If we consider a list [* * * *] to encode ((((A * B) * C) * D) * E), what we want is,
    (E * (D * (C * (B * A))))

    Obviously if you're shipping something you would need to actually deal with a tree,
    but that is straightforward (if more naturally recursive)
    """
    vecs = list(reversed(copy.deepcopy(vecs)))
    linear_mats = copy.deepcopy(linear_op_mats)
    ops = list(reversed(copy.deepcopy(ops)))
    curr_summation_idx = 0

    for vec_idx, curr_vec in enumerate(vecs):
        if vec_idx == len(vecs) - 1:
            curr_op = ops[-1]
        else:
            curr_op = ops[vec_idx]

        if curr_op == "*":
            for res_idx, res in enumerate(our_res):
                # If it's already "colored in" don't mess with it!
                if res[vec_idx] is None:
                    res[vec_idx] = vecs[vec_idx].dot(linear_mats[vec_idx])
        elif curr_op == "+":
            for res_idx, res in enumerate(our_res):
                """
                We are sort of "coloring in" the remaining portion of the row of the table
                w/ the results of linear mat w/ vec idx, with the linear mat w/ superposition
                """
                if res_idx == curr_summation_idx:
                    for curr_vec_idx in range(vec_idx, len(vecs)):
                        res[curr_vec_idx] = np.array([1, 1]).dot(
                            linear_mats[curr_vec_idx]
                        )
                    res[vec_idx] = vecs[vec_idx].dot(linear_mats[vec_idx])

                else:
                    # Coloring in the column here
                    res[vec_idx] = np.array([1, 1]).dot(linear_mats[vec_idx])
            curr_summation_idx += 1
        else:
            raise Exception("wrong op")

    # Undo reversion
    return list(map(lambda x: list(reversed(x)), our_res))


def test_entangled_linear_op():
    assert not IN_GF2, "Needs to be not in GF2, futz with the constant"
    prop_len = 3
    linear_op = [np.array([[1, 1], [0, 1]]) for _ in range(prop_len)]
    vecs, ops = rand_proplin(prop_len)
    expand_res = expand(vecs, ops)
    mult = functools.reduce(np.kron, linear_op)

    our_res = apply_linear_op_entangled(linear_op, vecs, ops)
    allprod = ["*" for _ in range(prop_len - 1)]
    total_res = sum(map(lambda x: expand(x, allprod), our_res))
    if not np.allclose(expand_res.dot(mult), total_res):
        print("Failed entangled linear operator")
        print("correct res: ", mult.dot(expand_res))
        print("total our res: ", total_res)
        print("our res, piece by piece: ")
        for res in our_res:
            print(res, expand(res, allprod))
        print(vecs)
        print(ops)
        raise Exception()


def repeated_test_entangled_linear_op():
    print(
        """
    Trying 1000 random examples of entangled state to deal with in linear op,
    comparing our answer to fully expanded ansewr done simply with matrix multiplication.

    If you were shipping in production, you would never do this.

    Will bomb out if not correct.
    """
    )
    for idx in range(1000):
        if idx % 30 == 0:
            print(idx)
        test_entangled_linear_op()
    print("Done!")


if __name__ == "__main__":
    repeated_test_entangled_linear_op()
