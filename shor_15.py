import kron as kron_ops
import gates
import expansion
import random
import functools
import math
import numpy as np
import scipy.linalg as sci_lin
import itertools

"""
Tests w/ deterministic qubits vs. expand-all-the-amplitudes slow way
They should agree, though.

The deterministic qubits aren't any faster if you're expanding them to have all amplitudes

You need to apply XORSAT to actually have the speedup. This will happen in the speed test

This is all hardcoded for (7 ** x) mod 15.
There will be lots of tedious work to get it to be (a ** x) mod b, but I don't think anything that's not going off and reading a paper and translating it into propositional variables.

This also incorporates some but not all of the many optimizations and tricks the group of https://arxiv.org/abs/quant-ph/0112176 did to get their thing in 7 qubits

I'm not sure why this thing still works with 0000 rather than 0001, but it does. Maybe a substantive mistake
"""


def bells(assignments):
    # Bells state, like quantum state, entangling basically ends up being "having XOR's in it".
    # As opposed to product state which is composed of vars linked together w/ AND's only
    # Think of the AND w/ tautologies as noops, because they are
    return [1, 1][assignments["A"]] * (
        [0, 1][assignments["A"]] ^ [0, 1][assignments["B"]]
    )


def test_bells():
    init_state = [[1, 1], [0, 1]]
    expanded = kron_ops.expand(init_state, ["*"])
    expanded_bell = gates.expanded_cnot2().dot(expanded)
    print(expanded_bell)

    expanded_res = expansion.expand_res_by_function(bells, 2)
    print(expanded_res)


def toffoli(assignments):
    gate_1 = assignments["A"] ^ assignments["B"]
    gate_2 = (assignments["A"] & gate_1) ^ assignments["C"]
    res = [1, 1][assignments["A"]] * gate_1 * gate_2

    return res


def test_toffoli():
    init_state = [[1, 1], [0, 1], [0, 1]]
    expanded = kron_ops.expand(init_state, ["*", "*"])
    expanded_1 = gates.expanded_toffoli(0, 1, 2, len(init_state)).dot(expanded)
    expanded_2 = gates.expanded_cnot(0, 1, len(init_state)).dot(expanded_1)

    print("mat expansion 0: ", expanded)
    print("mat expansion 1: ", expanded_1)
    print("mat expansion 2: ", expanded_2)

    expanded_res = expansion.expand_res_by_function(toffoli, 3)
    print("elementwise expansion: ", expanded_res)


def shors_forwards(assignments):
    # shors mod exp:
    # counting from 1, not 0
    # cnot gate 3 => 5
    # cnot gate 3 => 6
    # cnot gate 4 => 6
    # toffoli gate 2,6 => 4
    # cnot gate 4 => 6
    # cnot gate 7 => 5
    # toffoli gate 2,5 => 7
    # cnot gate 7 => 5

    gate_1 = assignments["C"] ^ assignments["E"]
    gate_2 = assignments["C"] ^ assignments["F"]
    gate_3 = assignments["D"] ^ gate_2
    gate_4 = (assignments["B"] * gate_3) ^ assignments["D"]
    gate_5 = gate_4 ^ gate_3
    gate_6 = assignments["G"] ^ gate_1
    gate_7 = (assignments["B"] * gate_6) ^ assignments["G"]
    gate_8 = gate_7 ^ gate_6

    # Really, A, B, C are noops but included here for completeness
    res = (
        [1, 1][assignments["A"]]
        * [1, 1][assignments["B"]]
        * [1, 1][assignments["C"]]
        * (gate_4)
        * (gate_8)
        * (gate_5)
        * (gate_7)
    )

    return res


def test_shors_forwards():
    # Hadamard gate pre-applied to make tautologies
    init_state = [[1, 1], [1, 1], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]]
    expanded = kron_ops.expand(init_state, ["*" for _ in range(6)])
    # expanded_cnot counts from 0
    # noticed the reversed order of gates as opposed to functional thing
    expanded_shor1 = gates.expanded_cnot(6, 4, len(init_state)).dot(expanded)
    expanded_shor2 = gates.expanded_toffoli(1, 4, 6, len(init_state)).dot(
        expanded_shor1
    )
    expanded_shor3 = gates.expanded_cnot(6, 4, len(init_state)).dot(expanded_shor2)
    expanded_shor4 = gates.expanded_cnot(3, 5, len(init_state)).dot(expanded_shor3)
    expanded_shor5 = gates.expanded_toffoli(1, 5, 3, len(init_state)).dot(
        expanded_shor4
    )
    expanded_shor6 = gates.expanded_cnot(3, 5, len(init_state)).dot(expanded_shor5)
    expanded_shor7 = gates.expanded_cnot(2, 5, len(init_state)).dot(expanded_shor6)
    expanded_shor8 = gates.expanded_cnot(2, 4, len(init_state)).dot(expanded_shor7)

    print(expanded_shor8)

    expanded_res = expansion.expand_res_by_function(shors_forwards, 7)
    print(expanded_res)
    print(np.nonzero(expanded_res))

    print(np.allclose(expanded_shor8, expanded_res))


def shors_backwards(assignments):
    # shors mod exp:
    # counting from 1, not 0
    # cnot gate 3 => 5
    # cnot gate 3 => 6
    # cnot gate 4 => 6
    # toffoli gate 2,6 => 4
    # cnot gate 4 => 6
    # cnot gate 7 => 5
    # toffoli gate 2,5 => 7
    # cnot gate 7 => 5

    # now, this one goes backwards, so gate 1 here is gate 8 forwards, etc
    gate_1 = assignments["G"] ^ assignments["E"]
    gate_2 = (assignments["B"] * gate_1) ^ assignments["G"]
    gate_3 = gate_2 ^ gate_1
    gate_4 = assignments["D"] ^ assignments["F"]
    gate_5 = (assignments["B"] * gate_4) ^ assignments["D"]
    gate_6 = gate_5 ^ gate_4
    gate_7 = assignments["C"] ^ gate_6
    gate_8 = assignments["C"] ^ gate_3

    # Really, A, B, C are noops but included here for completeness
    res = (
        [1, 1][assignments["A"]]
        * [1, 1][assignments["B"]]
        * [1, 1][assignments["C"]]
        * (gate_5)
        * (gate_8)
        * (gate_7)
        * (gate_2)
    )

    return res


def test_shors_backwards():
    init_state = [[1, 1], [1, 1], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]]
    expanded = kron_ops.expand(init_state, ["*" for _ in range(6)])
    # reversed of reversed is forwards
    expanded_shor1 = gates.expanded_cnot(2, 4, len(init_state)).dot(expanded)
    expanded_shor2 = gates.expanded_cnot(2, 5, len(init_state)).dot(expanded_shor1)
    expanded_shor3 = gates.expanded_cnot(3, 5, len(init_state)).dot(expanded_shor2)
    expanded_shor4 = gates.expanded_toffoli(1, 5, 3, len(init_state)).dot(
        expanded_shor3
    )
    expanded_shor5 = gates.expanded_cnot(3, 5, len(init_state)).dot(expanded_shor4)
    expanded_shor6 = gates.expanded_cnot(6, 4, len(init_state)).dot(expanded_shor5)
    expanded_shor7 = gates.expanded_toffoli(1, 4, 6, len(init_state)).dot(
        expanded_shor6
    )
    expanded_shor8 = gates.expanded_cnot(6, 4, len(init_state)).dot(expanded_shor7)
    print(expanded_shor8)

    expanded_res = expansion.expand_res_by_function(shors_backwards, 7)
    print(expanded_res)
    print(np.nonzero(expanded_res))

    print(np.allclose(expanded_shor8, expanded_res))


def test_sat_shors():
    """
    We have assignments, we have outputs from some previous run and we're getting possible inputs for those outputs, and then looking at periodicity of inputs

    So we got clauses, set all those clauses equal to 1. that's what needs to happen for overall to eq 1.
    Then, we apply all gates in reverse order. We get to do this because of involutive character of Toffoli and CNOT gate. This is guaranteed to terminate and not have cycles because the forward way terminated and didn't have cycles. Basically sort of like Gaussian elimination for XORSAT but with entire (implicit) linear operators instead of elementary row operations.

    You might think this might lead to funny business with ordinary SAT
    (because Toffoli gate ends up being a classical NAND)
    Conversion to Toffoli gate is the explosion in doing ordinary SAT w/ this, I think? With having too many ancilla gates?
    I'm not sure. There might be funny business.

    Toffoli inversion may be a problem when I scale it up, because of that AND in there. That may be a reason to go and pull out the QFT again, or consider the fact that Toffoli gate can be represented with 6 other gates, all CNOTs and single-qubit gates.
    
    Takes the gates (these are baked in here) and an output (also baked in), return an assignment
    """
    assignments = {}
    assignments["D"] = 1
    assignments["E"] = 1
    assignments["F"] = 1
    assignments["G"] = 1

    # Gates 5, 8, 7, 2 are clauses set to 1, thinking in backwards way

    # Gate 5 ends up being, (B AND (D XOR F)) XOR F = 1
    # Ends up being (B AND 0) XOR 1 = 1, which is tautology 

    # Gate 8 ends up being, (C XOR E XOR (B AND (G XOR E))) = 1
    # C XOR 1 XOR (B AND 0) = 1
    # C = 0

    # Gate 7 ends up being, C XOR F XOR (B AND (G XOR E)) = 1
    # C XOR 1 XOR (B AND 0) = 1
    # C = 0

    # Gate 2 ends up being (B AND (G XOR F)) XOR G = 1
    # Also tautology

    # So we get to the superposition with C = 0 and A, B in superposition still.
    # Unlike in quantum systems proper, we can just get the elements

    # Obviously in actual implementation we're going to have to go and not do the above by hand

    # Defines the samples from which we can draw the inputs to get GCD's from, (not actually samples)
    # A and B are superposition (tautology)
    assignments["A"] = [0, 1]
    assignments["B"] = [0, 1]
    assignments["C"] = [0]

    # Also in actual implementation, we're not going to do a Cartesian product, we're going to sample in an analogous way to actual Shor's
    possible_inputs = set()
    for sample in itertools.product(
            assignments["A"],
            assignments["B"],
            assignments["C"],
            ):
        int_form = 0
        for bit in sample:
            int_form = (int_form << 1) | bit
        possible_inputs.add(int_form)
    print("possible inputs: ", possible_inputs)
    print("Read the code to see how I got them")
    collision_pairs = set()
    for _ in range(5):
        collision_pair = (random.choice(list(possible_inputs)), random.choice(list(possible_inputs)))
        while (collision_pair[0] - collision_pair[1]) == 0:
            collision_pair = (random.choice(list(possible_inputs)), random.choice(list(possible_inputs)))
        collision_pairs.add(abs(collision_pair[0] - collision_pair[1]))
    print("collision pairs: ", collision_pairs)
    res = functools.reduce(math.gcd, collision_pairs)

    # This looks bad because we're finding period of 7^x mod 15, but the claim is that this scales up way better than number field sieve
    print("Periodicity in amplitude now shown to be:", 2**3 / res)


if __name__ == "__main__":
    print("You're better off reading the code, this is basically a halfassed literate program")
    test_sat_shors()
