#!/usr/bin/env python3
#
# Please look for "TODO" in the comments, which indicate where you
# need to write your code.
#
# Part 2: Integer Negation and Subtraction Using NAND Gates (1 point)
#
# * Objective:
#   Implement a function that performs integer negation using only NAND
#   gates and use it to implement subtraction.
# * Details:
#   The description of the problem and the solution template can be
#   found in `hw1/p2.py`.
#
# From lecture `01w`, we learned that NAND is a universal gate, that
# any binary operations can be built by using only NAND gates.
# Following the lecture notes, we define the "NAND gate" as

def NAND(a, b):
    return 1 - (a & b)  # NOT (a AND b)

# Following the notes again, we define also other basic operations:

def NOT(a):
    return NAND(a, a)

def AND(a, b):
    return NOT(NAND(a, b))

def OR(a, b):
    return NAND(NOT(a), NOT(b))

def XOR(a, b):
    c = NAND(a, b)
    return NAND(NAND(a, c), NAND(b, c))

# We also implemented the half, full, and multi-bit adders:

def half_adder(A, B):
    S = XOR(A, B)  # Sum using XOR
    C = AND(A, B)  # Carry using AND
    return S, C

def full_adder(A, B, Cin):
    s, c = half_adder(A,   B)
    S, C = half_adder(Cin, s)
    Cout = OR(c, C)
    return S, Cout

def multibit_adder(A, B, carrybit=False):
    assert(len(A) == len(B))

    n = len(A)
    c = 0
    S = []
    for i in range(n):
        s, c = full_adder(A[i], B[i], c)
        S.append(s)
    if carrybit:
        S.append(c)  # add the extra carry bit
    return S

# Now, getting into the assignment, we would like to first implement a
# negative function.
#
# Please keep the following function prototype, otherwise the
# auto-tester would fail, and you will not obtain point for this
# assigment.

def multibit_negative(A):
    """Multi-bit integer negative operator

    This function take the binary number A and return negative A using
    two's complement.
    In other words, if the input
        A = 3 = 0b011,
    then the output is
        -A = -3 = 0b101.

    Args:
        A: input number in binary represented as a python list, with
           the least significant digit be the first.
           That is, the binary 0b011 should be given by [1,1,0].

    Returns:
        Negative A using two's complement represented as a python
        list, with the least significant digit be the first.

    """
    # TODO: implement the function here

# We are now ready to implement subtraction using multibit_adder() and
# multibit_negative().

def multibit_subtractor(A, B):
    """Multi-bit integer subtraction operator

    This function take the binary numbers A and B, and return A - B.
    Be careful on how the carrying bit is handled in multibit_adder().
    Make sure that when A == B, the result A - B should be zero.

    Args:
        A, B: input number in binary represented as a python list,
           with the least significant digit be the first.
           That is, the binary 0b011 should be given by [1,1,0].

    Returns:
        A - B represented as a python list, with the least significant
        digit be the first.

    """
    # TODO: implement the function here
