# This is a pytest script to test multibit_negative() and
# multibit_subtractor() functions from part 2.
# These functions should be imported from the solution file.

from hw1.p2 import multibit_negative, multibit_subtractor

# -----------------------------------------------------------------------------
# Helper functions to convert between integers and binary lists:
# * int2list: Converts a Python integer to a list of 0s and 1s with a
#             specified number of bits, where the least significant bit (LSB)
#             is at the beginning of the list.
# * list2int: Converts a list of 0s and 1s back to a Python integer, assuming
#             the least significant bit (LSB) is at the beginning of the list.

def int2list(a, nbits=4):
    """Convert an integer to a list of 0s and 1s (binary
    representation) with a specified number of bits, where the least
    significant bit (LSB) is at the beginning of the list.

    Args:
        a (int): The integer to convert.
        nbits (int): The number of bits for the output binary list.

    Returns:
        list[int]: The binary representation of the integer as a list
                   of 0s and 1s, with the least significant bit (LSB)
                   as the first element.

    """
    return [int(c) for c in bin(a)[2:].zfill(nbits)[::-1]]

def list2int(A, nbits=None):
    """Convert a binary list (list of 0s and 1s) to a signed integer, assuming
    the least significant bit (LSB) is at the beginning of the list.

    Args:
        A (list[int]): The binary list to convert, with the LSB as the
            first element.
        nbits (int, optional): The number of bits to interpret the
            binary list.  Defaults to the length of A.

    Returns:
        int: The signed integer representation of the binary list.

    """
    if nbits is None: nbits = len(A)
    a = int(''.join(str(c) for c in A[::-1]), 2)
    return a if a < 2**(nbits-1) else a - 2**nbits

# -----------------------------------------------------------------------------
# Test functions for multibit_negative() and multibit_subtractor().
# These functions use the helper functions for cleaner test code.

def test_negative():
    """Test the multibit_negative() function by verifying that it
    correctly computes the 2's complement (negation) of integers
    within the test range.

    The input and output are binary lists where the least significant
    bit (LSB) is at the beginning of the list.

    """
    for a in range(0, 16):
        A = int2list(a, nbits=5)

        C = multibit_negative(A)
        c = list2int(C)

        # Debugging output
        print(f"Input: {A} (LSB first), Negated: {C} -> {c}")

        # Assert the result matches the expected value
        assert c == -a, f"Failed for input {a}"

def test_subtractor():
    """Test the multibit_subtractor() function by verifying that it correctly
    computes the subtraction of two integers within the test range.

    The input and output are binary lists where the least significant
    bit (LSB) is at the beginning of the list.

    """
    for a in range(0, 8):
        for b in range(0, 8):
            A = int2list(a)
            B = int2list(b)

            C = multibit_subtractor(A, B)
            c = list2int(C)

            # Debugging output
            print(f"Input: {A} - {B} (LSB first), Result: {C} -> {c}")

            # Assert the result matches the expected value
            assert c == a - b, f"Failed for inputs {a} - {b}"
