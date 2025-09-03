# This is a pytest script to test quadratic() from part 3.
# This function should be imported from the solution file.

from hw1.p3 import quadratic

# -----------------------------------------------------------------------------
# Helper functions

from math import isclose

def agree(x, t):
    if t is None:
        return x is None
    else:
        return isclose(x, t)

def compare(s, t1, t2, a=1, b=2, c=1):
    x1, x2 = quadratic(a, b, c)

    print(f"[{s}] Input: {a} {b} {c}, {t1} {t2}, Result: {x1} {x2}")

    assert agree(x1, t1)
    assert agree(x2, t2)

# -----------------------------------------------------------------------------
# Test functions for quadratic().

def test_catastrophic_cancellation():
    """Test the quadratic() function by verifying that it can compute
    accurately both roots when catastrophic cancellation may happen.

    """
    for t1 in [-1e8, -1e7, -1e6, -1e5, -1e4, 1e4, 1e5, 1e6, 1e7, 1e8]:
        t2 = 1.0 / t1
        b  = -t1-t2
        compare("Catastrophic Cancellation", min(t1, t2), max(t1, t2), b=b)

def test_one_root():
    """Test the quadratic() function by verifying x2 == None when
    there is only one real root.

    """
    for t1 in range(10):
        b = -2 * t1
        c = t1 * t1
        compare("One Root", t1, None, b=b, c=c)

def test_no_root():
    """Test the quadratic() function by verifying x1 == x2 == None
    when there is no real root.

    """
    for b in [-1e-0, -1e-1, -1e-2, -1e-3, -1e-4,
               1e-4,  1e-3,  1e-2,  1e-1,  1e-0]:
        compare("No Root", None, None, b=b)
