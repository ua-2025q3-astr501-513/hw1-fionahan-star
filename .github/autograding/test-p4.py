# This is a pytest script to test CoupledOscillators from part 4.
# This class should be imported from the solution file.

import pytest
import numpy as np

from hw1.p4 import CoupledOscillators

@pytest.mark.parametrize('params', [
    # Test cases with varying numbers of oscillators (N) and expected results (X1)
    {'N':2,  'X1':[0.17521471,0.09493644]},
    {'N':4,  'X1':[0.00060096,0.01698677,0.1758158 ,0.09494753]},
    {'N':8,  'X1':[5.36460880e-12,9.66325041e-10,1.25763949e-07,1.10888098e-05,
                   6.01088717e-04,1.69867701e-02,1.75815800e-01,9.49475297e-02]},
    {'N':16, 'X1':[1.32933324e-16,-3.11747666e-17,-4.10060005e-17,-4.85503047e-17,
                   8.74352901e-17, 5.24443263e-18, 1.21991521e-16, 2.24148784e-14,
                   5.36485433e-12, 9.66324942e-10, 1.25763949e-07, 1.10888098e-05,
                   6.01088717e-04, 1.69867701e-02, 1.75815800e-01, 9.49475297e-02]},
])
def test_Oscillators(params):
    """
    Test the CoupledOscillators class for correctness of displacements
    at a fixed time (t=1) for various numbers of oscillators (N).

    Args:
        params (dict): A dictionary containing:
            - 'N': Number of oscillators in the system.
            - 'X1': Expected displacements of the oscillators at t=1.
    """
    # Initial displacements: All zeros except the last oscillator
    X0 = np.zeros(params['N'])
    X0[-1] = 0.5

    # Initialize the CoupledOscillators system
    co = CoupledOscillators(X0)

    # Compute displacements at t=1
    X1 = co(1)

    # Assert that the computed displacements match the expected values
    assert np.allclose(X1, params['X1']), (
        f"Failed for N={params['N']}. "
        f"Expected: {params['X1']}, but got: {X1}"
    )
