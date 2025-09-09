#!/usr/bin/env python3
#
# Please look for "TODO" in the comments, which indicate where you
# need to write your code.
#
# Part 4: Solve the Coupled Simple Harmonic Oscillator Problem (1 point)
#
# * Objective:
#   Take the coupled harmonic oscillator problem we solved in class
#   and rewrite it using a well-structured Python class.
# * Details:
#   The description of the problem and the solution template can be
#   found in `hw1/p4.py`.
#
# From lecture `02w`, we solve systems of coupled harmonic oscillators
# semi-analytically by numerically solving eigenvalue problems.
# However, the code structure was not very clean, making the code hard
# to reuse.
# Although numerical analysis in general does not require
# object-oriented programming, it is sometime useful to package
# stateful caluation into classes.
# For this assignment, we will provide a template class.
# Your responsibility to implement the methods in the class.


import numpy as np


class CoupledOscillators:
    """A class to model a system of coupled harmonic oscillators.

    Attributes:
        Omega (np.ndarray): array of angular frequencies of the normal modes.
        V     (np.ndarray): matrix of eigenvectors representing normal modes.
        M0    (np.ndarray): initial amplitudes of the normal modes.

    """

    def __init__(self, X0=[-0.5, 0, 0.5], m=1.0, k=1.0):
        """Initialize the coupled harmonic oscillator system.

        Args:
            X0 (list or np.ndarray): initial displacements of the oscillators.
            m  (float):              mass of each oscillator (assumed identical for all oscillators).
            k  (float):              spring constant (assumed identical for all springs).

        """
        # TODO: Construct the stiffness matrix K
        self.dims = len(X0)
        K = np.zeros((self.dims, self.dims))
        for j in range(self.dims):
            if j != 0:
                K[j][j-1] = 1
                
            K[j][j] = -2
            
            if j != self.dims-1:
                K[j][j+1] = 1
                
        stiff_K = K
        
        # TODO: Solve the eigenvalue problem for K to find normal modes
        # Here we call on the numpy module to do the heavylift. 
        l, v = np.linalg.eig(stiff_K)
        
        # TODO: Store angular frequencies and eigenvectors
        # v contains the eigenvalues and l contains the eigenvectors. 
        # To translate between eigenvalues and angular frequencies, note that the eigenvalues we found
        # are actually m * omega**2. Retrieve the angular frequencies as such. 
        self.Omega = np.sqrt(abs(l) * k / m)
        self.V = v
        # The eigenvectors are unchanged. Although they are a matrix now.
        # We can use these eigenvectors to find the decompositions back to x (non-normal coordinates).
        
        # TODO: Compute initial modal amplitudes M0 (normal mode decomposition)
        # eig gives normalized eigenvectors so we should be able to do (x dot v)*v summed up?
        M0 = np.zeros(self.dims)
        for i in range(self.dims):
            # Dot product X0 and each of the eigenvectors for component on that direction. 
            temp = np.dot(X0, v[:, i])
            M0[i] = temp
        # and done. 
        self.M0 = M0

    def __call__(self, t):
        """Calculate the displacements of the oscillators at time t.

        Args:
            t (float): time at which to compute the displacements.

        Returns:
            np.ndarray: displacements of the oscillators at time t.

        """
        # TODO: Reconstruct the displacements from normal modes
        # To switch back to the original basis use the change of basis matrix. Which happened to be V. 
        cobmat = self.V
        # And first we evolve things in the normal coordinate.
        locarr = np.zeros(self.dims)
        for i in range(self.dims):
            locarr[i] = self.M0[i] * np.cos(self.Omega[i] * t)
        # Now we change of basis back:
        X = np.matmul(cobmat, locarr)
        return X

if __name__ == "__main__":

    # Initialize the coupled oscillator system with default parameters
    co = CoupledOscillators()

    # Print displacements of the oscillators at each time step
    print("Time(s)  Displacements")
    print("----------------------")
    for t in np.linspace(0, 10, num=101):
        X = co(t)             # compute displacements at time t
        print(f"{t:.2f}", X)  # print values for reference