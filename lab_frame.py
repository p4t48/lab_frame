from __future__ import division
import qutip as qt
import numpy as np
import matplotlib.pyplot as plt


# Setting parameters for the simulation

w0 = 2 * np.pi
w1 = 0.5 * 2 * np.pi
T1 = 0.3
T2 = 0.3
M0 = 1

# Definition of the steady-state solutions of the Bloch equations for magnetization M = [Mx, My, Mz]
# w is the frequency of the applied RF-field

def Mx(w):
    nominator = w1 * T2**2 * (w0 - w)
    denominator = 1 + w1**2 * T1 * T2 + (T2 * (w0 - w))**2

    return nominator/denominator * M0

def My(w):
    nominator = w1 * T2
    denominator = 1 + w1**2 * T1 * T2 + (T2*(w0 - w))**2
    
    return nominator/denominator * M0

def Mz(w):
    nominator = 1 + T2**2 * (w0 - w)**2
    denominator = 1 + w1**2 * T1 * T2 + (T2 * (w0 - w))**2

    return nominator/denominator * M0


wVals = np.linspace(-2 * np.pi, 6 * np.pi, 100)

# Solving the dynamics of a spin 1/2 in B0 and Brf in the lab frame

def H(w):
    """ Hamiltonian in the lab frame. """
    return w0 / 2 * qt.sigmaz() + w1 * np.cos(w) * qt.sigmax()

c1 = np.sqrt(1.0/T1) * qt.destroy(2)
c2 = np.sqrt(1.0/T2) / 2 * qt.destroy(2)
c_op_list = [c1,c2]


wHVals = np.linspace(-2 * np.pi, 6 * np.pi, 51)
MxH = [qt.expect(qt.sigmax(), qt.steadystate(H(w), c_op_list, maxiter=100)) for w in wHVals]
MyH = [qt.expect(qt.sigmay(), qt.steadystate(H(w), c_op_list, maxiter=100)) for w in wHVals]
MzH = [qt.expect(qt.sigmaz(), qt.steadystate(H(w), c_op_list, maxiter=100)) for w in wHVals]


# Plot results of the algebraic formula and the solution given by qutip

plt.ylim([-1.0, 1.0])

plt.plot(wVals, Mx(wVals), label='Mx')
plt.plot(wVals, My(wVals), label='My')
plt.plot(wVals, Mz(wVals), label='Mz')

plt.plot(wHVals, MxH, 'bo', label='MxH')
plt.plot(wHVals, MyH, 'go', label='MyH')
plt.plot(wHVals, MzH, 'ro', label='MzH')

plt.legend(loc=3, numpoints=1)
plt.show()



# Define parameters for the state and Hamiltonian
"""
w0 = 2 * np.pi
w1 = .1 * 2 * np.pi
w = 2 * np.pi



print H.eigenenergies()
"""
