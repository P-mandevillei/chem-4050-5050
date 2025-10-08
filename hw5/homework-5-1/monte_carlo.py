import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import expon, gamma

def psi_2p_z(x, y, z):
    """
    Computes the 2pz wavefunctions for a hydrogen atom.

    Parameters
    -----------------------------
    x, y, z (atomic unit): (float or np.array) coordinates

    Returns
    Hydrogen 2pz wavefunction at (x, y, z)
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    a0 = 1
    coeff = 1 / (4 * np.sqrt(2*np.pi) * (a0**(3/2)))
    # prevent divide by 0
    if type(r) == 'numpy.float64' and r == 0:
        r = 1e-12
    else:
        r[r == 0] = 1e-12
    angular = z/r
    radial = (r/a0) * (np.e**(-r/(2*a0)))
    return coeff * angular * radial

def uniform_mc(n, sep=2, length=20):
    """
    Uniform Monte Carlo sampling

    Parameters
    ------------------------------
    n: number of point estimators
    sep (atomic unit): separation distance of two hydrogens
    length (atomic unit): defines the cube (-L, L) over which to integrate
    
    Returns
    Monte Carlo estimation of the overlap of hydrogen 2pz orbitals
    """
    print(f"Monte Carlo with {n} estimators")
    memory_lim = int(1e5)
    rng = np.random.default_rng(114514)
    if n <= memory_lim:
        x, y, z = rng.uniform(low=0, high=length, size=(3, n)) # only need to integrate over the first quadrant due to symmetry
        integrands = psi_2p_z(x, y, z+sep/2)*psi_2p_z(x, y, z-sep/2)
        volume = (2*length)**3
        overlap = volume*np.mean(integrands)
        return overlap
    else:
        # prevent memory overflow!
        # MC in chunks
        chunks = n // memory_lim
        remainder = n % memory_lim
        integrands = []
        for chunk in tqdm(range(chunks)):
            x, y, z = rng.uniform(low=0, high=length, size=(3, memory_lim))
            integrands.append(psi_2p_z(x, y, z+sep/2)*psi_2p_z(x, y, z-sep/2))
        if remainder>0:
            x, y, z = rng.uniform(low=0, high=length, size=(3, remainder))
            integrands.append(psi_2p_z(x, y, z+sep/2)*psi_2p_z(x, y, z-sep/2))
        integrands = np.concat(integrands)
        volume = (2*length)**3
        overlap = volume*np.mean(integrands)
        return overlap

def expon_mc(n, sep=2, length=20):
    """
    Monte Carlo importance sampling

    Parameters
    ------------------------------
    n: number of point estimators
    sep (atomic unit): separation distance of two hydrogens
    length (atomic unit): defines the cube (-L, L) over which to integrate
    
    Returns
    Monte Carlo estimation of the overlap of hydrogen 2pz orbitals using importance sampling
    """
    rng = np.random.default_rng(114514)
    print(f"Monte Carlo with {n} estimators")
    memory_lim = int(1e5)
    if n <= memory_lim:
        x, y = expon.rvs(size=(2, n), scale=4, random_state=rng) # exponential distribution has a support of [0, +infty)
        z = gamma.rvs(size=n, a=3, loc=0, scale=2, random_state=rng) # also support of [0,infty)
        # only need to integrate over the first quadrant due to symmetry
        integrands = psi_2p_z(x, y, z+sep/2)*psi_2p_z(x, y, z-sep/2)
        g = expon.pdf(x, scale=4) * expon.pdf(y, scale=4) * gamma.pdf(z, a=3, loc=0, scale=2) # independent variables
        integrands /= g
        overlap = 8*np.mean(integrands) # eight quadrants
        return overlap
    else:
        # prevent memory overflow!
        # MC in chunks
        chunks = n // memory_lim
        remainder = n % memory_lim
        integrands = []
        for chunk in tqdm(range(chunks)):
            x, y = expon.rvs(size=(2, memory_lim), scale=4, random_state=rng)
            z = gamma.rvs(size=memory_lim, a=3, loc=0, scale=2, random_state=rng)
            wavefunctions = psi_2p_z(x, y, z+sep/2)*psi_2p_z(x, y, z-sep/2)
            g = expon.pdf(x, scale=4) * expon.pdf(y, scale=4) * gamma.pdf(z, a=3, loc=0, scale=2) 
            integrands.append(wavefunctions/g)
        if remainder>0:
            x, y = expon.rvs(size=(2, remainder), scale=4, random_state=rng)
            z = gamma.rvs(size=remainder, a=3, loc=0, scale=2, random_state=rng)
            wavefunctions = psi_2p_z(x, y, z+sep/2)*psi_2p_z(x, y, z-sep/2)
            g = expon.pdf(x, scale=4) * expon.pdf(y, scale=4) * gamma.pdf(z, a=3, loc=0, scale=2) 
            integrands.append(wavefunctions/g)
        integrands = np.concat(integrands)
        overlap = 8*np.mean(integrands)
        return overlap