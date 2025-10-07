from scipy import constants
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

class CanonicalEnsemble:
    """
    Defines a canonical ensemble with partition function Z
    """
    def __init__(self, z):
        """
        CanonicalEnsemble(z)

        Parameters
        ---------------------------
        z: partition function at a given temperature (K)
        """
        self.z = z

    def calc_u(self, temp_1, temp_2, n=100):
        """
        Calculates the internal energy of the system from temp_1 to temp_2 for n points

        Parameters
        ----------------------------
        temp_1 (K): starting temperature
        temp_2 (K): ending temperature
        n: number of datapoints

        Returns
        -----------------------------
        An np.array of calculated internal energies from temp_1 to temp_2 with length n (J)
        """
        temps = np.linspace(temp_1, temp_2, n)
        lnZs = np.log(np.array([self.z(temp) for temp in temps]))
        betas = 1/(constants.k*temps)
        us = -np.gradient(lnZs, betas)
        self.u = us
        return us

    def calc_f(self, temp_1, temp_2, n=100):
        """
        Calculates the free energy of the system from temp_1 to temp_2 for n points

        Parameters
        ----------------------------
        temp_1 (K): starting temperature
        temp_2 (K): ending temperature
        n: number of datapoints

        Returns
        -----------------------------
        An np.array of calculated free energies from temp_1 to temp_2 with length n (J)
        """
        temps = np.linspace(temp_1, temp_2, n)
        lnZs = np.log(np.array([self.z(temp) for temp in temps]))
        fs = -constants.k*temps*lnZs
        self.f = fs
        return fs

    def calc_s(self, temp_1, temp_2, n=100):
        """
        Calculates the entropy of the system from temp_1 to temp_2 for n points

        Parameters
        ----------------------------
        temp_1 (K): starting temperature
        temp_2 (K): ending temperature
        n: number of datapoints

        Returns
        -----------------------------
        An np.array of calculated entropies from temp_1 to temp_2 with length n (J/K)
        """
        temps = np.linspace(temp_1, temp_2, n)
        fs = self.calc_f(temp_1, temp_2, n)
        ss = -np.gradient(fs, temps)
        self.s = ss
        return ss

    def calc_thermo(self, temp_1, temp_2, n=100):
        """
        Calculates the internal energy, free energy and entropy of the system from temp_1 to temp_2 for n points
        Result is stored in attributes u, f and s, respectively.

        Parameters
        ----------------------------
        temp_1 (K): starting temperature
        temp_2 (K): ending temperature
        n: number of datapoints
        """
        self.calc_u(temp_1, temp_2, n)
        self.calc_f(temp_1, temp_2, n)
        self.calc_s(temp_1, temp_2, n)
        self.thermo_df = pd.DataFrame({
            'temperature': np.linspace(temp_1, temp_2, n),
            'internal_energy': self.u,
            'free_energy': self.f,
            'entropy': self.s
        })

def degenerate_partition(temp):
    """
    Calculates the partition function for an isolated Ce^{3+} ion with 14-fold degeneracy

    Returns:
    The partition function
    """
    return 14

isolatedCe = CanonicalEnsemble(degenerate_partition)
isolatedCe.calc_thermo(300, 2000, 100)

def soc_partition(temp):
    """
    Calculates the partition function for a Ce^{3+} ion with spin orbit coupling
    
    Parameters
    -----------------------------
    temp (K): temperature

    Returns:
    The partition function
    """
    e = 0.28*constants.eV # J
    exponent = -e/(constants.k*temp)
    return 6 + 8 * (np.e**exponent)

socCe = CanonicalEnsemble(soc_partition)
socCe.calc_thermo(300, 2000, 100)

def cfs_partition(temp):
    """
    Calculates the partition function for a Ce^{3+} ion with both spin orbit coupling and crystal field splitting
    
    Parameters
    -----------------------------
    temp (K): temperature
    
    Returns:
    The partition function
    """

    # convert eV to J
    def evToJ(ev):
        return ev*constants.eV
    # calculate the exponent
    def calc_exp(energy):
        # energy in eV
        return np.e**(-evToJ(energy) / (constants.k*temp))

    e0 = calc_exp(0)
    e1 = calc_exp(0.12)
    e2 = calc_exp(0.12+0.13)
    e3 = calc_exp(0.12+0.13+0.07)
    e4 = calc_exp(0.12+0.13+0.07+0.14)
    return 4*e0 + 2*e1 + 2*e2 + 4*e3 + 2*e4

cfsCe = CanonicalEnsemble(cfs_partition)
cfsCe.calc_thermo(300, 2000, 100)

# combine the calculated properties into a single df
combined_df = pd.concat([
    isolatedCe.thermo_df.assign(system='isolated'),
    socCe.thermo_df.assign(system='soc'),
    cfsCe.thermo_df.assign(system='soc and cfs')
], axis=0)
combined_df.to_csv('thermo_properties.csv', index=False)
print("Saved combined thermodynamic properties to thermo_properties.csv")

# plot internal energies
fig, ax = plt.subplots()
fig.set_size_inches(7,4)
sns.lineplot(combined_df, x='temperature', y='internal_energy', hue='system', ax=ax)
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Internal Energy (J)')
ax.set_title('Internal energies for different Ce$^{3+}$ systems vs. temperature')
fig.savefig('internal_energies.png', dpi=210)
print("Saved plot of internal energies to internal_energies.png")

# plot free energies
fig, ax = plt.subplots()
fig.set_size_inches(7,4)
sns.lineplot(combined_df, x='temperature', y='free_energy', hue='system', ax=ax)
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Free Energy (J)')
ax.set_title('Free energies for different Ce$^{3+}$ systems vs. temperature')
fig.savefig('free_energies.png', dpi=210)
print("Saved plot of free energies to free_energies.png")

# plot entropies
fig, ax = plt.subplots()
fig.set_size_inches(7,4)
sns.lineplot(combined_df, x='temperature', y='entropy', hue='system', ax=ax)
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Entropy (J/K)')
ax.set_title('Entropies for different Ce$^{3+}$ systems vs. temperature')
fig.savefig('entropies.png', dpi=210)
print("Saved plot of entropies to entropies.png")