# imports
from scipy.optimize import minimize
from matplotlib import pyplot as plt
import numpy as np
import math, warnings

def lennard_jones(r, epsilon=0.01, sigma=3.4):
    """
    Calculates the Lennard_Jones potential.

    Parameters:
        r (angstrom): interatomic distance
        epsilon (eV): depth of the potential well (minimum energy)
        sigma (angstrom): zero potential distance

    Returns:
        Lennard_Jones potential
    """
    V = 4*epsilon*((sigma/r)**12 - (sigma/r)**6)
    return V

# Optimize the geometry of Ar2
guess = 4
Vmin = minimize(
    lennard_jones,
    guess,
    method="Nelder-Mead",
    tol=1e-6
)
r = Vmin['x'][0]
print("Successfully optimized Ar2 geometry!")

# plot the potential energy curve
x = np.linspace(3, 6, 500)
y = lennard_jones(x)
fig, ax = plt.subplots(figsize=(8,5))

ax.plot(x, y, label="Lennard-Jones Potential")
ax.axvline(x = r, linestyle='--', linewidth=1, color='r', alpha=0.8, label=rf"Vmin Distance ({r:.2f}$\AA$)")
ax.set_xlabel(r'Interatomic Distance ($\AA$)')
ax.set_ylabel('Potential Energy (eV)')
ax.set_title("Lennard-Jones Potential of Ar$_2$")
ax.grid()
ax.legend()

fig.savefig('homework-2-1/argon_dimer_potential_vs_distance.png', dpi=200, bbox_inches='tight')

# report
print(f"Optimized atomic distances: r = {r:.2f}")

# storing geometry
# Assume y=z=0
geometry = f"""2
Argon dimer
Ar\t{0:.6f}\t{0:.6f}\t{0:.6f}
Ar\t{Vmin['x'][0]:.6f}\t{0:.6f}\t{0:.6f}
"""

with open('homework-2-1/argon_dimer.xyz', 'w') as file:
    file.write(geometry)
print("Saved homework-2-1/argon_dimer.xyz")