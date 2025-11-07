<h1 style="text-align: center;"> Project 1 Report </h1>

This readme is for display purposes only. Please refer to the [pdf report document](./Chem_4050_Project_1_Report.pdf) for the actual report.

- [source notebook](./project_1.ipynb)
- [report pdf](./Chem_4050_Project_1_Report.pdf)

[back](../README.md)

---------------------------------


# Example system
![example phase diagram](output_figures/example_phase_diagram.png)

# Exploring adsorption of $N_2$ and $H_2$
- Phase diagrams:
![Habor phase diagrams](output_figures/haber_exploration_phase_diagrams.png)
- Example lattice configurations:
![Habor example lattices](output_figures/haber_exploration_lattices.png)

# Optional Enhancement

## 1: particle swaps
- Phase diagrams using the same conditions for $H_2$ and $N_2$:

![Habor exploration phase diagrams enhanced](output_figures/haber_exploration_enhanced_phase_diagrams.png)

- Example lattice configurations:

![Habor exploration lattices enhanced](output_figures/haber_exploration_lattices_enhanced.png)

We see the same trends in the resulting phase diagrams!

## 2: Complex geometries and interactions
A 4% Ru-Ba-K/C catalyst is an emerging catalyst for the Haber-Bosch process. Here we consider the adsorption of $H_2$ and $N_2$ on ruthenium 001 surfaces.
Ru (001) surface is hexagonal.<br>
[Jacobi](https://doi.org/10.1002/(SICI)1521-396X(200001)177:1%3C37::AID-PSSA37%3E3.0.CO;2-Y) found the binding energy of nitrogen gas to be -0.5 eV, with the nitrogen binding on the top site of Ru atom.
![N2-Ru-binding.png](report_figures/N2-Ru-binding.png)
[Ungerer and Leeuw](https://pubs.rsc.org/en/content/articlelanding/2025/cp/d4cp04165h) found the total adsorption energy of molecular $H_2$, including its dissociation into elemental H*, to be around -1.34 eV.
![H2-Ru-binding.png](report_figures/H2-Ru-binding.png)
![H-Ru-Binding.png](report_figures/H-Ru-Binding.png)
**For simplicity's sake, we assume both $H_2$ and $N_2$ bind to the top sites, with $\epsilon_{H_2}=-1.34eV$ and $\epsilon_{N_2}=-0.5eV$.**
We model the interaction between $H_2$ and $N_2$ using Lennard-Jones potential.
$$ V(r) = 4\epsilon\left[ (\frac{\sigma}{r})^{12} - (\frac{\sigma}{r})^{6} \right] $$
With data from [Wang et al](https://pubs-acs-org.libproxy.washu.edu/doi/pdf/10.1021/acs.jctc.0c01132):
![H2_N2_LJ.png](report_figures/H2_N2_LJ.png)
Lorentz-Berthelot mixing rules: $\sigma_{12} = \frac{\sigma_1+\sigma_2}{2}$ and $\epsilon_{12} = \sqrt{\epsilon_{1}\epsilon_{2}} $.<br>
We will use the following measures calculated using the first set of parameters (12-6 LJ):
| interaction | $\epsilon (eV)$ | $\sigma (Å)$ |
| ----------- | ------------------- | ------------ |
| $H_2$-$H_2$ | $6.6\times 10^{-4}$ |  2.918       |
| $N_2$-$N_2$ | $3.5\times 10^{-3}$ |  3.614       |
| $H_2$-$N_2$ | $1.5\times 10^{-3}$ |  3.266       |

![LJ comparisons](output_figures/lj_comparisons.png)

Ru bond length is calculated to be $2.645$ Å from the [mp-33 structure](https://next-gen.materialsproject.org/materials/mp-33?formula=Ru#crystal_structure).<br>
We will use a 4x4 lattice to model the adsorption process. As can be seen from the LJ curves above, the interactions are basically 0 when the particles are $2\times 2.645 = 5.29\AA$ away. This means using the minimum image convention is appropriate.<br>

We define the unit cell as a bounding parallelogram with $\alpha=60^{\cdot}$
![trigonal_lattice.jpg](report_figures/trigonal_lattice.jpg)
<br>
Simulation result:
![Ru H2-N2 binding](output_figures/ru_h2_n2_binding.png)

## 3: Effect of lattice size

Since we are using minimum image convention, I expect the lattice size to have no influence on adsorption behavior when it is more than 2 times the interaction cutoff. For the square lattice with neighbor-neighbor interactions (cutoff=1), this means the lattice should at least be 3x3.
Here is a comparison between the phase diagrams using different lattice sizes:
![effect of lattice size](output_figures/effect_of_size_phase_diagrams.png)
As expected, when n<2, the resulting phase diagram deviates from the others. When n>2, the resulting phase diagrams are consistent with each other.

## 4: Competition between 3 species
We include $NH_3$ in our model.<br>
[Danielson et al](https://www.sciencedirect.com/science/article/abs/pii/0039602878904508#:~:text=Abstract,100%20K%20will%20be%20discussed.) found that the desorption energy of $NH_3$ molecules from Ru(001) at low temperature (100K) can be 0.32 or 0.46 eV depending on its molecular states. For simplicity, we take the average and let its adsorption energy be $\epsilon_{NH_3}=-0.39eV$.<br>
For Lennard-Jones, we use $\sigma_{NH_3}=2.900$ and $\epsilon=4.8\times 10^{-2}eV$, using data from [Poling et al](https://personal.ems.psu.edu/~radovic/LennardJones_1.pdf).
| interaction | $\epsilon (eV)$ | $\sigma (Å)$ |
| ----------- | ------------------- | ------------ |
| $H_2$-$H_2$ | $6.6\times 10^{-4}$ |  2.918       |
| $N_2$-$N_2$ | $3.5\times 10^{-3}$ |  3.614       |
| $H_2$-$N_2$ | $1.5\times 10^{-3}$ |  3.266       |
|$NH_3$-$NH_3$| $4.8\times 10^{-2}$ |  2.900       |
|$H_2$-$NH_3$ | $5.6\times 10^{-3}$ |  2.909       |
|$N_2$-$NH_3$ | $1.3\times 10^{-2}$ |  3.257       |

Set $\mu_{NH_3}$ to -0.1 eV.<br>
![LJ comparisons](output_figures/lj_comparisons.png)

Resulting phase diagram:
![Ru H2-N2-NH3 binding](output_figures/ru_h2_n2_nh3_binding.png)

## 5: Animation
Animation with $\mu_{H2}=0$ and $T=116K$, saving every 100 frames for a total simulation time of 10000 frames:
![animation of competitive sorption](animations/ru_h2_n2_nh3_binding.gif)

<br>

[back](../README.md)
