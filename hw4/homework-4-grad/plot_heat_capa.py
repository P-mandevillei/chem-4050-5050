from lj_particles import CanonicalEnsemble, lj_partition
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

if __name__ == '__main__':
	length = 10 # angstrom
	n = 10
	temp_1 = 10
	temp_2 = 1000
	grid_n = 200
	
	# construct the ensemble
	ljParticles = CanonicalEnsemble(
	    lambda temps: lj_partition(temps, length=length, n=n),
	    temp_1,
	    temp_2
	)
	
	# calculate thermo properties
	ljParticles.calc_thermo(grid_n)

	# extract dissociation temperature
	df = ljParticles.thermo_df
	max_row = df['constant_volume_heat_capacity_J/K'].idxmax(axis=0)
	dissoc_temp = df.iloc[max_row]['temperature']

	# plot
	fig, ax = plt.subplots()
	fig.set_size_inches(6, 4)
	sns.lineplot(
	    ljParticles.thermo_df,
	    x = 'temperature',
	    y = 'constant_volume_heat_capacity_J/K'
	)
	ax.axvline(x=dissoc_temp, 
	           label=f'Dissociation Temperature = {dissoc_temp:.2f} K',
	           linestyle='--', color='r', alpha=0.6)
	ax.set_xlabel('Temperature (K)')
	ax.set_ylabel('$C_v$ (J/K)')
	ax.set_title("$C_v$ vs. temperature for a pair of Lennard-Jones Ar particles")
	ax.legend()
	fig.savefig("cv_vs_temperature.png", dpi=210)
	print("Saved plot to cv_vs_temperature.png")