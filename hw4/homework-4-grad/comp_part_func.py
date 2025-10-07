from lj_particles import CanonicalEnsemble, lj_partition
import pandas as pd
import numpy as np

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
	
	# calculate the partition function
	ljParticles.calc_z(grid_n)
	
	# export as csv
	df = pd.DataFrame({
		'temperature_K': np.linspace(temp_1, temp_2, grid_n),
		'z': ljParticles.z_vals
	})
	df.to_csv('partitions.csv', index=False)
	print("Saved partition functions to partitions.csv")