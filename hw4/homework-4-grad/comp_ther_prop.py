from lj_particles import CanonicalEnsemble, lj_partition

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
	
	# calculate the thermo properties
	ljParticles.calc_thermo(grid_n)
	
	# export as csv
	ljParticles.thermo_df.rename(columns={'temperature': "temperature_K"}).to_csv("thermo_properties.csv", index=False)
	print("Saved thermo properties to thermo_properties.csv")