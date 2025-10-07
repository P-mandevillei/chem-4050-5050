from ideal_gas import IdealGas
import numpy as np
import pandas as pd

class IsothermalIdealGas(IdealGas):
	"""
	Defines an ideal gas under isothermal conditions
	"""
	def pressure_at(self, vol): # override
		"""
		Calculates the pressure at the given volume.

		Parameters:
		-----------------------------
		vol (m^3): current volume.

		Returns:
		Pressure (Pa) at the current volume
		"""
		return self.n*self.R*self.temp_i/vol

def compute_work_isothermal(vol_i, vol_f, n=1, temp=298):
    """
	Computes the work done on an ideal gas undergoing isothermal expansion

	Parameters:
	-----------------------------
	vol_i (m^3): initial volume
	vol_f (m^3): final volume
	n (mol): number of moles
	temp (K): temperature

	Returns:
	Work done from vol_i to vol_f (J)
	"""
    isotherm = IsothermalIdealGas(n, temp, vol_i)
    return isotherm.work(vol_i, vol_f)

if __name__ == "__main__":
	# compute works
	n = 1 # mol
	temp_i = 300 # K
	v_i = 0.1 # m^3
	gamma = 1.4
	vols = np.linspace(v_i, 3*v_i, 100)
	works = pd.DataFrame({'final_volume': vols})
	works['isothermal'] = works['final_volume'].apply(lambda v_f: compute_work_isothermal(v_i, v_f, n, temp_i))
	works.to_csv('isothermal_work.csv')
	print("Saved computed works to isothermal_work.csv")