from ideal_gas import IdealGas
import numpy as np
import pandas as pd

class AdiabaticIdealGas(IdealGas):
    """ 
    Defines an ideal gas under adiabatic conditions 
    """
    def __init__(self, n, temp_i, vol_i, gamma): # override 
        """ 
        AdiabaticIdealGas(n, temp_i, vol_i, gamma) 
        Parameters:
        ---------------------------- 
        n (mol): number of moles of ideal gas
        temp_i (K): initial temperature 
        vol_i (m^3): initial volume 
        gamma: adiabatic index of the ideal gas 
        """ 
        super().__init__(n, temp_i, vol_i) 
        self.gamma = gamma 
        p_i = n*self.R*temp_i/vol_i # compute initial pressure 
        self.constant = p_i * vol_i**gamma 
    
    def pressure_at(self, vol): # override 
        """ 
        Calculates the pressure at the given volume. 
        
        Parameters: 
        ----------------------------- 
        vol (m^3): current volume. 
        
        Returns: 
        Pressure (Pa) at the current volume 
        """ 
        return self.constant / (vol**self.gamma)

def compute_work_adiabatic(vol_i, vol_f, gamma, n=1, temp_i=298):
    """
	Computes the work done on an ideal gas undergoing adiabatic expansion

	Parameters:
	-----------------------------
	vol_i (m^3): initial volume
	vol_f (m^3): final volume
    gamma: the adiabatic index of the ideal gas
	n (mol): number of moles
	temp_i (K): initial temperature

	Returns:
	Work done from vol_i to vol_f (J)
	"""
    adiabat = AdiabaticIdealGas(n, temp_i, vol_i, gamma)
    return adiabat.work(vol_i, vol_f)

if __name__ == "__main__":
	# compute works
	n = 1 # mol
	temp_i = 300 # K
	v_i = 0.1 # m^3
	gamma = 1.4
	vols = np.linspace(v_i, 3*v_i, 100)
	works = pd.DataFrame({'final_volume_m^3': vols})
	works['adiabatic_J'] = works['final_volume_m^3'].apply(lambda v_f: compute_work_adiabatic(v_i, v_f, gamma, n, temp_i))
	works.to_csv('adiabatic_work.csv')
	print("Saved computed works to adiabatic_work.csv")