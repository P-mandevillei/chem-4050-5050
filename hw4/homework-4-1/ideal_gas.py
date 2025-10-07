class IdealGas():
    """
    Base class for defining an ideal gas and computing its pressure-volume work
    """
    def __init__(self, n, temp_i, vol_i):
        """
        IdealGas(n, temp_i, vol_i)
        
        Parameters:
        ----------------------------
        n (mol): number of moles of ideal gas
        temp_i (K): initial temperature
        vol_i (m^3): initial volume
        """
        self.n = n
        self.temp_i = temp_i
        self.vol_i = vol_i
        self.R = 8.314 # J/mol-K
    
    def pressure_at(self, **kwargs):
        """
        Calculates the pressure with the given conditions
        
        Parameters:
        -----------------------------
        **kwargs: default calculation assumes an isotherm and needs argument `vol`. Default value is 1 m^3.
        
        Returns:
        Pressure (Pa) at the given condition
        """
        return self.n*self.R*self.temp_i/kwargs.get('vol', 1)
        
    def work(self, vol_i, vol_f):
        """
        Computes the work done from vol_i to vol_f
        
        Parameters:
        -----------------------------
        vol_i (m^3): initial volume
        vol_f (m^3): final volume
        
        Returns:
        Work done from vol_i to vol_f (J)
        """
        from scipy.integrate import trapezoid
        import numpy as np
        x = np.linspace(vol_i, vol_f, 1000)
        y = np.array([self.pressure_at(vol) for vol in x])
        work = -trapezoid(y, x)
        return work