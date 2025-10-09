import numpy as np
from tqdm import tqdm
from scipy import stats

# Define the 1s orbital and its laplacian
def psi_1s(x, y, z, Z=1, a0=1):
    """
    1s wavefunction.

    Parameters
    -----------------------------
    x, y, z (atomic units): floats or 1d np.arrays of coordinates
    Z: atomic number
    a0 (atomic unit): Bohr radius

    Returns
    Float of 1d ap.array of 1s wavefunctions
    """
    a = a0/Z
    r = np.sqrt(x**2 + y**2 + z**2)
    coeff = 1/np.sqrt(np.pi*(a**3))
    radial = np.e ** (-r/a)
    return coeff * radial

def laplacian_psi_1s(x, y, z, Z=1, a0=1):
    """
    Laplacian of 1s orbital.

    Parameters
    -----------------------------
    x, y, z (atomic units): floats or 1d np.arrays of coordinates
    Z: atomic number
    a0 (atomic unit): Bohr radius

    Returns
    Float of 1d ap.array of Laplacian of 1s wavefunctions
    """
    a = a0/Z
    r = np.sqrt(x**2 + y**2 + z**2)
	# prevent divide by 0
    if type(r) == 'numpy.float64' and r == 0:
        r = 1e-12
    else:
        r[r == 0] = 1e-12
    wavefunctions = psi_1s(x, y, z, Z=Z, a0=a0)
    term_1 = (1/a0)**2 * wavefunctions
    term_2 = 2/r * (-1/a0)*wavefunctions
    return term_1 + term_2

# define the interface for how to draw from a distribution in Monte Carlo
class Drawer:
    """
    Monte Carlo sampling
    """
    def __init__(self, distribution, seed=114514):
        """
        Drawer(distribution, seed)

        Parameters
        ----------------------------
        distribution: a distribution that implements .draw(n, rng) and .pdf(points) where points are the return value of .draw
        seed: random seed
        """
        self.distribution = distribution
        self.seed = seed
        self.rng = np.random.default_rng(seed)

    def reseed(self):
        """
        Reseeds the random number generator instance.
        """
        self.rng = np.random.default_rng(self.seed)

    def draw(self, n):
        """
        Draws from the distribution
        """
        return self.distribution.draw(n, rng=self.rng)

    def pdf(self, points):
        """
        Computes the pdf values using input points.

        Parameters
        ------------------------------
        points: should be the returns of self.draw()
        """
        return self.distribution.pdf(points)

# Monte Carlo class with dependency injection
class MonteCarlo:
    """
    Monte Carlo integration.
    """

    def __init__(self, integrand, drawer):
        """
        MonteCarlo(integrand, drawer)

        Parameters
        -----------------------------------------
        integrand: a function that computes the integrand value. Must have a signature of integrand(points) where points are the output of drawer.draw()
        drawer: a Drawer instance
        """
        self.integrand = integrand
        self.drawer = drawer

    def integrate(self, n, chunk_size=int(1e5)):
        """
        Integrate using samples from the drawer instance and the integrand function.

        Parameters
        -----------------------------------------
        n: number of estimators
        chunk_size: the size of chunks to integrate at the same time. Prevents memory overflow

        Returns:
        Result of integration
        """
        print(f"Monte Carlo on {n} estimators")
        self.drawer.reseed()
        n_chunks = n // chunk_size
        remainder = n % chunk_size
        computed = []

        # sample by chunks
        for _ in tqdm(range(n_chunks)):
            samples = self.drawer.draw(chunk_size)
            integrands = self.integrand(samples) / self.drawer.pdf(samples)
            computed.append(integrands)
        samples = self.drawer.draw(remainder)
        integrands = self.integrand(samples) / self.drawer.pdf(samples)
        computed.append(integrands)
        # combine and average
        computed = np.concat(computed)
        integral = np.mean(computed)
        return integral

# uniform distribution using the adapter pattern for use in the Drawer class
class UniformDistribution:
    """
    Adapter class for np.random.uniform for use in the Drawer instance.
    """
    def __init__(self, L=7):
        """
        UniformDistribution(L=7)

        Parameters
        --------------------------------
        L: half side length of the cubix box over which to draw samples
        """
        self.L = L

    def draw(self, n, rng=None):
        if not rng:
            return None
        return rng.uniform(low=-self.L, high=self.L, size=(3, n))

    def pdf(self, points):
        return 1/((2*self.L)**3)

# Gaussian distribution using the adapter pattern for use in the Drawer class
class GaussianDistribution:
    """
    Adapter class for np.random.normal for use in the Drawer instance.
    """
    def __init__(self, loc=0.0, scale=1.0):
        """
        GaussianDistribution(loc=0.0, scale=1.0)

        Parameters
        ----------------------------------
        loc: center of the Gaussian
        scale: standard deviation of the Gaussian
        """
        self.loc = loc
        self.scale = scale

    def draw(self, n, rng=None):
        if not rng:
            return None
        return rng.normal(loc=self.loc, scale=self.scale, size=(3, n))

    def pdf(self, points):
        x, y, z = points
        x_pdf = stats.norm.pdf(x)
        y_pdf = stats.norm.pdf(y)
        z_pdf = stats.norm.pdf(z)
        return x_pdf * y_pdf * z_pdf

def diagonalKIntegrand(points):
    """
    Adapter class for calculating the integrand for the diagonal kinetic energy matrix elements.

    Parameters
    ----------------------------------
    points: size (3, n), coordinates x, y, z produced by a Drawer instance
    """
    x, y, z = points
    term_1 = psi_1s(x, y, z)
    term_2 = laplacian_psi_1s(x, y, z)
    return -1/2 * term_1 * term_2

def offDiagonalKIntegrand(Rz=1.4):
    """
    Adapter class for calculating the integrand for the off-diagonal kinetic energy matrix elements.

    Parameters
    ----------------------------------
    Rz (atomic units): separation along the z axis
    """

    # closure
    def calc_integrand(points):
        """
        Calculates the integrand for the off-diagonal kinetic energy matrix elements.
    
        Parameters
        ----------------------------------
        points: size (3, n), coordinates x, y, z produced by a Drawer instance
        """
        x, y, z = points
        term_1 = psi_1s(x, y, z+Rz/2)
        term_2 = laplacian_psi_1s(x, y, z-Rz/2)
        return -1/2 * term_1 * term_2
        
    return calc_integrand # wrapper function