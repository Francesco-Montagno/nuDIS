import numpy as np

class Muon_Width:
    E_mu = None
    beta = None

    @classmethod
    def configure(cls, E_mu_val, beta_val):
        cls.E_mu = E_mu_val
        cls.beta = beta_val

    @classmethod
    def dGamma_dE_nu(cls, E_nu):    
        if cls.E_mu is None: raise RuntimeError("Muon_Width non configurato!")

        E_max = (1 + cls.beta) * cls.E_mu / 2

        dGamma_normalised = 1 / ( 3 * E_max  ) * (5 - 9 * (  E_nu / E_max ) **2 + 4 * (  E_nu / E_max )**3 ) 
        
        dGamma_normalised = dGamma_normalised * np.where((E_nu > 0) & (E_nu < E_max), 1, 0) 

        return dGamma_normalised 
    
    @classmethod
    def dGamma_dE_nubar(cls, E_nu_bar):
        if cls.E_mu is None: raise RuntimeError("Muon_Width non configurato!")

        E_max = (1 + cls.beta) * cls.E_mu / 2

        dGamma_normalised = 2 / ( E_max  ) * (1 - 3 * (  E_nu_bar / E_max ) **2 + 2 * (  E_nu_bar / E_max )**3 ) 

        dGamma_normalised = dGamma_normalised * np.where((E_nu_bar > 0) & (E_nu_bar < E_max), 1, 0) 

        return dGamma_normalised