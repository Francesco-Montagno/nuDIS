import numpy as np 
from src.pdfs.pdf_set import *
from src.particles.pdg_id import pdg_id
from src.couplings.couplings import *
from src.energy_spectrum.muon_width import Muon_Width
from src.cuts.theory_cut import *



class Neutrino_Muon_Cross_Section():
    def __init__(self, replica):
        self.replica = replica 
        
    # Proton
    def sigma_d_Nu_p(self, variables) : 
        # Variables
        x, Q2, E_nu = variables   
        y           = Q2 / (x * 2 * E_nu * mp)
        W2          = (mp**2 +Q2 * (1-x)/x)
        
        # PDFs
        if self.replica == 0 :
            x_dx    =  p_central.xfxQ2(pdg_id("d"), x, Q2)
        else :
            x_dx = pdfs[self.replica].xfxQ2(pdg_id("d"), x, Q2) - p_central.xfxQ2(pdg_id("d"), x, Q2)
            
        # Spectrum
        dGamma  =   Muon_Width.dGamma_dE_nu(E_nu)
        
        # Cuts
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1) * np.where(Q2>x*2*E_nu*mp,0,1)
        
        # Cross Section
        xsec =   Gfsqpi    * x_dx / x * theory_cut * theory_limits /(  (1+ Q2 / mW**2)**2 )

        return xsec * dGamma

    def sigma_ubar_Nu_p(self, variables) : 
        # Variables
        x, Q2, E_nu = variables
        y           = Q2 / (x * 2 * E_nu * mp)
        W2          = (mp**2 +Q2 * (1-x)/x)

        # PDFs
        if self.replica == 0:
            x_ubarx     =   p_central.xfxQ2(pdg_id("ub"), x, Q2)
        else:
            x_ubarx = pdfs[self.replica].xfxQ2(pdg_id("ub"), x, Q2) - p_central.xfxQ2(pdg_id("ub"), x, Q2)

        # Spectrum
        dGamma      =   Muon_Width.dGamma_dE_nu(E_nu)

        # Cuts
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)
        
        # Cross Section
        xsec =   Gfsqpi  * x_ubarx/x * (1-y)**2 * theory_cut * theory_limits /(  (1+ Q2 / mW**2)**2 ) 

        return xsec * dGamma

    # Neutron (using Isospin)
    def sigma_d_Nu_n(self, variables) : 
        # Variables
        x, Q2, E_nu = variables   
        y           = Q2 / (x * 2 * E_nu * mp)
        W2          = (mp**2 +Q2 * (1-x)/x)
        
        # PDFs
        if self.replica == 0 :
            x_dx    =  p_central.xfxQ2(pdg_id("u"), x, Q2)
        else :
            x_dx    =  pdfs[self.replica].xfxQ2(pdg_id("u"), x, Q2) - p_central.xfxQ2(pdg_id("u"), x, Q2)

        # Spectrum
        dGamma  =   Muon_Width.dGamma_dE_nu(E_nu)

        # Cuts
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1) * np.where(Q2>x*2*E_nu*mp,0,1)
        
        # Cross Section
        xsec =   Gfsqpi    * x_dx / x * theory_cut * theory_limits /(  (1+ Q2 / mW**2)**2 )

        return xsec * dGamma

    def sigma_ubar_Nu_n(self, variables) : 
        # Variables
        x, Q2, E_nu = variables
        y           = Q2 / (x * 2 * E_nu * mp)
        W2          = (mp**2 +Q2 * (1-x)/x)

        # PDFs
        if self.replica == 0:
            x_ubarx     =   p_central.xfxQ2(pdg_id("db"), x, Q2)
        else:
            x_ubarx     =   pdfs[self.replica].xfxQ2(pdg_id("db"), x, Q2) - p_central.xfxQ2(pdg_id("db"), x, Q2)

        # Spectrum
        dGamma      =   Muon_Width.dGamma_dE_nu(E_nu)

        # Cuts
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)
        
        # Cross Section
        xsec =   Gfsqpi  * x_ubarx/x * (1-y)**2 * theory_cut * theory_limits /(  (1+ Q2 / mW**2)**2 ) 

        return xsec * dGamma

    # Sea Quarks
    def sigma_b_Nu_p(self, variables) : 
        # Variables
        x, Q2, E_nu = variables
        y           = Q2 / (x * 2 * E_nu * mp)
        W2          = (mp**2 +Q2 * (1-x)/x)
        
        # PDFs
        if self.replica == 0:
            x_bx        =   p_central.xfxQ2(pdg_id("b"), x, Q2)
        else:
            x_bx = pdfs[self.replica].xfxQ2(pdg_id("b"), x, Q2) - p_central.xfxQ2(pdg_id("b"), x, Q2)

        # Spectrum
        dGamma      =   Muon_Width.dGamma_dE_nu(E_nu)

        # Cuts
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)
        
        # Cross Section
        xsec =   Gfsqpi   * x_bx/x * theory_cut * theory_limits  /(  (1+ Q2 / mW**2)**2 )

        return xsec * dGamma

    def sigma_s_Nu_p(self, variables) : 
        # Variables
        x, Q2, E_nu = variables
        y           = Q2 / (x * 2 * E_nu * mp)
        W2          = (mp**2 +Q2 * (1-x)/x)
        
        # PDFs
        if self.replica == 0:
            x_sx        =   p_central.xfxQ2(pdg_id("s"), x, Q2)
        else :
            x_sx = pdfs[self.replica].xfxQ2(pdg_id("s"), x, Q2) - p_central.xfxQ2(pdg_id("s"), x, Q2)
        # Spectrum
        dGamma      =   Muon_Width.dGamma_dE_nu(E_nu)

        # Cuts
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)
        
        # Cross Section
        xsec =   Gfsqpi   * x_sx/x * theory_cut * theory_limits  /(  (1+ Q2 / mW**2)**2 )

        return xsec * dGamma

    def sigma_cbar_Nu_p(self, variables) : 
        # Variables
        x, Q2, E_nu = variables
        y           = Q2 / (x * 2 * E_nu * mp)
        W2          = (mp**2 +Q2 * (1-x)/x) 
        
        # PDFs
        if self.replica == 0:
            x_cbarx     =   p_central.xfxQ2(pdg_id("cb"), x, Q2)
        else:
            x_cbarx = pdfs[self.replica].xfxQ2(pdg_id("cb"), x, Q2) - p_central.xfxQ2(pdg_id("cb"), x, Q2)

        # Spectrum
        dGamma      =   Muon_Width.dGamma_dE_nu(E_nu)

        # Cuts
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)
        
        # Cross Section
        xsec =   Gfsqpi    * x_cbarx/x * (1-y)**2 * theory_cut * theory_limits  /(  (1+ Q2 / mW**2) **2 )
        
        return xsec * dGamma


