import numpy as np 
from src.pdfs.pdf_set import *
from src.particles.pdg_id import pdg_id
from src.couplings.couplings import *
from src.energy_spectrum.muon_width import Muon_Width
from src.cuts.theory_cut import *

class AntiNeutrino_Electron_Cross_Section():
    def __init__(self, replica):
        self.replica = replica

    # Proton
    def sigma_dbar_NuBar_p(self, variables) : 
        
        x, Q2, E_nu = variables   
        y = Q2 / (x * 2 * E_nu * mp)
        
        if self.replica == 0 :
            x_dx = p_central.xfxQ2(pdg_id("db"), x, Q2)
        else :
            x_dx = pdfs[self.replica].xfxQ2(pdg_id("db"), x, Q2) - p_central.xfxQ2(pdg_id("db"), x, Q2)

        dGamma = Muon_Width.dGamma_dE_nubar(E_nu)

        W2     = (mp**2 +Q2 * (1-x)/x) 
        
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1) * np.where(Q2>x*2*E_nu*mp,0,1)
        
        xsec =   Gfsqpi    * x_dx / x * theory_cut * theory_limits /(  (1+ Q2 / mW**2)**2 )

        return xsec * dGamma

    def sigma_u_NuBar_p(self, variables) : 
        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)

        if self.replica == 0:
            x_ux     =   p_central.xfxQ2(pdg_id("u"), x, Q2)
        else:
            x_ux = pdfs[self.replica].xfxQ2(pdg_id("u"), x, Q2) - p_central.xfxQ2(pdg_id("u"), x, Q2)

        dGamma = Muon_Width.dGamma_dE_nubar(E_nu)

        W2     = (mp**2 +Q2 * (1-x)/x) 
        
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)
        
        xsec =   Gfsqpi  * x_ux/x * (1-y)**2 * theory_cut * theory_limits /(  (1+ Q2 / mW**2)**2 ) 

        return xsec * dGamma

    # Neutron (Using Isospin)
    def sigma_dbar_NuBar_n(self, variables) : 
        
        x, Q2, E_nu = variables   
        y = Q2 / (x * 2 * E_nu * mp)
        
        if self.replica == 0 :
            x_dx = p_central.xfxQ2(pdg_id("ub"), x, Q2)
        else :
            x_dx = pdfs[self.replica].xfxQ2(pdg_id("ub"), x, Q2) - p_central.xfxQ2(pdg_id("ub"), x, Q2)
        
        dGamma = Muon_Width.dGamma_dE_nubar(E_nu)

        W2     = (mp**2 +Q2 * (1-x)/x) 
        
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1) * np.where(Q2>x*2*E_nu*mp,0,1)
        
        xsec =   Gfsqpi    * x_dx / x * theory_cut * theory_limits /(  (1+ Q2 / mW**2)**2 )

        return xsec * dGamma

    def sigma_u_NuBar_n(self, variables) : 
        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)

        if self.replica == 0:
            x_ux     =   p_central.xfxQ2(pdg_id("d"), x, Q2)
        else:
            x_ux = pdfs[self.replica].xfxQ2(pdg_id("d"), x, Q2) - p_central.xfxQ2(pdg_id("d"), x, Q2)

        dGamma = Muon_Width.dGamma_dE_nubar(E_nu)

        W2     = (mp**2 +Q2 * (1-x)/x) 
        
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)
        
        xsec =   Gfsqpi  * x_ux/x * (1-y)**2 * theory_cut * theory_limits /(  (1+ Q2 / mW**2)**2 ) 

        return xsec * dGamma

    # Sea Quarks

    def sigma_sbar_NuBar_p(self, variables) : 
        x, Q2, E_nu = variables

        y = Q2 / (x * 2 * E_nu * mp)
        
        if self.replica == 0:
            x_sx        =   p_central.xfxQ2(pdg_id("sb"), x, Q2)
        else:
            x_sx = pdfs[self.replica].xfxQ2(pdg_id("sb"), x, Q2) - p_central.xfxQ2(pdg_id("sb"), x, Q2)

        dGamma = Muon_Width.dGamma_dE_nubar(E_nu)

        W2     = (mp**2 +Q2 * (1-x)/x) 
        
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)
        
        xsec =   Gfsqpi   * x_sx/x * theory_cut * theory_limits  /(  (1+ Q2 / mW**2)**2 )

        return xsec * dGamma

    def sigma_bbar_NuBar_p(self, variables) : 
        x, Q2, E_nu = variables

        y = Q2 / (x * 2 * E_nu * mp)
        if self.replica == 0:
            x_bx        =   p_central.xfxQ2(pdg_id("bb"), x, Q2)
        else:
            x_bx = pdfs[self.replica].xfxQ2(pdg_id("bb"), x, Q2) - p_central.xfxQ2(pdg_id("bb"), x, Q2)

        dGamma = Muon_Width.dGamma_dE_nubar(E_nu)

        W2     = (mp**2 +Q2 * (1-x)/x) 
        
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)
        
        xsec =   Gfsqpi   * x_bx/x * theory_cut * theory_limits  /(  (1+ Q2 / mW**2)**2 )

        return xsec * dGamma

    def sigma_c_NuBar_p(self, variables) : 
        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)

        if self.replica == 0:
            x_cx     =   p_central.xfxQ2(pdg_id("c"), x, Q2)
        else:
            x_cx = pdfs[self.replica].xfxQ2(pdg_id("c"), x, Q2) - p_central.xfxQ2(pdg_id("c"), x, Q2)

        dGamma = Muon_Width.dGamma_dE_nubar(E_nu)

        W2     = (mp**2 +Q2 * (1-x)/x) 
        
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        xsec =   Gfsqpi    * x_cx/x * (1-y)**2 * theory_cut * theory_limits  /(  (1+ Q2 / mW**2) **2 )

        return xsec * dGamma
