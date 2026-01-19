import numpy as np 
from src.pdfs.pdf_set import *
from src.particles.pdg_id import pdg_id
from src.couplings.couplings import *
from src.energy_spectrum.muon_width import Muon_Width
from src.cuts.theory_cut import *

class Neutrino_NC_Cross_Section():
    def __init__(self,replica):
        self.replica = replica
    # nu u > nu u proton
    def dsigmah_dxdQ2_NC_nu_u(self, variables):    
        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        W2     = (mp**2 +Q2 * (1-x)/x) 
        
        
        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("u"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("u"), x, Q2)-pdfs[0].xfxQ2(pdg_id("u"), x, Q2)

        # Spectrum
        dGamma = Muon_Width.dGamma_dE_nu(E_nu)
        
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        sigma = Gfsqpi *  pdf/ x * (gZuLSQ + gZuRSQ * (1 - y)**2 )  * theory_cut * theory_limits /( (1+Q2/mZ**2)**2 )

        return dGamma * sigma

    # nu ub > nu ub proton
    def dsigmah_dxdQ2_NC_nu_ub(self, variables):    
        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        
        OneMinySQ = (1 - y)**2 # this is the factor (1-y)^2
        
        
        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("ub"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("ub"), x, Q2)-pdfs[0].xfxQ2(pdg_id("ub"), x, Q2)

        W2     = (mp**2 +Q2 * (1-x)/x) 
        
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        sigma = Gfsqpi * (pdf/ x * (gZuLSQ + gZuRSQ * OneMinySQ ) ) * theory_cut * theory_limits * (1/(1+Q2/mZ**2)**2)

        return Muon_Width.dGamma_dE_nu(E_nu) * sigma



    # nu u > nu u (neutron)
    def dsigmah_dxdQ2_NC_nu_u_neutron(self, variables):    
        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        
        OneMinySQ = (1 - y)**2 # this is the factor (1-y)^2


        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("d"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("d"), x, Q2)-pdfs[0].xfxQ2(pdg_id("d"), x, Q2)

        W2     = (mp**2 +Q2 * (1-x)/x) 
        
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        sigma = Gfsqpi * (pdf/ x * (gZuLSQ + gZuRSQ * OneMinySQ ) ) * theory_cut * theory_limits * (1/(1+Q2/mZ**2)**2)

        return Muon_Width.dGamma_dE_nu(E_nu) * sigma
    
    # nu ub > nu ub (neutron)
    def dsigmah_dxdQ2_NC_nu_ub_neutron(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2
        
        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("db"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("db"), x, Q2)-pdfs[0].xfxQ2(pdg_id("db"), x, Q2)

        W2     = (mp**2 +Q2 * (1-x)/x) 
        
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)
        
        sigma = Gfsqpi * (pdf/ x * (gZuRSQ + gZuLSQ * OneMinySQ ) ) * theory_cut * theory_limits * (1/(1+Q2/mZ**2)**2)

        return Muon_Width.dGamma_dE_nu(E_nu) * sigma

    #-------------------------------------
    # nu d > nu d
    def dsigmah_dxdQ2_NC_nu_d(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2

        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2
        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("d"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("d"), x, Q2)-pdfs[0].xfxQ2(pdg_id("d"), x, Q2)
        
        W2    = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)
        
        sigma = Gfsqpi * (pdf / x * (gZdLSQ + gZdRSQ * OneMinySQ ) ) * theory_cut * theory_limits /( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nu(E_nu) * sigma

    # nu db > nu db
    def dsigmah_dxdQ2_NC_nu_db(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2

        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2
        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("db"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("db"), x, Q2)-pdfs[0].xfxQ2(pdg_id("db"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        sigma = Gfsqpi * (pdf / x * (gZdRSQ + gZdLSQ * OneMinySQ ) ) * theory_cut * theory_limits /( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nu(E_nu) * sigma

    # nu d > nu d (neutron)
    def dsigmah_dxdQ2_NC_nu_d_neutron(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2

        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2
        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("u"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("u"), x, Q2)-pdfs[0].xfxQ2(pdg_id("u"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        
        sigma = Gfsqpi * (pdf / x * (gZdLSQ + gZdRSQ * OneMinySQ ) ) * theory_cut * theory_limits /( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nu(E_nu) * sigma

    # nu db > nu db (neutron)
    def dsigmah_dxdQ2_NC_nu_db_neutron(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2

        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2
        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("ub"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("ub"), x, Q2)-pdfs[0].xfxQ2(pdg_id("ub"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        sigma = Gfsqpi * (pdf /x * (gZdRSQ + gZdLSQ * OneMinySQ ) ) * theory_cut * theory_limits /( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nu(E_nu) * sigma

    #-------------------------------------
    # nu s > nu s
    def dsigmah_dxdQ2_NC_nu_s(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2

        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2
        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("s"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("s"), x, Q2)-pdfs[0].xfxQ2(pdg_id("s"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        
        sigma = Gfsqpi * (pdf / x * (gZdLSQ + gZdRSQ * OneMinySQ ) ) * theory_cut * theory_limits /( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nu(E_nu) * sigma

    # nu sb > nu sb
    def dsigmah_dxdQ2_NC_nu_sb(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2

        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2
        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("sb"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("sb"), x, Q2)-pdfs[0].xfxQ2(pdg_id("sb"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        
        sigma = Gfsqpi * (pdf / x * (gZdRSQ + gZdLSQ * OneMinySQ ) ) * theory_cut * theory_limits /( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nu(E_nu) * sigma


    #-------------------------------------
    # nu c > nu c
    def dsigmah_dxdQ2_NC_nu_c(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2
        
        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("c"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("c"), x, Q2)-pdfs[0].xfxQ2(pdg_id("c"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        
        sigma = Gfsqpi * (pdf / x* (gZuLSQ + gZuRSQ * OneMinySQ ) ) * theory_cut * theory_limits /( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nu(E_nu) * sigma

    # nu cb > nu cb
    def dsigmah_dxdQ2_NC_nu_cb(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2
        
        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("cb"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("cb"), x, Q2)-pdfs[0].xfxQ2(pdg_id("cb"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        
        sigma = Gfsqpi * (pdf /x * (gZuRSQ + gZuLSQ * OneMinySQ ) ) * theory_cut * theory_limits /( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nu(E_nu) * sigma


    #-------------------------------------
    # nu b > nu b
    def dsigmah_dxdQ2_NC_nu_b(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2
        
        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("b"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("b"), x, Q2)-pdfs[0].xfxQ2(pdg_id("b"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)
        
        sigma = Gfsqpi * (pdf / x * (gZdLSQ + gZdRSQ * OneMinySQ ) ) * theory_cut * theory_limits /( (1+Q2/mZ**2)**2 ) 
        return Muon_Width.dGamma_dE_nu(E_nu) * sigma

    # nu bb > nu bb
    def dsigmah_dxdQ2_NC_nu_bb(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2
        
        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("bb"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("bb"), x, Q2)-pdfs[0].xfxQ2(pdg_id("bb"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        
        sigma = Gfsqpi * (pdf / x * (gZdRSQ + gZdLSQ * OneMinySQ ) )  * theory_cut * theory_limits /( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nu(E_nu) * sigma



    # #-------------------------------------
    # nub u > nub u
    def dsigmah_dxdQ2_NC_nub_u(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2
        
        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("u"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("u"), x, Q2)-pdfs[0].xfxQ2(pdg_id("u"), x, Q2)

        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)


        sigma = Gfsqpi * (pdf / x * (gZuRSQ + gZuLSQ * OneMinySQ ) ) * theory_cut * theory_limits /( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nubar(E_nu) * sigma

    # nub ub > nub ub
    def dsigmah_dxdQ2_NC_nub_ub(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2
        
        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("ub"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("ub"), x, Q2)-pdfs[0].xfxQ2(pdg_id("ub"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        
        sigma = Gfsqpi * (pdf / x * (gZuLSQ + gZuRSQ * OneMinySQ ) ) * theory_cut * theory_limits /( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nubar(E_nu) * sigma



    # nub u > nub u (neutron)
    def dsigmah_dxdQ2_NC_nub_u_neutron(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2
        
        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("d"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("d"), x, Q2)-pdfs[0].xfxQ2(pdg_id("d"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        
        sigma = Gfsqpi * (pdf / x * (gZuRSQ + gZuLSQ * OneMinySQ ) ) * theory_cut * theory_limits /( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nubar(E_nu) * sigma

    # nub ub > nub ub (neutron)
    def dsigmah_dxdQ2_NC_nub_ub_neutron(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2
        
        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("db"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("db"), x, Q2)-pdfs[0].xfxQ2(pdg_id("db"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)
        
        sigma = Gfsqpi * (pdf / x * (gZuLSQ + gZuRSQ * OneMinySQ ) ) * theory_cut * theory_limits /( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nubar(E_nu) * sigma



    #-------------------------------------
    # nub d > nub d
    def dsigmah_dxdQ2_NC_nub_d(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2
        
        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("d"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("d"), x, Q2)-pdfs[0].xfxQ2(pdg_id("d"), x, Q2)

        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        
        sigma = Gfsqpi * (pdf / x * (gZdRSQ + gZdLSQ * OneMinySQ ) ) * theory_cut * theory_limits /( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nubar(E_nu) * sigma

    # nub db > nub db
    def dsigmah_dxdQ2_NC_nub_db(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2
        
        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("db"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("db"), x, Q2)-pdfs[0].xfxQ2(pdg_id("db"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)
        
        sigma = Gfsqpi * (pdf / x * (gZdLSQ + gZdRSQ * OneMinySQ ) ) * theory_cut * theory_limits /( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nubar(E_nu) * sigma

    # nub d > nub d (neutron)
    def dsigmah_dxdQ2_NC_nub_d_neutron(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2
        
        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("u"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("u"), x, Q2)-pdfs[0].xfxQ2(pdg_id("u"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        sigma = Gfsqpi * (pdf / x  * (gZdRSQ + gZdLSQ * OneMinySQ ) ) * theory_cut * theory_limits / ( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nubar(E_nu) * sigma

    # nub db > nub db (neutron)
    def dsigmah_dxdQ2_NC_nub_db_neutron(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2
        
        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("ub"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("ub"), x, Q2)-pdfs[0].xfxQ2(pdg_id("ub"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        sigma = Gfsqpi * (pdf / x * (gZdLSQ + gZdRSQ * OneMinySQ ) ) * theory_cut * theory_limits / ( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nubar(E_nu) * sigma

    #-------------------------------------
    # nub s > nub s
    def dsigmah_dxdQ2_NC_nub_s(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2

        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("s"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("s"), x, Q2)-pdfs[0].xfxQ2(pdg_id("s"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)
        
        sigma = Gfsqpi * (pdf  / x * (gZdRSQ + gZdLSQ * OneMinySQ ) ) * theory_cut * theory_limits / ( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nubar(E_nu) * sigma

    # nub sb > nub sb
    def dsigmah_dxdQ2_NC_nub_sb(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2
        
        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("sb"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("sb"), x, Q2)-pdfs[0].xfxQ2(pdg_id("sb"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        sigma = Gfsqpi * (pdf / x * (gZdLSQ + gZdRSQ * OneMinySQ ) ) * theory_cut * theory_limits / ( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nubar(E_nu) * sigma


    #-------------------------------------
    # nub c > nub c
    def dsigmah_dxdQ2_NC_nub_c(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2
        
        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("c"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("c"), x, Q2)-pdfs[0].xfxQ2(pdg_id("c"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        sigma = Gfsqpi * (pdf / x * (gZuRSQ + gZuLSQ * OneMinySQ ) ) * theory_cut * theory_limits / ( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nubar(E_nu) * sigma

    # nub cb > nub cb
    def dsigmah_dxdQ2_NC_nub_cb(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2
        
        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("cb"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("cb"), x, Q2)-pdfs[0].xfxQ2(pdg_id("cb"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        sigma = Gfsqpi * (pdf / x * (gZuLSQ + gZuRSQ * OneMinySQ ) ) * theory_cut * theory_limits / ( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nubar(E_nu) * sigma


    #-------------------------------------
    # nub b > nub b
    def dsigmah_dxdQ2_NC_nub_b(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2
        
        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("b"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("b"), x, Q2)-pdfs[0].xfxQ2(pdg_id("b"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)
        
        sigma = Gfsqpi * (pdf / x * (gZdRSQ + gZdLSQ * OneMinySQ ) ) * theory_cut * theory_limits / ( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nubar(E_nu) * sigma

    # nub bb > nub bb
    def dsigmah_dxdQ2_NC_nub_bb(self, variables):    

        x, Q2, E_nu = variables
        y = Q2 / (x * 2 * E_nu * mp)
        s = 2 * E_nu * mp + mp**2
        
        OneMinySQ = (1 - Q2 / (x * s))**2 # this is the factor (1-y)^2

        if self.replica == 0:
            pdf = pdfs[self.replica].xfxQ2(pdg_id("bb"), x, Q2)
        else :
            pdf = pdfs[self.replica].xfxQ2(pdg_id("bb"), x, Q2)-pdfs[0].xfxQ2(pdg_id("bb"), x, Q2)
        
        W2   = (mp**2 +Q2 * (1-x)/x)
        theory_limits   = np.where(Q2 /(4* E_nu**2*(1-y))>1,0,1 )
        theory_cut      = np.where(W2 < W2_min, 0, 1) * np.where(Q2<Q2_min,0,1)*np.where(Q2>x*2*E_nu*mp,0,1)

        sigma = Gfsqpi * (pdf / x * (gZdLSQ + gZdRSQ * OneMinySQ ) ) *theory_cut * theory_limits / ( (1+Q2/mZ**2)**2 )

        return Muon_Width.dGamma_dE_nubar(E_nu) * sigma

    #-------------------------------------
    # nu nbar inclusive processes 

    def dsigma_MuC_dxdQ2_nunub_u(self, variables):
        return self.dsigmah_dxdQ2_NC_nu_u(variables) + self.dsigmah_dxdQ2_NC_nub_u(variables)
        
    def dsigma_MuC_dxdQ2_nunub_ub(self, variables):
        return self.dsigmah_dxdQ2_NC_nu_ub(variables) + self.dsigmah_dxdQ2_NC_nub_ub(variables)
        
    def dsigma_MuC_dxdQ2_nunub_d(self, variables):
        return self.dsigmah_dxdQ2_NC_nu_d(variables) + self.dsigmah_dxdQ2_NC_nub_d(variables)

    def dsigma_MuC_dxdQ2_nunub_db(self, variables):
        return self.dsigmah_dxdQ2_NC_nu_db(variables) + self.dsigmah_dxdQ2_NC_nub_db(variables)

    def dsigma_MuC_dxdQ2_nunub_s(self, variables):
        return self.dsigmah_dxdQ2_NC_nu_s(variables) + self.dsigmah_dxdQ2_NC_nub_s(variables)

    def dsigma_MuC_dxdQ2_nunub_sb(self, variables):
        return self.dsigmah_dxdQ2_NC_nu_sb(variables) + self.dsigmah_dxdQ2_NC_nub_sb(variables)

    def dsigma_MuC_dxdQ2_nunub_c(self, variables):
        return self.dsigmah_dxdQ2_NC_nu_c(variables) + self.dsigmah_dxdQ2_NC_nub_c(variables)

    def dsigma_MuC_dxdQ2_nunub_cb(self, variables):
        return self.dsigmah_dxdQ2_NC_nu_cb(variables) + self.dsigmah_dxdQ2_NC_nub_cb(variables)

    def dsigma_MuC_dxdQ2_nunub_b(self, variables):
        return self.dsigmah_dxdQ2_NC_nu_b(variables) + self.dsigmah_dxdQ2_NC_nub_b(variables)

    def dsigma_MuC_dxdQ2_nunub_bb(self, variables):
        return self.dsigmah_dxdQ2_NC_nu_bb(variables) + self.dsigmah_dxdQ2_NC_nub_bb(variables)

    def dsigma_MuC_dxdQ2_nunub_u_Nucleon(self, variables):
        return self.dsigmah_dxdQ2_NC_nu_u_Nucleon(variables) + self.dsigmah_dxdQ2_NC_nub_u_Nucleon(variables)

    def dsigma_MuC_dxdQ2_nunub_ub_Nucleon(self, variables):
        return self.dsigmah_dxdQ2_NC_nu_ub_Nucleon(variables) + self.dsigmah_dxdQ2_NC_nub_ub_Nucleon(variables)

    def dsigma_MuC_dxdQ2_nunub_d_Nucleon(self, variables):
        return self.dsigmah_dxdQ2_NC_nu_d_Nucleon(variables) + self.dsigmah_dxdQ2_NC_nub_d_Nucleon(variables)

    def dsigma_MuC_dxdQ2_nunub_db_Nucleon(self, variables):
        return self.dsigmah_dxdQ2_NC_nu_db_Nucleon(variables) + self.dsigmah_dxdQ2_NC_nub_db_Nucleon(variables)

    def dsigma_MuC_dxdQ2_nunub_u_neutron(self, variables):
        return self.dsigmah_dxdQ2_NC_nu_u_neutron(variables) + self.dsigmah_dxdQ2_NC_nub_u_neutron(variables)

    def dsigma_MuC_dxdQ2_nunub_ub_neutron(self, variables):
        return self.dsigmah_dxdQ2_NC_nu_ub_neutron(variables) + self.dsigmah_dxdQ2_NC_nub_ub_neutron(variables)

    def dsigma_MuC_dxdQ2_nunub_d_neutron(self, variables):
        return self.dsigmah_dxdQ2_NC_nu_d_neutron(variables) + self.dsigmah_dxdQ2_NC_nub_d_neutron(variables)

    def dsigma_MuC_dxdQ2_nunub_db_neutron(self, variables):
        return self.dsigmah_dxdQ2_NC_nu_db_neutron(variables) + self.dsigmah_dxdQ2_NC_nub_db_neutron(variables)
