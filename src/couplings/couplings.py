
import numpy as np 
GeV = 1
GeVtopb = 0.3894e9
GeVtofb = 0.3894e12



m_mu    = 0.10566   * GeV
mp      = 0.938272  * GeV
mW      = 80.369    * GeV
mZ      = 91.1876   * GeV



GF      = 1.1663787e-5 * GeV**(-2)
Gfsqpi = GeVtopb * GF**2 / (np.pi)



pb_to_cm2 = 1e-36
N_Avogadro = 6.022e23

N_neutrino_per_year = 9*1e16

sWSQ = 0.23129  # +0.00004 :  weak mixing angle from PDG 2024 avarage

gZuL = 1/2 - (2/3) * sWSQ
gZuR = - (2/3) * sWSQ

gZdL = -1/2 + (1/3) * sWSQ
gZdR = (1/3) * sWSQ

gZuLSQ = gZuL**2
gZuRSQ = gZuR**2

gZdLSQ = gZdL**2
gZdRSQ = gZdR**2