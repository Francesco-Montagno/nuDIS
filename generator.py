#--------------------------------------- Libraries ----------------------------------#
import vegas  
import numpy as np
from src.couplings.couplings import * 
# ---------------------------------------------------------------------------
# 1. Configuring card path
# ---------------------------------------------------------------------------
import argparse
parser = argparse.ArgumentParser(description="Neutrino DIS Event Generator")
parser.add_argument('-c', '--card', type=str, default='card/run_card.dat', 
                    help='Path to the run card (default: card/run_card.dat)')
args = parser.parse_args()

# ---------------------------------------------------------------------------
# 2. Loading and Reading the Card
# ---------------------------------------------------------------------------
from src.utils.read_card import * 
print(f"Reading configuration from: {args.card}")
card = InputCard(args.card)

# --- Process and Beam ---
process_name   = card.get("process")       # str: eg "cbar_p"
E_muon         = card.get("E_muon")        # float: 1500.0

# --- PDF Settings ---
replica_id     = card.get("replica")       # int: 0

# --- Theoretical Cuts ---
Q2_min_theory  = card.get("Q2_min_theory") # float: 0.0
W2_min_theory  = card.get("W2_min_theory") # float: 0.0

# --- VEGAS Integration ---
n_iter         = card.get("vegas_iter")    # int: 10
n_eval         = card.get("vegas_eval")    # int: 1000

# --- Binning: X ---
num_bins_x     = card.get("num_bins_x")    # int: 1
x_min_bin      = card.get("x_min_bin")     # float: 0.001
x_max_bin      = card.get("x_max_bin")     # float: 1.0

# --- Binning: Q2 ---
num_bins_Q2    = card.get("num_bins_Q2")   # int: 1
Q2_min_bin     = card.get("Q2_min_bin")    # float: 1.0
Q2_max_bin     = card.get("Q2_max_bin")    # float: 3000.0

# --- Binning: Energy ---
num_bins_E     = card.get("num_bins_E")    # int: 100
E_min_bin      = card.get("E_min_bin")     # float: 0.0
E_max_bin      = card.get("E_max_bin")     # float: -1.0 to be handled


#------------------------------------- Config ----------------------------------#

iterations  = n_iter
evaluations = n_eval

E_mu    = E_muon    
η       =  np.log(E_mu/m_mu * (1+ np.sqrt(1- (m_mu/E_mu)**2)))             # muon rapidity
β       = np.sqrt(1- (m_mu/E_mu)**2)                                       # muon velocity

Emax = (1+β)*E_mu/2
E_min_bin    = E_min_bin         

if E_max_bin == -1:
    E_max_bin    = Emax      

from src.energy_spectrum.muon_width import * 
Muon_Width.configure(E_mu, β)

from src.cuts.theory_cut import *
#---------------------------------- Cross Section Dictionary -------------------------------#
from src.cross_sections.cc.nu_e.antineutrino_positron import *
from src.cross_sections.cc.nu_mu.neutrino_muon import *
from src.cross_sections.nc.neutral_current_xsec import *
from datetime import datetime
from pathlib import Path

nc = Neutrino_NC_Cross_Section(replica_id)
cc_mu = Neutrino_Muon_Cross_Section(replica_id)
cc_e = AntiNeutrino_Electron_Cross_Section(replica_id)

crossSection = {
        # nu p
        'd_p'   : cc_mu.sigma_d_Nu_p,
        's_p'   : cc_mu.sigma_s_Nu_p,
        'b_p'   : cc_mu.sigma_b_Nu_p,
        'ubar_p': cc_mu.sigma_ubar_Nu_p,
        'cbar_p': cc_mu.sigma_cbar_Nu_p,
        #--------------------------#
        # nubar p
        'dbar_p'  : cc_e.sigma_dbar_NuBar_p,
        'sbar_p'  : cc_e.sigma_sbar_NuBar_p,
        'bbar_p'  : cc_e.sigma_bbar_NuBar_p,
        'u_p'   : cc_e.sigma_u_NuBar_p,
        'c_p'   : cc_e.sigma_c_NuBar_p,
        #-------------------------#
        # nu n (neutron)
        'd_n'   : cc_mu.sigma_d_Nu_n,
        'ubar_n': cc_mu.sigma_ubar_Nu_n,
        #--------------------------#
        # nubar n (neutron)
        'dbar_n'  : cc_e.sigma_dbar_NuBar_n,
        'u_n'   : cc_e.sigma_u_NuBar_n,
        #-------------------------#
        #---- Neutral Current ----#
        #--------------------------------#
        #----------- nu -----------------#
        'nu_u'   : nc.dsigmah_dxdQ2_NC_nu_u,
        'nu_ub'  : nc.dsigmah_dxdQ2_NC_nu_ub,
        'nu_d'   : nc.dsigmah_dxdQ2_NC_nu_d,
        'nu_db'  : nc.dsigmah_dxdQ2_NC_nu_db,
        'nu_s'   : nc.dsigmah_dxdQ2_NC_nu_s,
        'nu_sb'  : nc.dsigmah_dxdQ2_NC_nu_sb,
        'nu_c'   : nc.dsigmah_dxdQ2_NC_nu_c,
        'nu_cb'  : nc.dsigmah_dxdQ2_NC_nu_cb,
        'nu_b'   : nc.dsigmah_dxdQ2_NC_nu_b,
        'nu_bb'  : nc.dsigmah_dxdQ2_NC_nu_bb,
        #-------------------------------#
        #----------- nubar -------------#
        'nub_u'   : nc.dsigmah_dxdQ2_NC_nub_u,
        'nub_ub'  : nc.dsigmah_dxdQ2_NC_nub_ub,
        'nub_d'   : nc.dsigmah_dxdQ2_NC_nub_d,
        'nub_db'  : nc.dsigmah_dxdQ2_NC_nub_db,
        'nub_s'   : nc.dsigmah_dxdQ2_NC_nub_s,
        'nub_sb'  : nc.dsigmah_dxdQ2_NC_nub_sb,
        'nub_c'   : nc.dsigmah_dxdQ2_NC_nub_c,
        'nub_cb'  : nc.dsigmah_dxdQ2_NC_nub_cb,
        'nub_b'   : nc.dsigmah_dxdQ2_NC_nub_b,
        'nub_bb'  : nc.dsigmah_dxdQ2_NC_nub_bb,
        #---------------------------------#
        #------------ nu-nubar -----------#
        'nunub_u' : nc.dsigma_MuC_dxdQ2_nunub_u,
        'nunub_ub': nc.dsigma_MuC_dxdQ2_nunub_ub,
        'nunub_d' : nc.dsigma_MuC_dxdQ2_nunub_d,
        'nunub_db': nc.dsigma_MuC_dxdQ2_nunub_db,
        'nunub_s' : nc.dsigma_MuC_dxdQ2_nunub_s,
        'nunub_sb': nc.dsigma_MuC_dxdQ2_nunub_sb,
        'nunub_c' : nc.dsigma_MuC_dxdQ2_nunub_c,
        'nunub_cb': nc.dsigma_MuC_dxdQ2_nunub_cb,
        'nunub_b' : nc.dsigma_MuC_dxdQ2_nunub_b,
        'nunub_bb': nc.dsigma_MuC_dxdQ2_nunub_bb,
        
        'nunub_u_neutron' : nc.dsigma_MuC_dxdQ2_nunub_u_neutron,
        'nunub_ub_neutron': nc.dsigma_MuC_dxdQ2_nunub_ub_neutron,
        'nunub_d_neutron' : nc.dsigma_MuC_dxdQ2_nunub_d_neutron,
        'nunub_db_neutron': nc.dsigma_MuC_dxdQ2_nunub_db_neutron,
}

process = process_name
#----------------------------------- Binning ------------------------------------#
x_bins  = np.logspace(np.log10(x_min_bin), np.log10(x_max_bin), num_bins_x + 1)      # Logarithmic bins for x
Q2_bins = np.logspace(np.log10(Q2_min_bin), np.log10(Q2_max_bin), num_bins_Q2 + 1)  # Logarithmic bins for y
E_bins  = np.linspace(E_min_bin, E_max_bin, num_bins_E + 1)                          # Linear bins for E

#----------------------------------- Integration ------------------------------------#

sigma_bin = np.zeros((num_bins_x, num_bins_Q2, num_bins_E))
sigma_bin_error = np.zeros((num_bins_x, num_bins_Q2, num_bins_E))

def compute_bin(m):
    i, j, k = m  
    integral = vegas.Integrator([[x_bins[i], x_bins[i + 1]], [Q2_bins[j], Q2_bins[j + 1]], [E_bins[k], E_bins[k + 1]]])
    
    sigma = integral(crossSection[process], nitn=iterations, neval=evaluations, adapt = True,adapt_to_errors=False)
    print(f"Bin {i}, {j}, {k}: integral = {sigma.mean} ± {sigma.sdev}")
    return (i, j, k, 
            sigma.mean, sigma.sdev,)

bin_indices = [(i, j, k) for i in range(num_bins_x) for j in range(num_bins_Q2) for k in range(num_bins_E)]


results = [compute_bin(bin) for bin in bin_indices]

folder = Path("data","output") / datetime.now().strftime("%Y%m%d_%H%M%S")
folder.mkdir(parents=True, exist_ok=True)
print(f"Created folder: {folder}")

with open(f'{folder}/{process}.txt', 'w') as f:
    for result in results:
        i, j, k, xsec, xsec_dev = result
        f.write(f"{i}\t{j}\t{k}\t{xsec}\n")
        




















