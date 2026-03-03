#--------------------------------------- Libraries ----------------------------------#
import vegas  
import numpy as np
import matplotlib.pyplot as plt
from src.couplings.couplings import * 
# ---------------------------------------------------------------------------
# 1. Configuring card path
# ---------------------------------------------------------------------------
import argparse
parser = argparse.ArgumentParser(description="Neutrino DIS Event Generator")
parser.add_argument('-c', '--card', type=str, default='card/events_card.dat', 
                    help='Path to the run card (default: card/events_card.dat)')
parser.add_argument('-f', '--folder', type=str, default=None,
                    help='Folder in which to save the output (overrides name from card)')
parser.add_argument('--n_events', type=int, nargs='?', default=2e6,
                    help='Number of events to generate (overrides n_events from card)')
parser.add_argument('--process', type=str, nargs='?', default=None,
                    help='Process to simulate (overrides process from card)')
parser.add_argument('--E_muon', type=float, nargs='?', default=None,
                    help='Muon energy (overrides E_muon from card)')
parser.add_argument('--replica', type=int, nargs='?', default=None,
                    help='PDF replica ID (overrides replica from card)')
parser.add_argument('--Q2_min_theory', type=float, nargs='?', default=None,
                    help='Minimum Q2 for theoretical cuts (overrides Q2_min_theory from card)')
parser.add_argument('--W2_min_theory', type=float, nargs='?', default=None,
                    help='Minimum W2 for theoretical cuts (overrides W2_min_theory from card)')
parser.add_argument('--vegas_iter', type=int, nargs='?', default=None,
                    help='Number of VEGAS iterations (overrides vegas_iter from card)')
parser.add_argument('--vegas_eval', type=int, nargs='?', default=None,
                    help='Number of VEGAS evaluations (overrides vegas_eval from card)')
parser.add_argument('--ncores', type=int, nargs='?', default=None,
                    help='Number of CPU cores to use (overrides ncores from card)')
parser.add_argument('--x_min_bin', type=float, nargs='?', default=None,
                    help='Minimum x for binning (overrides x_min_bin from card)')
parser.add_argument('--x_max_bin', type=float, nargs='?', default=None,
                    help='Maximum x for binning (overrides x_max_bin from card)')
parser.add_argument('--Q2_min_bin', type=float, nargs='?', default=None,
                    help='Minimum Q2 for binning (overrides Q2_min_bin from card)')
parser.add_argument('--Q2_max_bin', type=float, nargs='?', default=None,
                    help='Maximum Q2 for binning (overrides Q2_max_bin from card)')
parser.add_argument('--E_min_bin', type=float, nargs='?', default=None,
                    help='Minimum energy for binning (overrides E_min_bin from card)')
parser.add_argument('--E_max_bin', type=float, nargs='?', default=None,
                    help='Maximum energy for binning (overrides E_max_bin from card)')
args = parser.parse_args()

# ---------------------------------------------------------------------------
# 2. Loading and Reading the Card
# ---------------------------------------------------------------------------
from src.utils.read_card import * 
print(f"Reading configuration from: {args.card}")
card = InputCard(args.card)


# --- Name ---
if args.folder is None:
    name           = card.get("name")          # str: "My Neutrino DIS Run"
else:
    name          = args.folder         # str: "My Neutrino DIS Run"
# --- Process and Beam ---
if args.process is  None:
    process_name   = card.get("process")       # str: eg "cbar_p"
else:
    process_name   = args.process       # str: eg "cbar_p"
    
if args.E_muon is None:
    E_muon         = card.get("E_muon")        # float: 1500.0
else:
    E_muon         = args.E_muon        # float: 1500.0

    # --- PDF Settings ---
if args.replica is None:
    replica_id     = card.get("replica")       # int: 0
else:
    replica_id     = args.replica       # int: 0

    # --- Theoretical Cuts ---
if args.Q2_min_theory is None:
    Q2_min_theory  = card.get("Q2_min_theory") # float: 0.0
else:
    Q2_min_theory  = args.Q2_min_theory # float: 0.0

if args.W2_min_theory is None:
    W2_min_theory  = card.get("W2_min_theory") # float: 0.0
else:
    W2_min_theory  = args.W2_min_theory # float: 0.0

    # --- VEGAS Integration ---
if args.vegas_iter is None:
    n_iter         = card.get("vegas_iter")    # int: 10
else: 
    n_iter         = args.vegas_iter    # int: 10
if args.vegas_eval is None:
    n_eval         = card.get("vegas_eval")    # int: 1000
else:
    n_eval         = args.vegas_eval    # int: 1000
if args.ncores is None:
    ncores         = card.get("ncores")        # int: number of cores to use
else:
    ncores         = args.ncores        # int: number of cores to use
if args.n_events is None:
    n_events       = card.get("n_events")      # int: number of events to generate
else:
    n_events       = args.n_events      # int: number of events to generate
    # --- Binning: X ---
if args.x_min_bin is  None:
    x_min_bin      = card.get("x_min_bin")     # float: 0.001
else:
    x_min_bin      = args.x_min_bin     # float: 0.001
if args.x_max_bin is  None:
    x_max_bin      = card.get("x_max_bin")     # float: 1.0
else:
    x_max_bin      = args.x_max_bin     # float: 1.0

    # --- Binning: Q2 ---
if args.Q2_min_bin is  None:
    Q2_min_bin     = card.get("Q2_min_bin")    # float: 1.0
else:
    Q2_min_bin     = args.Q2_min_bin    # float: 1.0
if args.Q2_max_bin is  None:
    Q2_max_bin     = card.get("Q2_max_bin")    # float: 3000.0
else:
    Q2_max_bin     = args.Q2_max_bin    # float: 3000.0
    # --- Binning: Energy ---
if args.E_min_bin is  None:
    E_min_bin      = card.get("E_min_bin")     # float: 0.0
else: 
    E_min_bin      = args.E_min_bin     # float: 0.0
if args.E_max_bin is  None:
    E_max_bin      = card.get("E_max_bin")     # float: -1.0 to be handled
else:
    E_max_bin      = args.E_max_bin     # float: -1.0 to be handled

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

import src.cuts.theory_cut as cuts
cuts.configure(card)
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

PID_Process = {
        # nu p
        'd_p'   : (pdg_id("d"), 2212),
        's_p'   : (pdg_id("s"), 2212),
        'b_p'   : (pdg_id("b"), 2212),
        'ubar_p': (pdg_id("ub"), 2212),
        'cbar_p': (pdg_id("cb"), 2212),
        #--------------------------#
        # nubar p
        'dbar_p'  : (pdg_id("db"), 2212),
        'sbar_p'  : (pdg_id("sb"), 2212),
        'bbar_p'  : (pdg_id("bb"), 2212),
        'u_p'   : (pdg_id("u"), 2212),
        'c_p'   : (pdg_id("c"), 2212),
        #-------------------------#
        # nu n (neutron)
        'd_n'   : (pdg_id("d"), 2112),
        'ubar_n': (pdg_id("ub"), 2112),
        #--------------------------#
        # nubar n (neutron)
        'dbar_n'  : (pdg_id("db"), 2112),
        'u_n'   : (pdg_id("u"), 2112),
        #-------------------------#
        #---- Neutral Current ----#
        #--------------------------------#
        #----------- nu -----------------#
        'nu_u'   : (pdg_id("u"), 2212),
        'nu_ub'  : (pdg_id("ub"), 2212),
        'nu_d'   : (pdg_id("d"), 2212),
        'nu_db'  : (pdg_id("db"), 2212),
        'nu_s'   : (pdg_id("s"), 2212),
        'nu_sb'  : (pdg_id("sb"), 2212),
        'nu_c'   : (pdg_id("c"), 2212),
        'nu_cb'  : (pdg_id("cb"), 2212),
        'nu_b'   : (pdg_id("b"), 2212),
        'nu_bb'  : (pdg_id("bb"), 2212),
        #-------------------------------#
        #----------- nubar -------------#
        'nub_u'   : (pdg_id("ub"), 2212),
        'nub_ub'  : (pdg_id("ub"), 2212),
        'nub_d'   : (pdg_id("db"), 2212),
        'nub_db'  : (pdg_id("db"), 2212),
        'nub_s'   : (pdg_id("sb"), 2212),
        'nub_sb'  : (pdg_id("sb"), 2212),
        'nub_c'   : (pdg_id("cb"), 2212),
        'nub_cb'  : (pdg_id("cb"), 2212),
        'nub_b'   : (pdg_id("b"), 2212),
        'nub_bb'  : (pdg_id("bb"), 2212),
        #---------------------------------#
        #------------ nu-nubar -----------#
        'nunub_u' : (pdg_id("u"), 2212),
        'nunub_ub': (pdg_id("ub"), 2212),
        'nunub_d' : (pdg_id("d"), 2212),
        'nunub_db': (pdg_id("db"), 2212),
        'nunub_s' : (pdg_id("s"), 2212),
        'nunub_sb': (pdg_id("sb"), 2212),
        'nunub_c' : (pdg_id("c"), 2212),
        'nunub_cb': (pdg_id("cb"), 2212),
        'nunub_b' : (pdg_id("b"), 2212),
        'nunub_bb': (pdg_id("bb"), 2212),

        'nunub_u_neutron' : (pdg_id("u"), 2112),
        'nunub_ub_neutron': (pdg_id("ub"), 2112),
        'nunub_d_neutron' : (pdg_id("d"), 2112),
        'nunub_db_neutron': (pdg_id("db"), 2112),
}


if __name__ == '__main__':
    process = process_name
    from datetime import datetime

    if args.folder is not None:
        folder = Path("data","output",args.folder)
    else:
        if name == "":
            folder = Path("data","output") / datetime.now().strftime("%Y%m%d_%H%M%S")
        else:
            folder = Path("data","output") / name

    folder.mkdir(parents=True, exist_ok=True)
    print(f"Created folder: {folder}")
    
    integral = vegas.Integrator([[x_min_bin, x_max_bin], [Q2_min_bin, Q2_max_bin], [E_min_bin,E_max_bin]])

    # First do adaptive integration to get a good grid

    sigma = integral(crossSection[process], nitn=iterations, neval=evaluations, adapt = True,adapt_to_errors=False, nproc = ncores)
    # Now do the final integration with fixed grid
    sigma = integral(crossSection[process], nitn=2, neval=n_events, adapt = False,adapt_to_errors=False, nproc = ncores)
    print(sigma.summary())


#--------------------------- Generating and Saving Events -----------------------------#

#--------------------------- Generating and Saving Events -----------------------------#
    x_ev, Q2_ev, E_nu_ev, weight_ev = [], [], [], []

    for x, wgt in integral.random():
        f_val = crossSection[process]((x[0], x[1], x[2]))
        x_ev.append(x[0])
        Q2_ev.append(x[1])
        E_nu_ev.append(x[2])
        weight_ev.append(wgt * f_val)
        
    # Convertiamo in numpy array per i calcoli
    x_ev = np.array(x_ev)
    Q2_ev = np.array(Q2_ev)
    E_nu_ev = np.array(E_nu_ev)
    weight_ev = np.array(weight_ev)

    # ---------------------

    with open(args.card, 'r') as f_card:
        card_content = f_card.read()
        
    with open(f'{folder}/{process}.txt', 'w') as file:
        file.write("# =============================================================\n")
        file.write("#                       RUN CONFIGURATION                      \n")
        file.write("# =============================================================\n")
        for line in card_content.splitlines():
            file.write(f"# {line}\n")
        file.write("# =============================================================\n\n")
        print("# =============================================================\n\n")
        file.write("# PID Parton\t PID Nucleon\t x\tQ2\tE\tWeight\n")
        print("# PID Parton\t PID Nucleon\t x\tQ2\tE\tWeight\n")
        # Salviamo solo gli eventi con peso > 0 per pulizia
        mask = weight_ev > 0
        PID = PID_Process[process]
        for x, Q2, E, w in zip(x_ev[mask], Q2_ev[mask], E_nu_ev[mask], weight_ev[mask]):
            file.write(f"{PID[0]}\t{PID[1]}\t{x}\t{Q2}\t{E}\t{w}\n")
            print(f"{PID[0]}\t{PID[1]}\t{x}\t{Q2}\t{E}\t{w}")


















