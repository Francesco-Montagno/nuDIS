import argparse
parser = argparse.ArgumentParser(description="Neutrino DIS Event Generator")
parser.add_argument('-c', '--card', type=str, default='card/run_card.dat', 
                    help='Path to the run card (default: card/run_card.dat)')
args = parser.parse_args()

from src.utils.read_card import * 
print(f"Reading configuration from: {args.card}")
card = InputCard(args.card)


# --- Theoretical Cuts ---
Q2_min_theory  = card.get("Q2_min_theory") # float: 0.0
W2_min_theory  = card.get("W2_min_theory") # float: 0.0

W2_min = W2_min_theory 

Q2_min = Q2_min_theory