
Q2_min = 0.0
W2_min = 0.0

def configure(card):
    global Q2_min, W2_min
    Q2_min = card.get("Q2_min_theory")
    W2_min = card.get("W2_min_theory")
    print(f"Cuts set to: Q2_min={Q2_min}, W2_min={W2_min}")