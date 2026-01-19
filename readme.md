# A Leading Order $\nu\text{DIS}$ Monte Carlo Event Generator

**nuDIS_Generator** is a Monte Carlo event generator designed for Neutrino-Nucleon Deep Inelastic Scattering (DIS). It simulates neutrino fluxes originating from muon decays in flight, integrating physical cross-sections with LHAPDF sets and the VEGAS integration algorithm.

---

## 🛠 Installation and Setup

To ensure all  dependencies (like LHAPDF) are correctly installed, follow these steps:

### 1. Create the Conda Environment
Use the provided `environment.yml` file to create a dedicated environment with the necessary C++ and Python libraries:
```bash
conda env create -f environment.yml
```

### 2. Activate the Environment
```bash
conda activate nuDIS
```

### 3. Install Project Dependencies
Once the environment is active, install the project in "editable" mode. This will read the \`pyproject.toml\` file, install any remaining Python dependencies, and allow you to run the generator from anywhere:
```bash
pip install -e .
```

### 4. Download PDF Sets
You must download the LHAPDF sets specified in your configuration (e.g., NNPDF3.1):
```bash
lhapdf install PDF4LHC21_40
```

Notice that you can replace `PDF4LHC21_40` with any other LHAPDF set of your choice. You have to update the `pdf_set_name` parameter in `src/pdfs/pdf_set.py` accordingly.

---

## 🔬 Possible Processes

The following processes can be specified in the `process` field of the run card:

### Neutrino on Proton (mu- in final state)
* **d_p**: d + neutrino → mu- + h
* **s_p**: s + neutrino → mu- + h
* **b_p**: b + neutrino → mu- + h
* **ubar_p**: u~ + neutrino → mu- + h
* **cbar_p**: c~ + neutrino → mu- + h

### Antineutrino on Proton (mu+ in final state)
* **dbar_p**: d~ + antineutrino → mu+ + h
* **sbar_p**: s~ + antineutrino → mu+ + h
* **bbar_p**: b~ + antineutrino → mu+ + h
* **u_p**: u + antineutrino → mu+ + h
* **c_p**: c + antineutrino → mu+ + h

### Neutrino on Neutron (mu- in final state)
* **d_n**: d + neutrino → mu- + h
* **ubar_n**: u~ + neutrino → mu- + h

### Antineutrino on Neutron (mu+ in final state)
* **dbar_n**: d~ + antineutrino → mu+ + h
* **u_n**: u + antineutrino → mu+ + h

---

## ⚙️ Run Card Configuration

The simulation is controlled via `card/run_card.dat`. Below is a description of the parameters:

### Process and Beam
* **process**: The name of the process from the list above.
* **E_muon**: Energy of the parent muon beam in GeV.

### PDF Settings
* **replica**: The PDF replica ID. Use **0** for the central member, or N for specific error variations.

### Theoretical Cuts
* **Q2_min_theory**: Minimum momentum transfer squared ($Q^2$) allowed in the integration [GeV²].
* **W2_min_theory**: Minimum hadronic invariant mass squared ($W^2$) allowed [GeV²].

### VEGAS Parameters
* **vegas_iter**: Number of adaptive iterations for the VEGAS algorithm.
* **vegas_eval**: Number of function evaluations per iteration.

### Binning Configuration
* **num_bins_x / Q2 / E**: Number of bins for each observable.
* **x_min_bin / x_max_bin**: Range for Bjorken x.
* **Q2_min_bin / Q2_max_bin**: Range for $Q^2$ [GeV²].
* **E_min_bin / E_max_bin**: Range for the Neutrino Energy [GeV].
  * *Note*: If **E_max_bin** is set to **-1**, the generator automatically calculates the maximum kinematic limit based on the muon beam energy.

---

## 💻 Usage

Run the generator by passing the path to your run card:

```bash
python generator.py --card card/run_card.dat
```