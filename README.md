# 💊 Clinical Drug Level Modeling App
### *Pharmacokinetic Simulation (PK/PD) & Therapeutic Drug Monitoring (TDM) Engine*

[![Python](https://img.shields.io/badge/Python-3.12%20%7C%203.11-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Streamlit](https://img.shields.io/badge/Streamlit-1.42+-FF4B4B?logo=streamlit&logoColor=white)](https://streamlit.io/)
[![Plotly](https://img.shields.io/badge/Plotly-5.20+-3F4F75?logo=plotly&logoColor=white)](https://plotly.com/python/)
[![Tests](https://img.shields.io/badge/Pytest-64%20passed-success?logo=pytest&logoColor=white)](tests/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An open-source clinical pharmacokinetic simulation platform built for healthtech engineers, pharmacologists, and clinical researchers. The system implements exact one-compartment pharmacokinetic models with oral absorption (Bateman function) and intravenous (IV) bolus administration, steady-state accumulation analysis, missed-dose impact simulations, and interactive Therapeutic Drug Monitoring (TDM) with patient variability bounds.

---

## 🔬 Scientific & Mathematical Foundation

The core modeling engine (`core/pk_engine.py`) models drug disposition using linear, time-invariant one-compartment pharmacokinetics.

### 1. Oral Absorption – Bateman Function
For extravascular (oral) administration with first-order absorption and first-order elimination:

$$C(t) = \frac{F \cdot D \cdot k_a}{V_d \, (k_a - k_e)} \left( e^{-k_e t} - e^{-k_a t} \right)$$

Where:
* $D$ = Administered dose ($\text{mg}$)
* $F$ = Systemic bioavailability fraction ($0 < F \le 1$)
* $k_a$ = Absorption rate constant ($\text{h}^{-1}$)
* $k_e$ = Elimination rate constant ($\text{h}^{-1}$), derived from terminal half-life: $k_e = \frac{\ln 2}{t_{1/2}}$
* $V_d$ = Volume of distribution ($\text{L}$ or $\text{L/kg} \times \text{body weight}$)

#### Flip-Flop Kinetics & Continuity
When $k_a \to k_e$, the Bateman formula encounters numerical cancellation. The engine uses Taylor-expanded and continuous limit evaluations:

$$\lim_{k_a \to k_e} C(t) = \frac{F \cdot D \cdot k_a}{V_d} \cdot t \cdot e^{-k_e t}$$

And the time of maximum peak concentration ($T_{\max}$):

$$T_{\max} = \frac{\ln(k_a / k_e)}{k_a - k_e}, \quad \lim_{k_a \to k_e} T_{\max} = \frac{1}{k_e}$$

### 2. Intravenous (IV) Bolus Kinetics
Instantaneous systemic availability without an absorption phase ($F = 1$, $k_a \to \infty$):

$$C_{IV}(t) = \frac{D}{V_d} \, e^{-k_e t}, \quad T_{\max} = 0$$

### 3. Multiple Dosing & Steady-State Accumulation
Repeated administration at regular dosing interval $\tau$ leads to accumulation until drug intake equals clearance per interval:

* **Accumulation Index ($R$)**:
  $$R = \frac{1}{1 - e^{-k_e \tau}}$$

* **Steady-State Peak ($C_{ss,\max}$)** and **Trough ($C_{ss,\min}$)**:
  $$C_{ss,\min} = C(t = \tau) \quad \text{at steady state}$$

* **Average Steady-State Concentration ($C_{ss,\text{avg}}$)**:
  $$C_{ss,\text{avg}} = \frac{F \cdot D}{CL \cdot \tau} = \frac{\text{AUC}_{\tau}}{ \tau}$$
  Where total systemic clearance is $CL = k_e \cdot V_d$.

* **Time to Steady State**: Reached after $4\text{--}5 \times t_{1/2}$ ($\approx 93.75\%\text{--}96.88\%$ of asymptotic equilibrium).

### 4. Area Under the Curve ($\text{AUC}_{0\text{--}24}$)
Calculated using the linear trapezoidal integration method across discrete time intervals:

$$\text{AUC}_{t_1 - t_2} = \sum_{i=1}^{N-1} \frac{C_i + C_{i+1}}{2} (t_{i+1} - t_i)$$

---

## 🎯 Key Capabilities

| Feature | Description |
| :--- | :--- |
| **Exact PK Engine** | Bateman 1-compartment, IV bolus, and numerical ODE-validated profiles. |
| **Steady-State Analytics** | Computes $C_{ss,\max}$, $C_{ss,\min}$, $C_{ss,\text{avg}}$, accumulation factor $R$, and peak-trough fluctuation $\%$. |
| **Therapeutic Drug Monitoring (TDM)** | Real-time classification into **Toxic** ($> \text{MTC}$), **Therapeutic Window** ($\text{MEC}\text{--}\text{MTC}$), and **Sub-therapeutic** ($< \text{MEC}$). |
| **Missed Dose Simulation** | Evaluates clinical impact when a patient skips a dose or takes it delayed (+X hours). Calculates time to recover within $\pm 5\%$ of planned profile. |
| **Patient Variability Bands** | Shaded confidence envelope reflecting clinical $t_{1/2}$ inter-individual differences. |
| **Linear & Semi-Log Plotting** | Instant toggle between linear scale and semi-logarithmic ($\log_{10} C$ vs. $t$) view to inspect first-order elimination kinetics. |
| **SÚKL Registry Integration** | Live query into the Czech State Institute for Drug Control (SÚKL) open database with auto-dosing and direct registry links. |
| **Layperson Variable Guide** | Comprehensive glossary and interactive tooltips explaining all PK variables ($C_{\max}$, $C_{ss}$, $t_{1/2}$, $AUC$, etc.). |
| **Preset Drug Library** | Pre-configured with validated literature parameters for high-impact molecules. |
| **Audit & Protocol Export** | Download raw time-series to CSV and generate human-readable Markdown clinical dosing protocols. |
| **Bilingual Interface** | Seamless live switching between **English** and **Čeština**. |

---

## 💊 Clinical Drug Presets Library

The application includes curated pharmacokinetic profiles in [`data/drugs.json`](data/drugs.json):

| Drug | Clinical Category | $t_{1/2}$ (h) | $V_d$ (L/kg) | $F$ | $k_a$ (h⁻¹) | MEC | MTC | Unit | Key Clinical Context |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **Ibuprofen** | NSAID / Analgesic | $2.0$ | $0.15$ | $0.90$ | $2.00$ | $10.0$ | $50.0$ | $\mu\text{g/mL}$ | Rapid absorption; short half-life requires dosing every 6–8 h. |
| **Paracetamol** | Analgesic / Antipyretic | $2.5$ | $0.90$ | $0.80$ | $2.00$ | $10.0$ | $30.0$ | $\mu\text{g/mL}$ | High therapeutic margin; toxic metabolite threshold guarded by MTC. |
| **Amoxicillin** | $\beta$-lactam Antibiotic | $1.2$ | $0.30$ | $0.80$ | $1.50$ | $2.0$ | $40.0$ | $\mu\text{g/mL}$ | Time-dependent killing ($T > \text{MIC}$); rapid renal clearance. |
| **Lithium** | Mood Stabilizer | $24.0$ | $0.70$ | $0.95$ | $1.00$ | $0.6$ | $1.5$ | $\text{mmol/L}$ | **Narrow Therapeutic Index (NTI)**: severe toxicity $> 1.5\,\text{mmol/L}$. |
| **Diazepam** | Long-acting Benzodiazepine | $40.0$ | $1.10$ | $0.93$ | $2.00$ | $0.2$ | $0.8$ | $\mu\text{g/mL}$ | Prolonged elimination; pronounced steady-state accumulation. |
| **Livizux** | ADHD CNS Stimulant (lisdexamfetamin) | $11.0$ | $4.00$ | $0.95$ | $0.65$ | $0.02$ | $0.10$ | $\mu\text{g/mL}$ | Prodrug of d-amphetamine; once-daily morning dosing. AGNP 20–80 ng/mL. |
| **Levetiracetam** | Antiepileptic (Keppra) | $7.0$ | $0.60$ | $1.00$ | $1.80$ | $12.0$ | $46.0$ | $\mu\text{g/mL}$ | Complete oral absorption; AGNP reference range 12–46 µg/mL. |
| **Pregabalin** | Antiepileptic / Anxiolytic / Pain | $6.3$ | $0.56$ | $0.90$ | $1.80$ | $2.0$ | $8.0$ | $\mu\text{g/mL}$ | Linear renal clearance; low protein binding. AGNP target 2–5 µg/mL. |
| **Propranolol** | Non-selective $\beta$-blocker | $4.0$ | $4.00$ | $0.25$ | $1.50$ | $0.02$ | $0.20$ | $\mu\text{g/mL}$ | High first-pass metabolism ($F \approx 25\%$), large $V_d$ ($4\,\text{L/kg}$). |
| **Rivotril** | Benzodiazepine (clonazepam) | $35.0$ | $2.00$ | $0.90$ | $1.20$ | $0.02$ | $0.08$ | $\mu\text{g/mL}$ | High accumulation at steady-state. AGNP range 20–70 ng/mL. |
| **Mirtazapin** | Antidepressant (NaSSA) | $26.0$ | $1.50$ | $0.50$ | $1.00$ | $0.03$ | $0.15$ | $\mu\text{g/mL}$ | Once-daily dosing; AGNP therapeutic window 30–80 ng/mL. |
| **Risperidon** | Atypical Antipsychotic | $20.0$ | $1.20$ | $0.70$ | $1.50$ | $0.02$ | $0.08$ | $\mu\text{g/mL}$ | Modeled as total active moiety (risperidone + 9-OH); target 20–60 ng/mL. |
| **Custom Drug** | User-defined | *custom* | *custom* | *custom* | *custom* | *custom* | *custom* | *custom* | Fully configurable for novel molecules and in-vitro data. |

---

## 🏗️ Architecture & Codebase Structure

The project adheres to clean architecture principles, decoupling numerical PK algorithms from UI presentation:

```
druglevelmodelingapp/
├── app.py                      # Main Streamlit application and layout orchestration
├── core/                       # Pure numerical and pharmacological engines (zero UI dependencies)
│   ├── __init__.py
│   ├── pk_engine.py            # Bateman function, IV bolus, superposition, AUC, Steady-State
│   ├── drug_database.py        # Preset loader, validation, and PK model factory
│   ├── metrics.py              # Clinical metrics, TDM classification, crossing times & recovery
│   └── sukl_client.py          # SÚKL public REST API client, dosage parser, and ATC mapping
├── ui/                         # Presentation layer
│   ├── __init__.py
│   ├── charts.py               # Interactive Plotly chart builders with shaded TDM zones
│   ├── reporting.py            # Tabular data exports and clinical dosing protocol generation
│   └── i18n.py                 # Multi-language translations (EN, CS)
├── data/
│   └── drugs.json              # Validated clinical pharmacokinetic presets
├── tests/                      # Automated test suite (60 unit and smoke tests)
│   ├── __init__.py
│   ├── test_pk_engine.py       # Formula checks, ODE cross-validation, AUC and steady-state
│   ├── test_metrics.py         # Metric calculation, TDM classification and recovery
│   ├── test_drug_database.py   # Database consistency and unit conversion validation
│   ├── test_sukl_client.py     # SÚKL REST queries, dose extraction, and preset matching
│   ├── test_ui.py              # Plotly chart elements and i18n coverage
│   └── test_app_smoke.py       # Headless Streamlit smoke and UI interaction tests
├── docs/
│   └── wireframes.md           # Text wireframes and UI component layout specifications
├── .streamlit/
│   └── config.toml             # HealthTech theme styling and server configurations
├── .github/
│   └── workflows/ci.yml        # Automated GitHub Actions test pipeline
├── Dockerfile                  # Production containerization
├── requirements.txt            # Explicit dependency pinning
├── pytest.ini                  # Test discovery configuration
├── LICENSE                     # MIT License
└── README.md                   # Technical documentation
```

---

## 🚀 Quickstart

### Option A: Local Python Environment

1. **Clone the repository:**
   ```bash
   git clone https://github.com/sidd-1337/druglevelmodelingapp.git
   cd druglevelmodelingapp
   ```

2. **Create and activate a virtual environment:**
   ```bash
   # Linux / macOS
   python3 -m venv .venv
   source .venv/bin/activate

   # Windows (PowerShell)
   python -m venv .venv
   .\.venv\Scripts\Activate.ps1
   ```

3. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

4. **Launch the application:**
   ```bash
   streamlit run app.py
   ```
   Open `http://localhost:8501` in your browser.

---

### Option B: Docker Container

Run with a single command without installing Python locally:

```bash
# Build the Docker image
docker build -t druglevelmodelingapp .

# Run the container
docker run -d -p 8501:8501 --name druglevelmodelingapp druglevelmodelingapp
```

Then navigate to `http://localhost:8501`.

---

## 🧪 Testing & Verification

The mathematical formulas and UI components are covered by a suite of 54 automated pytest tests, including cross-validation against numerical ODE solvers (`scipy.integrate.solve_ivp`):

```bash
pytest -v
```

```
============================== test session starts ==============================
collected 54 items

tests/test_app_smoke.py::test_app_runs_without_exception PASSED            [  1%]
tests/test_app_smoke.py::test_app_all_drugs_both_languages PASSED          [  3%]
tests/test_app_smoke.py::test_app_missed_dose_and_semilog PASSED           [  5%]
tests/test_drug_database.py::test_all_presets_load_and_are_valid PASSED     [  7%]
tests/test_drug_database.py::test_specified_reference_values PASSED        [  9%]
tests/test_metrics.py::test_compute_metrics_single_dose PASSED             [ 11%]
tests/test_metrics.py::test_assess_toxic_reports_time PASSED               [ 12%]
tests/test_pk_engine.py::test_bateman_matches_ode_solution PASSED          [ 14%]
tests/test_pk_engine.py::test_auc_inf_equals_fd_over_cl PASSED             [ 16%]
tests/test_pk_engine.py::test_steady_state_matches_simulation PASSED       [ 18%]
...
============================== 54 passed in 4.5s ===============================
```

---

## 🛡️ Medical Disclaimer

> **IMPORTANT NOTICE:**  
> This software is intended strictly for educational, research, and exploratory simulation purposes. It is **not** a certified medical device and is **not** intended for primary clinical diagnosis, patient management, or definitive therapeutic dosing decisions. Real clinical dosing requires qualified physician and clinical pharmacist consultation, therapeutic drug monitoring assays, and comprehensive individual patient assessments (renal function, liver function, genetics, concomitant medications).

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).
