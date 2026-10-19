# Clinical Drug Level Modeling App
### *Pharmacokinetic Simulation (PK/PD) & Therapeutic Drug Monitoring (TDM) Engine*

[![Python](https://img.shields.io/badge/Python-3.12%20%7C%203.11-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Next.js](https://img.shields.io/badge/Next.js-15.5-000000?logo=next.js&logoColor=white)](https://nextjs.org/)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.115+-009688?logo=fastapi&logoColor=white)](https://fastapi.tiangolo.com/)
[![Streamlit](https://img.shields.io/badge/Streamlit-1.54+-FF4B4B?logo=streamlit&logoColor=white)](https://streamlit.io/)
[![Tests](https://img.shields.io/badge/Pytest-623%20passed-success?logo=pytest&logoColor=white)](tests/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An open-source clinical pharmacokinetic simulation platform built for healthtech engineers, pharmacologists, and clinical researchers. The system implements exact one-compartment pharmacokinetic models with oral absorption (Bateman function), intravenous (IV) bolus, intranasal (insufflated), and inhaled/smoked administration routes, steady-state accumulation analysis, missed-dose impact simulations, and interactive Therapeutic Drug Monitoring (TDM) with patient variability bounds.

It is paired with live Czech State Institute for Drug Control (SÚKL) open registries for pricing, reimbursement, drug shortage notifications (Section 33), and harm reduction guidance.

> **Languages**: English | [Čeština (README.cs.md)](README.cs.md)

---

## Table of Contents

- [Design Preview](#design-preview)
- [Quickstart & Getting Started](#quickstart--getting-started)
  - [Option A: Complete Stack via Docker Compose (Recommended)](#option-a-complete-stack-via-docker-compose-recommended)
  - [Option B: Local Development Environment](#option-b-local-development-environment)
- [Architecture & Codebase Structure](#architecture--codebase-structure)
- [Key Capabilities](#key-capabilities)
- [Scientific & Mathematical Foundation](#scientific--mathematical-foundation)
  - [1. Oral Absorption – Bateman Function](#1-oral-absorption--bateman-function)
  - [2. Intravenous (IV) Bolus Kinetics](#2-intravenous-iv-bolus-kinetics)
  - [3. Multiple Dosing & Steady-State Accumulation](#3-multiple-dosing--steady-state-accumulation)
  - [4. Area Under the Curve (AUC)](#4-area-under-the-curve-auc)
- [Clinical Drug Presets Library](#clinical-drug-presets-library)
- [Testing & Verification](#testing--verification)
- [Medical Disclaimer](#medical-disclaimer)
- [License](#license)

---

## Design Preview

The web portal ([`frontend/`](frontend/)) uses a quiet, Apple-style design: large type, hairline rows instead of cards, one blue accent and a single memorable element, the concentration curve. It is the logo mark, the line drawn on the home page and the interactive chart on every drug page. Light and dark themes follow the system setting, and the layout works down to phone width.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/preview/home-dark.png">
  <img src="docs/preview/home-light.png" alt="Home page with the headline "Co se s lékem děje v těle.", a drug search field and an illustrated dosing curve inside the therapeutic window">
</picture>

<table>
  <tr>
    <td width="50%" valign="top">
      <img src="docs/preview/drug-light.png" alt="Drug page for Ritalin with tabs and a list of registry data">
      <br><sub><b>Drug page.</b> Registry data as plain rows, four tabs: overview, price and copay, availability, body-level curve.</sub>
    </td>
    <td width="50%" valign="top">
      <picture>
        <source media="(prefers-color-scheme: dark)" srcset="docs/preview/kinetics-dark.png">
        <img src="docs/preview/kinetics-light.png" alt="Interactive concentration chart with the therapeutic window and dose, interval and weight sliders">
      </picture>
      <br><sub><b>Body-level curve.</b> Drag, touch or use the arrow keys to read the concentration at any time.</sub>
    </td>
  </tr>
</table>

<p align="center">
  <img src="docs/preview/home-mobile.png" width="240" alt="Home page on a phone">
  &nbsp;&nbsp;&nbsp;
  <img src="docs/preview/kinetics-mobile.png" width="240" alt="Concentration chart on a phone">
</p>

<sub>Screenshots show the full stack (Next.js and FastAPI) for the Ritalin 10 mg package. The portal itself is in Czech.</sub>

**Where the design lives**

- Colour, type and spacing are CSS variables in [`frontend/app/globals.css`](frontend/app/globals.css), mapped in [`frontend/tailwind.config.ts`](frontend/tailwind.config.ts).
- Shared UI primitives (Slider, Segmented, Tabs, Status, SpecList) are in [`frontend/components/ui/`](frontend/components/ui/).
- Typography: SF Pro on Apple devices, bundled Inter Variable elsewhere. Czech non-breaking spaces and decimal commas are handled in [`frontend/lib/typo.ts`](frontend/lib/typo.ts) and [`frontend/lib/format.ts`](frontend/lib/format.ts).

---

## Quickstart & Getting Started

### Option A: Complete Stack via Docker Compose (Recommended)

Run the complete multi-service production stack (Next.js web portal, FastAPI REST backend, and Streamlit deep laboratory) with a single command:

```bash
# 1. Clone the repository
git clone https://github.com/sidd-1337/druglevelmodelingapp.git
cd druglevelmodelingapp

# 2. Build and launch all services in background
docker compose up -d --build
```

#### Services overview:
* **Web Portal (Next.js 15, quiet Apple-style UI, see [Design Preview](#design-preview))**: [http://localhost:3000](http://localhost:3000)
  * Drug autocomplete search, SÚKL DLP data, SPC/PIL PDFs, pricing and copay protection calculators, shortage alerts, and harm reduction calculator.
* **REST API (FastAPI)**: [http://localhost:8000](http://localhost:8000)
  * Interactive OpenAPI Swagger documentation available at [http://localhost:8000/docs](http://localhost:8000/docs).
* **Clinical Laboratory (Streamlit)**: [http://localhost:8501](http://localhost:8501)
  * Deep pharmacometric modeling: missed dose simulations, semi-logarithmic axes, confidence envelopes, and raw CSV exports. Light and night themes follow the system setting; switch manually in the app menu (⋮ → Settings).

To stop the containers:
```bash
docker compose down
```

---

### Option B: Local Development Environment

#### 1. Python Backend & Core Setup

```bash
# Create and activate virtual environment
# Linux / macOS:
python3 -m venv .venv
source .venv/bin/activate

# Windows (PowerShell):
python -m venv .venv
.\.venv\Scripts\Activate.ps1

# Install Python dependencies
pip install -r requirements.txt

# Run FastAPI backend (port 8000)
uvicorn backend.app.main:app --host 0.0.0.0 --port 8000 --reload

# In a separate terminal: Run Streamlit simulator (port 8501)
streamlit run app.py
```

#### 2. Next.js Frontend Setup

```bash
cd frontend
npm install
npm run dev
```

Open [http://localhost:3000](http://localhost:3000) in your browser.

---

## Architecture & Codebase Structure

The project follows a decoupled, clean architecture separating pure mathematical modeling from API orchestration and modern presentation layers:

```
druglevelmodelingapp/
├── backend/                    # FastAPI production REST service
│   ├── app/
│   │   ├── main.py             # FastAPI entrypoint, CORS configuration
│   │   ├── routers/            # Endpoint routers: /drugs, /pk
│   │   └── schemas/            # Pydantic v2 validation models
│   └── tests/                  # Integration tests for REST endpoints
├── core/                       # Pure numerical and pharmacological core (zero UI dependencies)
│   ├── pk_engine.py            # Bateman function, IV bolus, superposition, AUC, Steady-State
│   ├── harm_reduction.py       # Safe redosing thresholds, dose stacking analysis, phase timings
│   ├── drug_database.py        # Preset loader, validation, and PK model factory
│   ├── metrics.py              # Clinical metrics, TDM classification, crossing times & recovery
│   └── sukl_client.py          # SÚKL public REST API client, dosage parser, and ATC mapping
├── frontend/                   # Next.js 15 web portal (quiet Apple-style design, see Design Preview)
│   ├── app/                    # App Router: /, /leky/[kod], /slovnik
│   ├── components/             # ui/ primitives, HeroCurve, PkChart (SVG), CopayCalculator, ThemeToggle
│   └── lib/                    # Client-side API client, PK calculation engine, SÚKL cache, Czech formatting
├── ui/                         # Streamlit presentation components
│   ├── charts.py               # Interactive Plotly chart builders with shaded TDM zones
│   ├── reporting.py            # Tabular data exports and clinical dosing protocol generation
│   └── i18n.py                 # Multi-language translations (EN, CS)
├── data/
│   └── drugs.json              # Validated clinical pharmacokinetic presets
├── tests/                      # Automated test suite (109 unit and integration tests)
├── docker-compose.yml          # Multi-container orchestration (backend, frontend, streamlit)
├── Dockerfile                  # Streamlit container definition
├── Dockerfile.backend          # FastAPI container definition
└── README.md                   # Technical documentation
```

---

## Key Capabilities

| Feature | Description |
| :--- | :--- |
| **Exact PK Engine** | Bateman 1-compartment, IV bolus, and numerical ODE-validated profiles. |
| **Multiple Routes of Administration** | Supports Oral ($p.o.$), Intravenous ($i.v.$ bolus), Insufflated (intranasal), and Smoked/Inhaled. |
| **Steady-State Analytics** | Computes $C_{ss,\max}$, $C_{ss,\min}$, $C_{ss,\text{avg}}$, accumulation factor $R$, and peak-trough fluctuation $\%$. |
| **Therapeutic Drug Monitoring (TDM)** | Real-time classification into **Toxic** ($> \text{MTC}$), **Therapeutic Window** ($\text{MEC}\text{--}\text{MTC}$), and **Sub-therapeutic** ($< \text{MEC}$). |
| **Co-medication (multi-drug)** | Several drugs on one chart, each as its position in its own therapeutic window; a shared enzyme a drug destroys (MDMA→CYP2D6) is simulated as one pool. |
| **Interaction time course** | Each interaction keeps its measured or derived factor but builds and fades with the perpetrator's pharmacokinetics; inhibitors over their half-life, inducers over ~a week; start/stop days supported. |
| **Pharmacodynamic clash check** | Flags serotonergic, CNS-depressant, QT-prolonging or sympathomimetic effects shared by two drugs that are active at once, from primary labels. |
| **Harm Reduction Calculator** | Safe redosing countdown timer, dose stacking prevention alerts, and phase timeline estimators. |
| **Missed Dose Simulation** | Evaluates clinical impact when a patient skips a dose or takes it delayed (+X hours). Calculates recovery time. |
| **Patient Variability Bands** | Shaded confidence envelope reflecting clinical $t_{1/2}$ inter-individual differences. |
| **Linear & Semi-Log Plotting** | Instant toggle between linear scale and semi-logarithmic ($\log_{10} C$ vs. $t$) view to inspect first-order elimination kinetics. |
| **SÚKL Registry Integration** | Live query into the Czech State Institute for Drug Control (SÚKL) open database with auto-dosing, shortage status, and reimbursement tiers. |
| **Copay Limit Calculator** | Annual out-of-pocket statutory protection limit calculator (5 000 / 1 000 / 500 CZK). |
| **Design Aesthetic** | Quiet Apple-style interface with Light & Night mode (system preference auto-detection, zero-FOUC initialization), SF Pro / Inter typography, keyboard-operable controls and a responsive layout down to phone width. |

---

## Scientific & Mathematical Foundation

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

### 4. Area Under the Curve (AUC)
Calculated using the linear trapezoidal integration method across discrete time intervals:

$$\text{AUC}_{t_1 - t_2} = \sum_{i=1}^{N-1} \frac{C_i + C_{i+1}}{2} (t_{i+1} - t_i)$$

---

## Clinical Drug Presets Library

The application includes curated pharmacokinetic profiles in [`data/drugs.json`](data/drugs.json). Values describe a 70 kg adult with normal renal function; $V_d$ and $k_a$ are fitted so that AUC / steady-state level, Cmax and Tmax match published data in [`data/reference_pk.json`](data/reference_pk.json). Each preset also carries the between-patient CV of clearance (and bioavailability), the renally cleared fraction, two-compartment distribution and absorption lag where published data support them, selectable formulations (propranolol immediate vs. prolonged release), and for MDMA mechanism-based CYP2D6 inactivation.

| Drug | Clinical Category | $t_{1/2}$ (h) | $V_d$ (L/kg) | $F$ | $k_a$ (h⁻¹) | MEC | MTC | Unit | Key Clinical Context |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **Ibuprofen** | NSAID / Analgesic | $2.0$ | $0.13$ | $0.90$ | $2.40$ | $10.0$ | $50.0$ | $\mu\text{g/mL}$ | Rapid absorption; short half-life requires dosing every 6–8 h. |
| **Paracetamol** | Analgesic / Antipyretic | $2.5$ | $0.70$ | $0.80$ | $3.80$ | $10.0$ | $30.0$ | $\mu\text{g/mL}$ | High therapeutic margin; toxic metabolite threshold guarded by MTC. |
| **Amoxicillin** | $\beta$-lactam Antibiotic | $1.0$ | $0.28$ | $0.80$ | $0.55$ | $2.0$ | $40.0$ | $\mu\text{g/mL}$ | Time-dependent killing ($T > \text{MIC}$); rapid renal clearance. |
| **Lithium** | Mood Stabilizer | $24.0$ | $0.70$ | $0.95$ | $1.00$ | $0.6$ | $1.5$ | $\text{mmol/L}$ | **Narrow Therapeutic Index (NTI)**: severe toxicity $> 1.5\,\text{mmol/L}$. |
| **Diazepam** | Long-acting Benzodiazepine | $40.0$ | $1.20$ | $0.93$ | $4.00$ | $0.2$ | $0.8$ | $\mu\text{g/mL}$ | Prolonged elimination; pronounced steady-state accumulation. |
| **Livizux** | ADHD CNS Stimulant (lisdexamfetamin) | $11.0$ | $3.35$ | $0.95$ | $0.67$ | $0.02$ | $0.10$ | $\mu\text{g/mL}$ | Prodrug of d-amphetamine; once-daily morning dosing. AGNP 20–80 ng/mL. |
| **Levetiracetam** | Antiepileptic (Keppra) | $7.0$ | $0.60$ | $1.00$ | $4.50$ | $12.0$ | $46.0$ | $\mu\text{g/mL}$ | Complete oral absorption; AGNP reference range 12–46 µg/mL. |
| **Pregabalin** | Antiepileptic / Anxiolytic / Pain | $6.3$ | $0.59$ | $0.90$ | $3.90$ | $2.0$ | $8.0$ | $\mu\text{g/mL}$ | Linear renal clearance; low protein binding. AGNP target 2–5 µg/mL. |
| **Propranolol** | Non-selective $\beta$-blocker | $4.0$ | $2.70$ | $0.25$ | $1.60$ | $0.02$ | $0.20$ | $\mu\text{g/mL}$ | High first-pass metabolism ($F \approx 25\%$, CV of exposure about 50 %); apparent $V_d$ fitted to published Cmax/AUC. |
| **Rivotril** | Benzodiazepine (clonazepam) | $35.0$ | $2.30$ | $0.90$ | $1.50$ | $0.02$ | $0.08$ | $\mu\text{g/mL}$ | High accumulation at steady-state. AGNP range 20–70 ng/mL. |
| **Mirtazapin** | Antidepressant (NaSSA) | $26.0$ | $10.50$ | $0.50$ | $1.00$ | $0.03$ | $0.15$ | $\mu\text{g/mL}$ | Once-daily dosing; AGNP therapeutic window 30–80 ng/mL. |
| **Risperidon** | Atypical Antipsychotic | $20.0$ | $1.70$ | $0.70$ | $2.00$ | $0.02$ | $0.08$ | $\mu\text{g/mL}$ | Modeled as total active moiety (risperidone + 9-OH); target 20–60 ng/mL. |
| **Ritalin** | CNS Stimulant / ADHD (methylphenidate) | $2.8$ | $3.10$ | $0.15$ | $1.30$ | $0.005$ | $0.030$ | $\mu\text{g/mL}$ | Rapid absorption; short $t_{1/2}$ requires dosing 2–3× daily ($\tau = 8\text{ h}$). Target 5–20 ng/mL. |
| **Custom Drug** | User-defined | *custom* | *custom* | *custom* | *custom* | *custom* | *custom* | *custom* | Fully configurable for novel molecules and in-vitro data. |

---

## Testing & Verification

The test suite covers numerical formulas, ODE cross-validation, SÚKL registry data extraction, and FastAPI REST endpoints across 623 passing automated tests (the SÚKL live-lookup tests need network access). `tests/test_reference_pk.py` checks every preset against published Cmax, Tmax, AUC and steady-state levels; `python -m core.validation` prints the predicted/observed table (exposure 1.08x, Cmax 1.14x, Tmax 1.10x, MDMA redosing 1.06x geometric mean fold error). The reference values were collected from search excerpts of the cited labels and papers (`verified: false`) and still need checking against the full documents.

```bash
# Run test suite
pytest -q
```

---

## Medical Disclaimer

> **IMPORTANT NOTICE:**  
> This software is intended strictly for educational, research, and exploratory simulation purposes. It is **not** a certified medical device and is **not** intended for primary clinical diagnosis, patient management, or definitive therapeutic dosing decisions. Real clinical dosing requires qualified physician and clinical pharmacist consultation, therapeutic drug monitoring assays, and comprehensive individual patient assessments (renal function, liver function, genetics, concomitant medications).

---

## License

This project is licensed under the [MIT License](LICENSE).
