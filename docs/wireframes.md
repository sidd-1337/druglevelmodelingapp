# UI Architecture & Wireframes

This document details the interface layout, information architecture, visual hierarchy, and component interaction flows for the **Drug Level Modeling App (PK/PD Simulation & Therapeutic Drug Monitoring)**.

---

## 1. High-Level Layout Architecture

The application is built on a responsive, dual-column structure:
- **Left Sidebar (`280px–340px`)**: Parameter configuration, drug selection, dosing regimen controls, display toggles, and missed-dose scenario controls.
- **Main Canvas**: KPI metric ribbon, contextual clinical alert banner, interactive Plotly concentration-time visualization, and detailed analytical tabs.

```
+----------------------------------------------------------------------------------------------------+
| 💊 Drug Level Modeling – PK/PD Simulation & TDM                                                    |
| One-compartment pharmacokinetics with oral absorption (Bateman), multiple dosing & steady state    |
+-------------------+--------------------------------------------------------------------------------+
| SIDEBAR           | MAIN CANVAS                                                                    |
|                   |                                                                                |
| [Lang: EN | CS]   | [ C_max ]  [ T_max ]  [ AUC_0-24 ]  [ Clearance ]  [ t1/2 ]  [ Time to SS ]  |
|                   | +----------------------------------------------------------------------------+ |
| 1. Drug & Patient | | ⚠️ CLINICAL ALERT BANNER (e.g. Warning: C exceeds MTC at T=6.5 h)          | |
| - Select Drug     | +----------------------------------------------------------------------------+ |
| - Body Weight     |                                                                                |
| - Editable PK:    | +----------------------------------------------------------------------------+ |
|   t1/2, Vd, F, ka | | INTERACTIVE PLOTLY CHART (Pan / Zoom / Hover Tooltips)                     | |
|   MEC, MTC        | |  Concentration (ug/mL or mmol/L) vs Time (hours)                           | |
|                   | |  [/// Red Hatched Zone: Toxic (> MTC)                                  ///]| |
| 2. Regimen        | |  [    Green Shaded Zone: Therapeutic Window (MEC - MTC)                 ]  | |
| - Route (Oral/IV) | |  [... Grey Shaded Zone: Sub-therapeutic (< MEC)                       ...]  | |
| - Dose (mg)       | |  --- Solid Blue Curve: Concentration C(t)                                  | |
| - Interval tau (h)| |  --- Dashed Grey Curve: Planned Ref (when simulating missed dose)          | |
| - Duration (days) | |  ░░░ Light Blue Band: Patient t1/2 variability (min - max range)           | |
|                   | |  ▲ Dose administration markers  |  ✕ Missed dose marker                    | |
| 3. Display        | +----------------------------------------------------------------------------+ |
| - Linear / SemiLog|                                                                                |
| - Toggles (Zones) | [ Tab 1: Steady State ] [ Tab 2: Missed Dose ] [ Tab 3: Data & Export ] [ Tab 4: Model ]
|                   | +----------------------------------------------------------------------------+ |
| 4. Missed Dose    | | TAB CONTENT AREA                                                           | |
| - Enable checkbox | | - Steady-state analytical metrics & single-interval profile                | |
| - Dose index      | | - Missed dose trough impact & recovery time calculation                    | |
| - Skip vs Delay   | | - Dataframe preview, CSV export, Clinical Protocol Markdown export        | |
+-------------------+ +----------------------------------------------------------------------------+ |
|                   | ⚕️ Medical Disclaimer: Educational and modeling tool only.                      |
+-------------------+--------------------------------------------------------------------------------+
```

---

## 2. Component Specifications

### 2.1 Sidebar Controls

1. **Language Selector (`st.radio`)**
   - Options: `English`, `Čeština`.
   - Real-time reactivity switching all labels, tooltips, warnings, and generated markdown protocols.

2. **Drug Selection & Presets (`st.selectbox`)**
   - Populated dynamically from `data/drugs.json` + `Custom drug`.
   - Populates validated default values:
     - **Ibuprofen**: $t_{1/2} = 2.0\,\text{h}$, $V_d = 0.15\,\text{L/kg}$, $F = 0.90$, $\text{MEC} = 10$, $\text{MTC} = 50\,\mu\text{g/mL}$.
     - **Paracetamol**: $t_{1/2} = 2.5\,\text{h}$, $V_d = 0.90\,\text{L/kg}$, $F = 0.80$, $\text{MEC} = 10$, $\text{MTC} = 30\,\mu\text{g/mL}$.
     - **Amoxicillin**: $t_{1/2} = 1.2\,\text{h}$, $V_d = 0.30\,\text{L/kg}$, $F = 0.80$, $\text{MEC} = 2.0$, $\text{MTC} = 40\,\mu\text{g/mL}$.
     - **Lithium**: $t_{1/2} = 24.0\,\text{h}$, $V_d = 0.70\,\text{L/kg}$, $F = 0.95$, $\text{MEC} = 0.6$, $\text{MTC} = 1.5\,\text{mmol/L}$.
     - **Diazepam**: $t_{1/2} = 40.0\,\text{h}$, $V_d = 1.10\,\text{L/kg}$, $F = 0.93$, $\text{MEC} = 0.2$, $\text{MTC} = 0.8\,\mu\text{g/mL}$.
     - **Custom**: Fully editable parameters.

3. **Editable Pharmacokinetic Drawer (`st.expander`)**
   - Numeric inputs for $t_{1/2}$, $V_d$, $F$, $k_a$, $\text{MEC}$, $\text{MTC}$, and population range $[t_{1/2,\min}, t_{1/2,\max}]$.

4. **Dosing Regimen**
   - Administration Route: `Oral` vs. `IV bolus`.
   - Dose ($D$ in mg), Interval ($\tau \in \{6, 8, 12, 24\}\,\text{h}$ or custom), and Duration ($N$ days).

5. **Missed-Dose Module**
   - Checkbox to toggle scenario.
   - Selection of scheduled dose index $i$.
   - Mode: `Skips dose` vs. `Takes late (+X hours)`.

---

### 2.2 KPI Dashboard Ribbon

Placed above the main chart for instant clinical triage:

```
+---------------+---------------+---------------+---------------+---------------+---------------+
| C_max         | T_max         | AUC_0-24      | Clearance     | t1/2          | Time to SS    |
| 28.4 ug/mL    | 1.62 h        | 342 ug*h/mL   | 3.64 L/h      | 2.0 h         | 9 h (4.5*t1/2)|
+---------------+---------------+---------------+---------------+---------------+---------------+
```

---

### 2.3 Clinical Alert Banner (`st.error` / `st.warning` / `st.success`)

Evaluated dynamically against therapeutic boundaries:
- **Toxic Exposure (`status = "toxic"`)**:
  `⚠️ Pozor: Odhadovaná hladina překračuje bezpečný limit MTC (50 ug/mL) v čase T=6.5 h! Čas nad MTC: 2.1 h, maximum 58.2 ug/mL.`
- **Sub-therapeutic Troughs (`status = "trough_below_mec"`)**:
  `⚠️ Údolní hladiny (6.4 ug/mL) klesají pod MEC (10 ug/mL) – riziko subterapeutické expozice.`
- **Ineffective Regimen (`status = "never_effective"`)**:
  `⚠️ Koncentrace nikdy nedosáhne MEC (10 ug/mL); maximum je jen 4.2 ug/mL – schéma je neúčinné.`
- **Therapeutic Equilibrium (`status = "within_window"`)**:
  `✅ Dávkování udržuje stabilní koncentraci v terapeutickém okně.`

---

### 2.4 Interactive Plotly Visualization

Features:
- **Zones**:
  - `> MTC`: Diagonal red hatching (`fillpattern=dict(shape="/")`).
  - `MEC – MTC`: Semi-transparent green band (`rgba(22, 163, 74, 0.16)`).
  - `< MEC`: Sub-therapeutic soft slate grey (`rgba(100, 116, 139, 0.22)`).
- **Traces**:
  - Main simulation curve $C(t)$ (Solid Royal Blue, width 3).
  - Planned regimen reference curve (Dashed Grey, width 2) shown during missed-dose simulations.
  - Inter-patient variability envelope (Shaded band based on $[t_{1/2,\min}, t_{1/2,\max}]$).
  - Upward triangle markers on x-axis representing dose administrations ($\blacktriangle$).
  - Red cross marker ($\mathbf{\times}$) indicating missed or delayed dose event.
- **Axes**:
  - Instant toggle between Linear scale and Semi-logarithmic scale ($\log_{10} C$ vs. $t$).
  - Unified hover tooltip with high-precision time and concentration readings.

---

### 2.5 Tabbed Deep-Dive Views

1. **Tab 1: Steady State Analysis**
   - Metrics: $C_{ss,\max}$, $C_{ss,\min}$ (trough), $C_{ss,\text{avg}}$, Accumulation index $R = \frac{1}{1 - e^{-k_e \tau}}$, Peak-trough fluctuation $\%$, and analytical $\text{AUC}_{\tau,ss}$.
   - Single steady-state cycle plot ($[0, \tau]$) showing exact asymptotic cyclic behavior.

2. **Tab 2: Missed Dose Analysis**
   - Focused timeline window around the missed dose event ($T_{\text{missed}} - 1.5\tau$ to $T_{\text{missed}} + 4\tau$).
   - Delta metrics: Minimum drop below planned trough, hours spent below MEC, hours spent above MTC (if compensatory double dose was taken), and time to recover within $\pm 5\%$ of planned profile.

3. **Tab 3: Data & Export**
   - Single-click CSV download containing time grid, concentrations, variability bounds, and zone tags.
   - Comprehensive Dosing Protocol export in GitHub-flavored Markdown.
   - Formatted interactive data table.

4. **Tab 4: Mathematical Model Reference**
   - Formal mathematical definitions: Bateman oral absorption, IV bolus monoexponential decline, clearance formulas, accumulation index, and flip-flop kinetics.
