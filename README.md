# LimoSim

LimoSim is an online COBRA simulator built by the CCU-Taiwan iGEM team on the same genome-scale *E. coli* model used in MATLAB, `iEC1356_Bl21DE3`, letting anyone estimate limonene production without needing a COBRA background — and letting COBRA users run full FBA/dFBA analyses through a web front end instead of hand-typed MATLAB scripts.

**Live demo (Simple mode):** https://ccu-igem-software.onrender.com

## Overview

- **Simple mode:** Runs entirely in the browser (phone or computer) via a pre-fitted quick-estimate equation. Adjust each enzyme's `Etot` (DXS, IDI, GPPS, LS) to instantly get an estimated limonene concentration, termination classification, and estimated end time — no installation, no simulation background required.
- **Advanced mode:** Opens the full COBRA console, with all adjustable bounds (medium uptake limits, ATPM, pathway bounds, induction timing, etc.) and four analysis tabs. Runs live COBRApy solves, so it needs more memory/CPU than Simple mode and is meant to be run locally.
  - **FBA:** Two-stage solve — maximize growth rate first, then fix growth at a lower bound and maximize limonene secretion flux within it. Results include max growth rate, production-phase growth rate, limonene flux, per-enzyme fluxes, and ATP/NADPH rates.
  - **dFBA:** Repeatedly re-solves FBA over time (fixed or adaptive time steps) to simulate biomass, glycerol, limonene, and growth rate before/after induction, producing the characteristic flat-then-rising limonene curve. Runs are saved to a history log and can be overlaid in Compare mode.
  - **Etot Analysis:** Three sub-modes — manual scan (re-solves dFBA across a value series for one enzyme, log-scale chart), auto-optimize (Nelder–Mead search for the limonene-maximizing enzyme value), and quick estimate (same surrogate equation as Simple mode).
  - **Multi-Enzyme Scan:** Joint response surface across 2–4 enzymes, via exhaustive grid or Latin-Hypercube/Sobol sampling, shown as heatmaps.
  - All results can expand into a raw data table and be exported as CSV or an image; interface switches instantly between Chinese/English with a light/dark theme toggle.

### Methods

| Function | Computes | Method |
|---|---|---|
| FBA | Growth / limonene flux at steady state | `maximize cᵀv`, s.t. `Sv=0`, `lb ≤ v ≤ ub` (Orth, Thiele & Palsson, 2010) |
| Enzyme-capacity constraint | Caps a reaction's flux by user `Etot` | `v ≤ kcat · Etot`; kcat (1/h): DXS 4.1, GPPS 1.4, IDI 0.33, LS 0.0905 |
| dFBA | Time-course accumulation after induction | Static optimization approach (Mahadevan, Edwards & Doyle, 2002) |
| Surrogate quick estimate | Instant limonene estimate | Curve pre-fitted offline against real dFBA runs |
| Auto-optimize | Limonene-maximizing enzyme value | Nelder–Mead simplex (Nelder & Mead, 1965) |
| Multi-enzyme scan | Joint response surface | Latin-Hypercube (McKay, Beckman & Conover, 1979) / Sobol (Sobol, 1967) sampling |

Manual scan, auto-optimize, and multi-enzyme scan all run through a shared worker pool (`ProcessPoolExecutor`).

### Other pages

- **Related Information page:** Why COBRA is used, where the model/parameters come from, what FBA and dFBA each answer, methods not yet covered (gene knockouts, OptKnock, flux sampling), validation findings, how to read result limits, and a changelog from Limonene001 to the current build.
- **About Us page:** Team introduction with links to GitHub, the live tool, Instagram, the school site, YouTube, and the team's other projects — [NoFold](https://leilong20060926.github.io/CCU-iGEM-game/), the [pest outbreak advisory system](https://leilong20060926.github.io/CCU-iGEM-model/), and the [Taiwan Organic Agriculture Map](https://igem.xn--hrr.tw/).

---

## Development

Technical notes for anyone working on the codebase.

### Getting Started

#### Use Simple Mode Online
No installation needed — just open:
https://ccu-igem-software.onrender.com

#### Run Advanced Mode Locally

Advanced mode runs full FBA/dFBA/scan solves through live COBRApy, which needs more memory and CPU than Render's free tier reliably provides — so for the full analysis console, run it on your own machine.

##### Step 1: Install the required tools
You need **Git** and **Python 3.10+**. Check whether you already have them by opening a terminal
(Windows: *Command Prompt* or *PowerShell*; macOS: *Terminal*; Linux: your terminal app) and running:

```bash
git --version
python3 --version
```

> On Windows, use `python --version` instead of `python3 --version`.

If either command shows an error, install it first:
- Git: https://git-scm.com/downloads
- Python 3: https://www.python.org/downloads/
  (Windows users: check **"Add Python to PATH"** during installation.)

##### Step 2: Download the project

**Option A — Using Git (recommended)**
```bash
git clone https://github.com/Leilong20060926/CCU-iGEM-software.git
cd limosim
```

**Option B — Without Git**
1. Go to https://github.com/Leilong20060926/CCU-iGEM-software.git
2. Click the **Code** button → **Download source code** → **zip**
3. Unzip the file, then open a terminal inside the unzipped folder (the one containing `main.py`).

##### Step 3: Install Python dependencies
```bash
python3 -m venv venv
source venv/bin/activate
```
> On Windows, use `venv\Scripts\activate`.
```bash
pip install -r requirements.txt
```

##### Step 4: Start the server
```bash
python3 main.py
```
You should see a message like `Running on http://127.0.0.1:5001`. Keep this terminal window open while using the tool.

##### Step 5: Open the console
Open your browser and go to:
```
http://localhost:5001
```

##### Step 6: Stop the server
When you're done, go back to the terminal and press **Ctrl + C**.

#### Troubleshooting
- **"Address already in use"** — port 5001 is taken. Use another port, e.g. `PORT=5002 python3 main.py`, then open `http://localhost:5002`.
- **Model fails to load / "model not found"** — make sure you ran the command *inside* the project folder (the one containing `main.py`).
- **Solves are slow on my machine** — LimoSim defaults to the HiGHS solver and parallelizes scans across CPU cores. On a low-memory machine, cap worker processes with `LIMONENE_MAX_POOL_WORKERS=2 python3 main.py`.
- **Buttons don't respond / blank page** — hard-refresh the page (Ctrl/Cmd+Shift+R) to clear a stale cached script.

### References

- COBRApy — https://opencobra.github.io/cobrapy/
- COBRA Toolbox v3.0 (MATLAB) — https://opencobra.github.io/cobratoolbox/stable/
- Orth, Thiele & Palsson, 2010, *Nature Biotechnology* — What is flux balance analysis?
- Mahadevan, Edwards & Doyle, 2002, *Biophysical Journal* — Dynamic Flux Balance Analysis of Diauxic Growth in *Escherichia coli*
- Nelder & Mead, 1965, *Computer Journal* — A Simplex Method for Function Minimization
- McKay, Beckman & Conover, 1979, *Technometrics* — A Comparison of Three Methods for Selecting Values of Input Variables...
- Sobol, 1967, *USSR Computational Mathematics and Mathematical Physics* — On the distribution of points in a cube and the approximate evaluation of integrals
- HiGHS solver — https://highs.dev/
- Plotly.js — https://plotly.com/javascript/
