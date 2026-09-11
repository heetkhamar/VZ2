# VZ2 – Verzinnungsanlage Thermal Simulation

Physics-based thermal simulation of a tin-coating line (Verzinnungsanlage).
The model solves the energy balance of the tin bath (heating, radiation/convection
losses, strip heat extraction, air-knife and water cooling) via forward-Euler
integration.

- `simulation.py` – simulation model, run loop, and plotting
- `parameters.py` – all physical/operational parameters
- `*.txt` – furnace measurement data (CUSN6, CUSN8, KHP…)
- `simulation_results.png` – output plot (regenerated on each run)

## Run it yourself on the web (GitHub Codespaces)

This repo ships a devcontainer so you can run the simulation in the browser
with no local setup — and **without spending Claude tokens** on the runs
themselves.

1. On the GitHub repo page, click **Code ▸ Codespaces ▸ Create codespace**.
2. Wait for the container to build. `numpy` and `matplotlib` are installed
   automatically from `requirements.txt`.
3. In the terminal, run:
   ```bash
   python simulation.py
   ```
4. Open `simulation_results.png` in the file explorer — VS Code renders the
   plot inline.

### Claude as an assistant (only when you want changes)

The devcontainer preloads the **Claude Code** VS Code extension. Sign in when
you want Claude to edit the code (e.g. change parameters, add an analysis,
plot the furnace data). Claude tokens are used **only** while you use the
extension — running `python simulation.py` yourself never costs tokens.

## Configure a run

Edit the `SIMULATION INPUTS` block at the bottom of `simulation.py` (duration,
heating schedule, strip parameters, air-knife pressures, water cooling, initial
bath temperature), or import `run_simulation()` from your own script. Physical
constants and lookup tables live in `parameters.py`.

## Run locally instead

```bash
pip install -r requirements.txt
python simulation.py
```
