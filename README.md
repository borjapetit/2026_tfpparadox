<br>**Replication codes**

**Firm Entry and the TFP Paradox of Productivity-Dependent Distortions**<br>
Borja Petit, Miguel Almunia, Juan F. Jimeno, and David López-Rodríguez.<br>

Compile with `make model` and run `./model` from this folder. The menu asks for:

1. the economy: Hopenhayn (1992) or Lucas (1978);
2. the exercise: solve, calibrate, symmetric experiment or its sensitivity analysis, for the baseline model or its alternative versions (no capital; fixed and entry costs in labor units, Hopenhayn only);
3. the profit tax rate: 30% or 0%;
4. in the Hopenhayn economy, the exogenous exit probability: 0%, 2%, or no endogenous exit (exit probability equal to the exit rate in the data, 8%, and no fixed cost).

File names (`<economy>` is `hopenhayn_tau<tau>_delta<delta>` or `lucas_tau<tau>`, with tau and delta in %):

| | |
|---|---|
| parameters (input) | `txtfiles/calib_<economy>_<version>.txt` |
| calibration | `results/<date>_calib_<economy>_<version>.txt` |
| solution | `results/<date>_solve_<economy>_<version>.txt` |
| experiments | `results/<economy>_<version>.txt`, `results/<economy>_<exercise>.txt` |

`<version>` is `baseline`, `nocapital` or `costs_labor`. The sensitivity exercises are `gamma_up5`, `gamma_down5`, `sigma0_up10`, `sigma0_down10`, `frisch_0-01`, `frisch_3` (Hopenhayn) and `gamma_up5`, `gamma_down5`, `xi_up10`, `xi_down10` (Lucas). A calibration is not used automatically: to use it, copy it to `txtfiles/` and remove the date from its name.

`run_paper.sh` runs the symmetric experiment and the sensitivity analysis for all the economies of the paper in parallel.

If you find any error or have a suggestion, do not hesitate to [contact](mailto:bpetit@cunef.edu) us.

_This code is distributed under the [MIT Licence](LICENSE)_
