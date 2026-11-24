# Validated results — 2026-11-24

All values below were generated from the repository with MATLAB R2025b. The
method-level validation plan passes its structural audit, and the regression
suite rejects a planted wrong sign on the reaction term.

## Main findings

1. **The implementation reproduces the paper.** At the four paper resolutions,
   the call prices are `3.1676`, `3.7501`, `3.8871`, and `3.9189`; the cash-call
   prices are `50.489`, `49.577`, `49.300`, and `49.241` after rounding.
2. **The paper's printed exact call benchmark is wrong.** The standard formula
   gives `3.929276`, not `3.9123`. With the corrected oracle, the call errors
   decrease with observed orders `2.09`, `2.09`, and `2.02`.
3. **The centered spatial operator is second-order consistent.** Its max-norm
   errors over $S\in[60,140]$ are `7.97e-3`, `1.99e-3`, `4.99e-4`, and
   `1.25e-4` as $h$ halves, giving orders `2.000`, `2.000`, and `2.000`.
4. **A fixed-$N_t$ full-solver h scan is not an isolated spatial study.** Its
   errors are non-monotone because sweep/time and spatial errors couple and
   cancel. That negative result is retained in `space_coupled_solver_scan.csv`.
5. **The put extension passes two independent checks.** The put target errors
   have the same refinement sequence as the call, and put--call parity residuals
   remain between about `2e-14` and `6e-12`.
6. **ADE survives the declared stress ladder while Forward Euler fails.** For
   $h=0.5$, $S_{max}=400$, and $N_t=20,40,80,160,320,500$, every ADE run remains
   finite and nonnegative. Forward Euler is classified as `BLOWUP` or
   `NONFINITE` in every case. This is empirical evidence on one parameter
   ladder, not a general stability proof.

## Paper reproduction with corrected call oracle

| $N_t$ | ADE call | Correct exact | Absolute error | Observed order |
|---:|---:|---:|---:|---:|
| 160 | 3.1676 | 3.929276 | 0.76164 | — |
| 320 | 3.7501 | 3.929276 | 0.17916 | 2.0879 |
| 640 | 3.8871 | 3.929276 | 0.042143 | 2.0879 |
| 1280 | 3.9189 | 3.929276 | 0.010415 | 2.0166 |

## Claim status

| Claim | Status | Scope |
|---|---|---|
| Paper tables are reproduced and corrected analytical formulas are matched under refinement | **VERIFIED** | Stated call and cash-call parameters and four resolutions |
| Approximately second-order time behavior and second-order spatial-operator consistency | **CONDITIONAL** | Declared ladders; not a proof of full-solver spatial order or unconditional stability |
| European put extension is correct | **VERIFIED** | Stated parameters, four resolutions, analytical price, boundary value, and put--call parity |
| ADE is unconditionally stable | **UNTESTED** | The stress study is finite and narrow; no variable-coefficient proof was supplied |
| ADE is generally faster than the 2018 method | **UNTESTED** | No fair implementation-level timing comparison with the stretched-grid 2018 method yet |

## Presentation-safe conclusion

The strongest defensible claim is: *the implementation reproduces the 2025
paper's numerical sequences, the corrected target errors and the spatial
operator are compatible with second-order behavior over the tested ranges,
and ADE is dramatically more robust than Forward Euler on the declared stress
ladder. The experiments do not establish unconditional stability or a fair
speed advantage over the stretched-grid 2018 scheme.*

