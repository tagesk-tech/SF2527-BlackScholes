---
type: project
course: sf2527-numerical-methods-for-differential-equations-i
status: active
created: 2026-11-24
updated: 2026-11-24
deadline: 2026-12-07
github: pending-authentication
source: "[[../../Aristotle/Courses/Period 2/sf2527-numerical-methods-for-differential-equations-i/resources/2026-11-14-black-scholes-fd-project-resource|audited project resource]]"
---

# SF2527 Black--Scholes project

> [!goal]
> Implement and audit the shrinking-grid, two-sweep Saul'yev/ADE method in Wang et al. (2025), reproduce the main numerical results, and defend what the evidence does and does not establish in a 15-minute group presentation.

## Live project state

| Field | Current state |
|---|---|
| Problem | Does the paper's method price call, digital-call, and put options accurately near the target asset price without a prescribed far-field value? |
| Given | The paper, its appendix, the corrected analytical Black--Scholes oracle, the course presentation rubric, and the group's comparison notes. |
| Unknown | Whether the implementation reproduces the tables; whether time and space refinement separately support order two; where ADE remains stable relative to Forward Euler; how the put extension behaves. |
| Conditions | MATLAB code; shrinking right grid; no artificial far-field value; all four members present; at most 15 minutes; submit PDF slides, code, and contribution list. |
| Success test | Tests pass; reproduction and separated-refinement tables are generated; failed cases remain visible; figures support a concise claim with stated limits. |
| Working idea | Treat the closed-form solution as the primary oracle and the paper tables as a benchmark, not as ground truth. |
| Present move | Review the validated claim table with the group, assign slide ownership, and rehearse the method explanation without source code. |
| Check | Put--call parity, target/profile errors, observed orders, domain sensitivity, Forward-Euler stress, and a planted wrong-sign reaction term. |

## Working surfaces

- Group document: [SF2527 Project - Black Scholes](https://docs.google.com/document/d/18yEwZ2_1pxxlADkSNC3Z_lzYYY7B7FqQHKqM57H_4gg/edit)
- Mathematical source note: [[../../Aristotle/Courses/Period 2/sf2527-numerical-methods-for-differential-equations-i/resources/2026-11-14-black-scholes-fd-project-resource|Black--Scholes without far-field boundary conditions]]
- Scope and work packages: [[docs/PROJECT_SCOPE]]
- Local code entry point: `experiments/run_all.m`

## Current evidence

- The paper's call and cash-call price sequences are reproduced to their printed precision.
- The paper's exact call benchmark $3.9123$ is **refuted**; the standard formula gives approximately $3.929276$.
- Corrected call errors have observed temporal orders $2.09$, $2.09$, and $2.02$.
- The centered spatial operator has observed max-norm orders $2.000$, $2.000$, and $2.000$.
- The European put extension satisfies the analytical benchmark, left boundary, and put--call parity.
- ADE stays finite and nonnegative on the declared stress ladder; Forward Euler blows up or becomes non-finite. This does not prove unconditional stability.

Full evidence and limits: [[docs/RESULTS|validated results]].

## Next group decision

Freeze which claims enter the slides and attach one named owner to each presentation segment using [[docs/CONTRIBUTIONS_TEMPLATE|the contribution template]]. No contribution claim should be inferred from the phone-number list or chat ordering.
