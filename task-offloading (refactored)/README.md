# Energy-Efficient Task Offloading Under E2E Latency Constraints

MATLAB simulation code for the joint and disjoint task-offloading and resource-allocation framework proposed in:

> M. Tajallifar, S. Ebrahimi, M. R. Javan, N. Mokari, and L. Chiaraviglio, "Energy-Efficient Task Offloading Under E2E Latency Constraints," *IEEE Transactions on Communications*, vol. 70, no. 3, pp. 1711–1725, Mar. 2022.
> DOI: [10.1109/TCOMM.2021.3132909](https://doi.org/10.1109/TCOMM.2021.3132909) · Preprint: [arXiv:1912.00187](https://arxiv.org/abs/1912.00187)

## Overview

The code compares two resource-management strategies for offloading user tasks from a radio access network to a set of NFV-enabled computing nodes, under a per-task end-to-end (E2E) latency constraint (transmission + execution + propagation delay):

- **JTO — Joint Task Offloading**: jointly optimizes transmit power, task placement (which node + route), and computational resource allocation.
- **DTO — Disjoint Task Offloading**: allocates radio resources and computational resources separately; used as the baseline for comparison.

Both are formulated as non-convex optimization problems and solved via a combination of the convex-concave procedure, CVX-based convex sub-problems, and a heuristic task-placement refinement step.

## Repository structure

```
├── src/                     MATLAB source
│   ├── main_1.m             Experiment: acceptance ratio vs. max. acceptable latency
│   ├── main_2.m             Experiment: acceptance ratio vs. RAN latency budget (T_Ratio sweep)
│   ├── main_3.m             Experiment: average delay breakdown vs. task data size
│   ├── JTO_function.m       Joint task offloading + resource allocation
│   ├── DTO_function.m       Disjoint task offloading + resource allocation
│   └── pathbetweennodes.m   Third-party graph helper (see Acknowledgments)
├── figures/                 Output figures (.fig editable / .png, .jpg preview / .eps for LaTeX)
├── paper/                   Author preprint (for reference)
├── LICENSE
├── CITATION.cff
└── README.md
```

## Requirements

- MATLAB (developed/tested with recent releases; update this line with the version you use)
- [CVX](http://cvxr.com/cvx/) with a supported solver (e.g., SDPT3, SeDuMi, MOSEK, or Gurobi), with `cvx_setup` already run
- Parallel Computing Toolbox (used by `parfor`/`parpool` in `main_1.m` and `main_2.m`)

## How to run

1. Install CVX and run `cvx_setup` once, if you haven't already.
2. Add `src/` to your MATLAB path, or `cd` into it.
3. Run `main_1.m`, `main_2.m`, or `main_3.m` directly — each script sets its own system/network parameters, calls `JTO_function`/`DTO_function`, and plots the corresponding figure from the paper.

Each script is self-contained: parameters such as the number of users, nodes, antennas, and the network adjacency matrix are defined at the top and can be edited to try other network configurations.

## Acknowledgments

`pathbetweennodes.m` is a third-party graph utility by Kelly Kearney (2014), originally distributed via MATLAB File Exchange, and is included unmodified. It is **not** covered by this repository's license — see the terms noted in the file's own header, and check MATLAB File Exchange for the current license before reusing it elsewhere.

## Citation

If this code is useful in your own work, please cite the paper above (a ready-to-use entry is in `CITATION.cff` — GitHub shows a "Cite this repository" button automatically once that file is present).

## License

Except for the third-party file noted above, this project is released under the MIT License — see `LICENSE`.
