# task-offloading
MATLAB code for the simulation of our paper entitled "**Energy-Efficient Task Offloading Under E2E Latency Constraints**"

How to cite:
M. Tajallifar, S. Ebrahimi,  M. Javan, N. Mokari, L. Chiaraviglio, "Energy-Efficient Task Offloading Under E2E Latency Constraints", IEEE Transactions on Communications, 2021.

Preprint on arXiv:
https://arxiv.org/abs/1912.00187

ResearchGate:
https://www.researchgate.net/publication/337654553_Energy-Efficient_Task_Offloading_Under_E2E_Latency_Constraints

IEEE Xplore:
https://ieeexplore.ieee.org/document/9638479
https://doi.org/10.1109/TCOMM.2021.3132909

IEEE DataPort:
https://ieee-dataport.org/open-access/energy-efficient-task-offloading-under-e2e-latency-constraints
https://doi.org/10.21227/w5tv-yz53


**Abstract**:

In this paper, we propose a novel resource management scheme that jointly allocates the transmit power and computational resources in a centralized radio access network architecture. The network comprises a set of computing nodes to which the requested tasks of different users are offloaded. The optimization problem minimizes the energy consumption of task offloading while takes the end-to-end-latency, i.e., the transmission, execution, and propagation latencies of each task, into account. We aim to allocate the transmit power and computational resources such that the maximum acceptable latency of each task is satisfied. Since the optimization problem is non-convex, we divide it into two sub-problems, one for transmit power allocation and another for task placement and computational resource allocation. Transmit power is allocated via the convex-concave procedure. In addition, a heuristic algorithm is proposed to jointly manage computational resources and task placement. We also propose a feasibility analysis that finds a feasible subset of tasks. Furthermore, a disjoint method that separately allocates the transmit power and the computational resources is proposed as the baseline of comparison. A lower bound on the optimal solution of the optimization problem is also derived based on exhaustive search over task placement decisions and utilizing Karush–Kuhn–Tucker conditions. Simulation results show that the joint method outperforms the disjoint method in terms of acceptance ratio. Simulations also show that the optimality gap of the joint method is less than 5%.


**Published in**: IEEE Transactions on Communications ( Volume: 70, Issue: 3, March 2022)

**Date of Publication**: 06 December 2021



**A typical task offloading example**:
![A typical task offloading example](https://github.com/sinaebrahimi/task-offloading/blob/main/Figures/1a.%20A%20typical%20task%20offloading%20example.png)

**System model**:
![System model](https://github.com/sinaebrahimi/task-offloading/blob/main/Figures/1b.%20System%20model.png)


**Usage Guide**:
**Notes on the simulation files**:

_DTO.m_ simulates the **disjoint task offloading (DTO)** method in the manuscript. This file receives the following parameters as its **inputs**:

	1. Number of single-antenna users, which is equal to the number of tasks, that is, K
	2. Maximum acceptable latency of tasks, that is, T=T_k, \forall k
	3. Ratio of RAN latency to the maximum acceptable latency, that is, T_{RAN}/T
	4. Computational load of each task, that is, L=L_k, \forall k
	5. Data size of each task, that is, D=D_k,  \forall k.

After receiving the parameters, _DTO.m_ **executes** the disjoint method and returns the outputs as in the following:

	1. Acceptance Ratio
	2. Radio Transmission latency of all tasks, that is, T_k^{tx}  \forall k
	3. Propagation latency of all tasks, that is, T_k^{prop}  \forall k
	4. Execution latency of all tasks, that is, T_k^{exe}  \forall k

_JTO.m_ simulates the **joint task offloading (JTO)** method in the manuscript. This file receives the following parameters as its **inputs**:

	1. Number of single-antenna users, which is equal to the number of tasks.
	2. Maximum acceptable latency of tasks, that is, T=T_k,  \forall k
	3. Computational load of each task, that is, L=L_k,  \forall k
	4. Data size of each task, that is, D=D_k,  \forall k.

After receiving the parameters, _JTO.m_ **executes** the disjoint method and returns the outputs as in the following:

	1. Acceptance Ratio
	2. Radio Transmission latency of all tasks, that is, T_k^{tx}  \forall k
	3. Propagation latency of all tasks, that is, T_k^{prop}  \forall k
	4. Execution latency of all tasks, that is, T_k^{exe}  \forall k




It is worth mentioning that the rest of the files function as follows:

	- pathbetweennodes.m includes a function that returns all the paths between two nodes of a graph (Copyright 2014 Kelly Kearney). It is used multiple times in both JTO.m and DTO.m.
	- channel.mat is a file used for initializing the channel in DTO.m.
	- 2021-06-22-TaskOffloading (+Complexity Analysis)-arXiv_v3.pdf: the complete version of the paper (with complexity analysis) submitted in arXiv.
	- 2021-06-23-TaskOffloading- Submitted Version to IEEE TCOM.pdf: the version of the paper submitted to the IEEE TCOM (also available on researchgate).
	- Figures: This directory includes all the figures (in PNG format) in the paper.


**Disclaimer**: We used the MOSEK solver and CVX package to solve all problems. Moreover, all simulation steps (including initialization, admission control mechanisms, and solving ILP problems with MOSEK toolbox) have been implemented in MATLAB software which is widely used to solve resource allocation problems.
