# task-offloading
MATLAB code for the simulation of our paper entitled "**Energy-Efficient Task Offloading Under E2E Latency Constraints**"

M. Tajallifar, S. Ebrahimi,  M. Javan, N. Mokari, L. Chiaraviglio, "Energy-Efficient Task Offloading Under E2E Latency Constraints", Submitted to IEEE Transactions on Communications, 2021.

Preprint on arXiv:
https://arxiv.org/abs/1912.00187

ResearchGate:
https://www.researchgate.net/publication/337654553_Energy-Efficient_Task_Offloading_Under_E2E_Latency_Constraints

**Abstract**:

In this paper, we propose a novel resource management scheme that jointly allocates the transmit power and computational resources in a centralized radio access network architecture. The network comprises a set of computing nodes to which the requested tasks of different users are offloaded. The optimization problem minimizes the energy consumption of task offloading while takes the end-to-end-latency, i.e., the transmission, execution, and propagation latencies of each task, into account. We aim to allocate the transmit power and computational resources such that the maximum acceptable latency of each task is satisfied. Since the optimization problem is non-convex, we divide it into two sub-problems, one for transmit power allocation and another for task placement and computational resource allocation. Transmit power is allocated via the convex-concave procedure. In addition, a heuristic algorithm is proposed to jointly manage computational resources and task placement. We also propose a feasibility analysis that finds a feasible subset of tasks. Furthermore, a disjoint method that separately allocates the transmit power and the computational resources is proposed as the baseline of comparison. A lower bound on the optimal solution of the optimization problem is also derived based on exhaustive search over task placement decisions and utilizing Karush–Kuhn–Tucker conditions. Simulation results show that the joint method outperforms the disjoint method in terms of acceptance ratio. Simulations also show that the optimality gap of the joint method is less than 5%.

**A typical task offloading example**:
![A typical task offloading example](https://github.com/sinaebrahimi/task-offloading/blob/main/Figures/1a.%20A%20typical%20task%20offloading%20example.png)

**System model**:
![System model](https://github.com/sinaebrahimi/task-offloading/blob/main/Figures/1b.%20System%20model.png)


**Usage Guide**:
**
Notes on the simulation files**:
DTO.m simulates the disjoint task offloading (DTO) method in the manuscript. This file receives the following parameters as its inputs:

	1. Number of single-antenna users, which is equal to the number of tasks, that is, K
	2. Maximum acceptable latency of tasks, that is, T=T_k, \forall k
	3. Ratio of RAN latency to the maximum acceptable latency, that is, T_{RAN}/T
	4. Computational load of each task, that is, L=L_k, \forall k
	5. Data size of each task, that is, D=D_k,  \forall k.

After receiving the parameters, DTO.m executes the disjoint method and returns the outputs as in the following:

	1. Acceptance Ratio
	2. Radio Transmission latency of all tasks, that is, T_k^{tx}  \forall k
	3. Propagation latency of all tasks, that is, T_k^{prop}  \forall k
	4. Execution latency of all tasks, that is, T_k^{exe}  \forall k

JTO.m simulates the joint task offloading (JTO) method in the manuscript. This file receives the following parameters as its inputs:

	1. Number of single-antenna users, which is equal to the number of tasks.
	2. Maximum acceptable latency of tasks, that is, T=T_k,  \forall k
	3. Computational load of each task, that is, L=L_k,  \forall k
	4. Data size of each task, that is, D=D_k,  \forall k.

After receiving the parameters, JTO.m executes the disjoint method and returns the outputs as in the following:

	1. Acceptance Ratio
	2. Radio Transmission latency of all tasks, that is, T_k^{tx}  \forall k
	3. Propagation latency of all tasks, that is, T_k^{prop}  \forall k
	4. Execution latency of all tasks, that is, T_k^{exe}  \forall k


Disclaimer: We used the MOSEK solver and CVX package to solve all problems. Moreover, all simulation steps (including initialization, admission control mechanisms, and solving ILP problems with MOSEK toolbox) have been implemented in MATLAB software which is widely used to solve resource allocation problems.
