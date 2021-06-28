# task-offloading
MATLAB code for the simulation of our paper entitled "**Energy-Efficient Task Offloading Under E2E Latency Constraints**"

M. Tajallifar, S. Ebrahimi,  M. Javan, N. Mokari, L. Chiaraviglio, "Energy-Efficient Task Offloading Under E2E Latency Constraints", Submitted to IEEE Transactions on Communications, 2021.

Preprint on arXiv:
https://arxiv.org/abs/1912.00187

ResearchGate:
https://www.researchgate.net/publication/337654553_Energy-Efficient_Task_Offloading_Under_E2E_Latency_Constraints

**Abstract**:
In this paper, we propose a novel resource management scheme that jointly allocates the transmit power and computational resources in a centralized radio access network architecture. The network comprises a set of computing nodes to which the requested tasks of different users are offloaded. The optimization problem minimizes the energy consumption of task offloading while takes the end-to-end-latency, i.e., the transmission, execution, and propagation latencies of each task, into account. We aim to allocate the transmit power and computational resources such that the maximum acceptable latency of each task is satisfied. Since the optimization problem is non-convex, we divide it into two sub-problems, one for transmit power allocation and another for task placement and computational resource allocation. Transmit power is allocated via the convex-concave procedure. In addition, a heuristic algorithm is proposed to jointly manage computational resources and task placement. We also propose a feasibility analysis that finds a feasible subset of tasks. Furthermore, a disjoint method that separately allocates the transmit power and the computational resources is proposed as the baseline of comparison. A lower bound on the optimal solution of the optimization problem is also derived based on exhaustive search over task placement decisions and utilizing Karush–Kuhn–Tucker conditions. Simulation results show that the joint method outperforms the disjoint method in terms of acceptance ratio. Simulations also show that the optimality gap of the joint method is less than 5%.

