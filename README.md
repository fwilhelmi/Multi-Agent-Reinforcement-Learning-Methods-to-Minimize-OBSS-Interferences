# Multi-Agent-Reinforcement-Learning-Methods-to-Minimize-OBSS-Interferences
Code used for the master thesis at MIIS (UPF)

In each of the sections it is included the code for performing the executions based on a very simple simulation of the CSMA/CA protocol.

1- OBSS Effects: it is compared the overall throughput and fairness of a given scenario according to the number of coexistent WLANs and one variable between number of channels, CCA level or TPC level.

2- OBSS Happiness: it is computed the overall time required to reach an overall compliance status by changing the channel of unsatisfied WLANs at each iteration.

3- Selfish actions: it is computed the overall throughput and fairness obtained by letting each WLAN is an OBSS to choose the sequentially the best configuration for its own benefit.

4- Decentralized Q-learning (local): it is applied the Q-learning algorithm to each WLAN in parallel with the goal of improving the overall throughput and fairness by choosing a better configuration. The reward in this case is the throughput experienced by the WLAN.

5- Decentralized Q-learning (shared): the same approach as (4) but reward in this case is shared, which is the overall throughput.

6- MA-MAB (e-greedy): it is applied the Multi Armed Bandits approach with the goal of improving the overall throughput and fairness by choosing a better configuration (e-greedy variation)

7- MA-MAB (EXP3): the same as 6, but using EXP3 variation.

8- MA-MAB (UCB): the same as 6, but using UCB variation.

9- Centralized Correlated Equilibria: it is applied a brute force algorithm to find the best possible configuration given an OBSS map that maximizes overall throughput * fairness