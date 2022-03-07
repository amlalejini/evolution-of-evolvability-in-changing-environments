# 2022-03-02 - Costly copy instruction exploratory experiment

From other experiments, we know that mutation rates evolve toward zero in constant environments.
In this experiment, we impose a cost on high-fidelity h-copy instruction operators. I.e., the lower the mutation rate, the higher the cost.

NOTE: These populations were not well-mixed (i.e., offspring were placed in local neighborhood instead of randomly in the population)

Conditions

- Cost of copy fidelity (6): 0.3, 0.1, 0.03, 0.01, 0.003, 0.0
- Initial mutation rate (3):  0.0001, 0.001, 0.0316
- Environment (2): const-a, const-b