# Receptor-Driven Cycloidal Optimization (RDCO)

This repository contains the official MATLAB source code for the **Receptor-Driven Cycloidal Optimization (RDCO)** algorithm, corresponding to the paper:
> **"Receptor-Driven Cycloidal Optimization (RDCO): A Novel Bio-Inspired Algorithm with Dynamic Geometric Search"**
> 
> *Currently under review.*

---

## Abstract

Effectively managing the trade-off between exploration and exploitation remains a critical challenge in metaheuristic algorithm design. Many existing algorithms rely on pre-scheduled control mechanisms and elementary search operators, which can limit their adaptability on complex fitness landscapes. To address these limitations, this paper introduces the Receptor-Driven Cycloidal Optimization (RDCO) algorithm. RDCO is based on a two-fold design: (1) a responsive control mechanism, derived from the Hill equation, that dynamically regulates the exploration-exploitation balance based on agent performance; and (2) a pair of purpose-built geometric search operators—epicycloids for global exploration and hypocycloids for local exploitation. The performance of RDCO is comprehensively evaluated on the CEC2017 and CEC2022 benchmark suites and three real-world constrained engineering problems, benchmarked against nine established and contemporary metaheuristics. Statistical analysis shows that RDCO consistently achieves the first overall rank in the Friedman tests across all evaluated scenarios. By integrating a biological feedback model with non-linear kinematics, RDCO provides an effective framework for solving complex optimization problems.

---

## How to Run the Experiment

To replicate the benchmark results from the paper:
1.  Ensure all files from this repository are in the same directory.
2.  Open MATLAB (tested on R2024a and newer).
3.  Run the `main.m` script.
    ```matlab
    >> main
    ```
---

## Citation

Full citation details will be provided here upon publication.

---

## Contact

For any questions, please contact Lin Yang at 15765057075@163.com.