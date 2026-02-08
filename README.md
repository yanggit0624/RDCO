# Receptor-driven cycloidal optimization (RDCO)

This repository contains the official MATLAB source code for the **Receptor-driven cycloidal optimization (RDCO)** algorithm.

This work has been published in **Applied Mathematical Modelling**.

> **Paper Title:** Receptor-driven cycloidal optimization: Integrating biological control models with geometric search trajectories for engineering design optimization  
> **Journal:** Applied Mathematical Modelling  
> **DOI:** [10.1016/j.apm.2026.116827]
---

## Abstract

The optimization of complex engineering systems presents persistent computational challenges involving high-dimensional spaces, nonlinear objectives, and multiple local optima. Metaheuristic algorithms have become essential tools for such problems, with their effectiveness fundamentally dependent on the balance between exploration and exploitation. However, many existing approaches employ control mechanisms with limited responsiveness to search state changes and operators lacking functional specialization. This paper introduces the **Receptor-Driven Cycloidal Optimization (RDCO)** algorithm to address these limitations. The algorithm integrates a Hill-equation-based control mechanism that enables threshold-sensitive behavioral switching, with cycloidal operators where bounded hypocycloids serve exploitation and expansive epicycloids serve exploration. Experimental evaluation on CEC2017 and CEC2022 benchmarks demonstrates competitive performance, with advantages most evident on multimodal and composition functions. Engineering design applications validate practical utility.

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

f you use this code or the RDCO algorithm in your research, please cite our paper:
"Lin Yang and Hongwen Xu, "Receptor-driven cycloidal optimization: Integrating biological control models with geometric search trajectories for engineering design optimization," Applied Mathematical Modelling, 2026. DOI: 10.1016/j.apm.2026.116827.

---

## Contact

For any questions, please contact Lin Yang at 15765057075@163.com.