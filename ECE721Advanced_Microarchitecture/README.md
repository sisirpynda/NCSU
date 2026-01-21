# RISC-V Out-of-Order Stride Value Predictor

## Project Overview
This project involved the implementation of a cycle-accurate **Stride Value Predictor (SVP)** and its associated microarchitectural support structures within the **721sim** (Professor's) simulator. The primary goal was to improve processor performance (IPC) by breaking data dependencies between producer and consumer instructions at the **Dispatch** stage using speculative execution.

The project was developed for **ECE 721: Advanced Microarchitecture** at North Carolina State University under the direction of **Prof. Eric Rotenberg**.

## Technical Implementation

### 1. Stride Value Predictor (SVP) Logic
I engineered a hardware-accurate SVP to predict the results of Integer ALU, Floating-Point ALU, and Load instructions.
* **Prediction Calculation:** The predictor uses a stride-based formula: `prediction = retired_value + (instance * stride)`.
* **Confidence Mechanism:** Implemented saturating confidence counters to ensure that only "confident" predictions are injected into the pipeline, reducing the frequency of costly pipeline squashes.
* **Tagging & Aliasing:** Developed support for full, partial, and no-tag configurations to explore the trade-offs between hardware storage and prediction aliasing.

### 2. Value Prediction Queue (VPQ)
A critical component of the architecture was the **Value Prediction Queue**, which manages the speculative state of in-flight instructions.
* **In-Order Training:** The VPQ facilitates non-speculative training of the SVP at the retirement stage.
* **Speculative Repair:** Implemented a "backward walk" algorithm within the VPQ to repair instance counters during pipeline rollbacks (e.g., branch mispredictions or load violations), ensuring the predictor remains synchronized with the architectural state.
* **Structural Hazard Modeling:** Modeled the VPQ as a structural hazard in the **Rename** stage, where the pipeline stalls if insufficient VPQ entries are available for the current rename bundle.

### 3. Pipeline Integration & Recovery
* **Data Dependency Breaking:** For confident predictions, the destination physical register's "ready bit" is set early in the Dispatch stage, allowing dependent instructions to issue speculatively.
* **Writeback Validation:** Implemented a verification step in the **Writeback** stage to compare computed values against predicted values.
* **Speculative Recovery:** Integrated with the simulator's **VR-1** recovery mechanism to perform a full-pipeline squash and restart execution upon detecting a value misprediction.
* **PRF Optimization:** Optimized Physical Register File (PRF) bandwidth by implementing logic to avoid redundant writes for correctly predicted values.


## Tools & Technologies
* **Language:** C++
* **Simulator:** 721sim (Microarchitecture Simulator)
* **ISA:** RISC-V
* **Concepts:** Speculative Execution, Stride Prediction, Out-of-Order Execution, Design Space Exploration.

## Filed edited
* dispatch.cc
* renamer.cc
* register_read.cc
* writeback.cc
* rename.cc
* retire.cc
* renamer.h
* squash.h

---

### Academic Integrity & Source Code
*In accordance with the academic integrity policy for ECE 721, the source code for this project is hosted in a [private repository](https://github.com/sisirpynda/NCSU_pvt/tree/3726cca04792b3f31e468a503d476219179cde63). The instructor prohibits the public distribution of project code to prevent unauthorized collaboration and protect copyrighted simulator frameworks. Access to the private repository can be provided to recruiters or hiring managers upon request.*
