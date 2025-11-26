# ⚛️ Quantum Solution for the Instant Insanity Puzzle (ECE/CSC 569 Final Project)

This repository contains the full implementation and documentation for the final project for **ECE/CSC 569**, focusing on solving the classic **Instant Insanity Puzzle** using near-term quantum algorithms.

## 📄 Project Overview

The Instant Insanity Puzzle is a classic combinatorial optimization problem involving stacking four colored cubes such that each of the four vertical sides of the column displays all four colors exactly once.

The goal of this project was to:
1.  Formulate the puzzle's complex constraints (Rules 1-3) into a **Quadratic Unconstrained Binary Optimization (QUBO) Hamiltonian**.
2.  Implement a **Quantum Approximate Optimization Algorithm (QAOA)** and **Variational Quantum Eigensolver (VQE)** hybrid approach to find the ground state solution of the Hamiltonian.
3.  Analyze the impact of classical optimization methods (like COBYLA and Nelder-Mead) on the quantum algorithm's convergence and accuracy.

---

## 🔑 Key Methodology & Findings

### 1. Problem Encoding
The solution leverages a **graph theory approach** to model the puzzle, mapping the cube configurations and constraints into a series of Boolean equations. These equations were then converted into a single **Ising Hamiltonian** (or QUBO), where the ground state corresponds to a valid solution.

### 2. Quantum Algorithm Implementation
We implemented a hybrid classical-quantum approach, utilizing:
* **QAOA/VQE:** To variationally explore the quantum state space and minimize the Hamiltonian expectation value.
* **Classical Optimizers:** Native **COBYLA** and a hybrid **COBYLA + Nelder-Mead** strategy were used to minimize the cost function.

### 3. Core Conclusions
* The **choice of classical minimization approach** significantly impacts the effective utility of QAOA.
* The system's convergence is highly sensitive to the initial parameters, and **weights and penalties** used in the Hamiltonian formulation must be carefully tailored to guide the system to the correct minimum.
* Further work could explore tailored Mixer architectures to limit the search space.

---

## 📂 Repository Contents

| File Name | Description | Key Content |
| :--- | :--- | :--- |
| `Final_Project_hybrid.ipynb` | **Primary Codebase (Jupyter Notebook)** | Contains the Hamiltonian formulation, QAOA/VQE implementation, objective function, and the results from the hybrid optimization runs (D-Wave's Ocean SDK and Qiskit components were utilized). |
| `CSC_ECE569_FinalProject_Presentation.pdf` | **Final Presentation Slides** | Detailed slides outlining the problem, encoding rules, classical modeling, quantum implementation steps, comparison of optimizer results, and final conclusions. |
| `ECE_569_Project_Description.pdf` | **Project Proposal** | Outlines the initial problem description, proposed methodology (QAOA, VQE), and project milestones. |

---

## 💻 Technology Stack

* **Languages:** Python
* **Quantum Frameworks:** Qiskit
* **Optimization/Modeling:** D-Wave Ocean SDK (for BQM formulation), `scipy.optimize` (COBYLA, Nelder-Mead), NumPy.
* **Visualization:** Matplotlib, NetworkX.

---

## 🤝 Team

This project was developed in collaboration with:

* **[Sisir Pynda]**
* **[Martin Carignan]**
