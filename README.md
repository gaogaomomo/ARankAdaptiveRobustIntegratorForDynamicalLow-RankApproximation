# Rank-Adaptive Robust Integrator for Dynamical Low-Rank Approximation

This repository contains the **MATLAB implementation** of the **rank-adaptive robust integrator** used in the numerical experiments of the thesis:  
**"A Rank-Adaptive Robust Integrator for Dynamical Low-Rank Approximation"**.

## 📌 Overview
This MATLAB code implements a **rank-adaptive integrator** that extends the robust fixed-rank integrator (Ceruti-Lubich 2021) by dynamically adjusting the rank based on the evolving structure of the matrix or tensor.

## 📂 Files
- `dlr` – Core implementation of the rank-adaptive integrator for discrete Schrödinger equation.
- `RK4` – Function of Runge-Kutta method.
- `err` - Computation for error estimation.
- `mainDLR` – Example script demonstrating how to use the integrator.

## 🛠️ Dependencies
- MATLAB R2023a
- Required toolboxes:
  - **MATLAB's built-in linear algebra functions** (`svd`, `qr`, etc.)

## 🚀 Usage
Run the following script to test the integrator on an example problem:

```matlab
mainDLR

