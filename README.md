## Overview
This repository contains MATLAB(2022b) scripts and data for the analysis of heat exchangers. It focuses on system identification using an evolving neuro-fuzzy system and fault detection through prediction intervals. The approach is demonstrated on a real plate heat exchanger with various faults introduced for analysis purposes.

## Contents
- **MATLAB Scripts**: These scripts handle tasks such as system identification, fault detection, data acquisition, and more.
- **Data Files**: Includes `.mat` files with both synthetic (PRBS and staircase excitation) and real-world data for system identification and fault detection.
- **LICENSE**: Specifies the terms under which this project can be used.

## Key Scripts
- `heat_exchanger_prediction_intervals_main.m`: The main script of the program.
- Other scripts (`heat_exchanger_prediction_intervals_*.m`): Various utilities for data processing, system identification, and fault analysis.

## Data Description
- **Identification Data**: Synthetic data (PRBS and staircase excitation) for system identification.
- **Fault Detection Data**: Staircase fault detection scenarios and real-world data (`measurement_steps_all_*.mat`) for testing and analysis.

## Getting Started
Clone the repository and run the `heat_exchanger_prediction_intervals_main.m` script in MATLAB. Ensure that MATLAB is up to date and that all necessary dependencies are installed.
