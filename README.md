# YADE Bearing Capacity Test Simulation

This repository contains a simple YADE (Yet Another Dynamic Engine) Python script to perform a preliminary plate load test simulation on borehole soils BH-5 (granular) and BH-6 (clayey/mixed). The goal is to validate approximate bearing capacity and settlement trends based on geotechnical borehole data.

---

## Objective

- Simulate a standard plate load test on a simplified soil domain.
- Estimate bearing capacity and settlement for two soil types:
  - **BH-5**: Granular soil with high stiffness and friction.
  - **BH-6**: Clayey/mixed soil with lower stiffness and higher compressibility.
- Support foundation design recommendations (shallow vs deep).

---

## Files

- `bearing_capacity_test.py`: YADE script to run the numerical simulation.
- Output files after running:
  - `bearing_capacity_results.txt`: Recorded displacement data.
  - `bearing_capacity_plot.png`: Plot of plate displacement vs simulation steps.

---

## Software Requirements

- [YADE DEM Simulator](https://yade-dem.org/) installed on Linux.

### Install YADE on Ubuntu/Debian:

```bash
sudo apt-get update
sudo apt-get install yade
