# Dynamic Analysis of the CityWave Canopy

Experimental Modal Analysis (EMA) and Finite Element Model (FEM) validation of a scaled wind-tunnel model of the **CityWave Canopy** structure (CityLife district, Milan).

* **Institution:** Politecnico di Milano 
* **Program:** M.Sc. in Mechanical Engineering (Mechatronics and Robotics)
* **Academic Year:** 2025/2026
* **Authors:** Caterina Brunetti, Alessia Cecchini, Federico Della Vittoria, Tareq El Kajjyouy

---

## Project Overview

The objective of this project is to identify the modal parameters (natural frequencies, non-dimensional damping ratios, and mode shapes) of the CityWave canopy from noisy experimental wind-tunnel measurements and validate them against a numerical FEM solution.

* **Experimental Setup:** Swept-sine excitation via modal shaker; dynamic response measured across **202 points** using a Laser Doppler Vibrometer (LDV).
* **Frequency Range of Interest:** $2 - 8\text{ Hz}$ (wind-tunnel scale), subsequently scaled to full-scale structural frequencies.

---

## Methodology & Key Results

1. **Eigenfrequency Identification:** Computed the spatial average of the magnitude across all 202 Frequency Response Functions ($H_1$) to isolate the first 4 structural modes despite high experimental noise.
2. **Damping Ratio Estimation ($\xi$):** Filtered anti-nodal sensors using an **80% peak-amplitude threshold** to exclude nodal/faulty channels, then applied the **Half-Power Bandwidth (Peak-Picking) method**:
   $$\xi_i = \frac{\omega_2^2 - \omega_1^2}{4\omega_i^2}$$
3. **Mode Shapes & Spatial Interpolation:** Extracted experimental mode shapes from the imaginary part of the FRFs and interpolated them onto the denser FEM mesh nodes.
4. **FEM vs. EMA Correlation (MAC & CoMAC):** Quantified modal correlation globally via the **Modal Assurance Criterion (MAC)** and locally at each degree of freedom via **CoMAC**.


## Repository Structure

```text
├── docs/
│   └── Dynamic Analysis of the CityWave Canopy _ Tareq El Kajjyouy.pdf  # Full technical report
├── src/
│   ├── main.m              # Main script: FRF processing, EMA identification, interpolation, MAC & CoMAC
│   ├── find_sensors.m      # Helper function: anti-nodal sensor filtering (80% peak-amplitude threshold)
│   ├── peak_picking.m      # Helper function: damping ratio estimation via Half-Power Bandwidth method
│   ├── students.m          # Visualization of numerical FEM mode shapes
│   ├── connectivity.mat    # FEM mesh node coordinates and surface connectivity
│   ├── modal_output.mat    # Numerical mode shapes extracted from FEM analysis
│   ├── punti_laser.mat     # Spatial coordinates of the 202 LDV measurement points
│   └── FRF_H1.mat          # Experimental Frequency Response Functions (202 sensors)
├── LICENSE
└── README.md
