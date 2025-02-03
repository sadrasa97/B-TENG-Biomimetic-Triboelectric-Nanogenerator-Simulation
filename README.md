# B-TENG (Biomimetic Triboelectric Nanogenerator) Simulation

## Overview
This project provides a MATLAB-based simulation of the **B-TENG (Biomimetic Triboelectric Nanogenerator)**. The simulation includes:
- **Energy harvesting calculations** based on governing equations.
- **Fluid–structure interaction and nonlinear film dynamics** modeling.
- **Electrical output performance analysis** in relation to wind velocity and load resistance.
- **Graphical representations** of dynamic behaviors and energy conversion efficiency.

## Features
- **Energy Harvesting Analysis** (Equations 7–11)
- **Fluid-Structure Interaction (FSI) Modeling** (Equations 1–6)
- **Electrical Output Simulation**
- **FFT-based Frequency Analysis**
- **MATLAB-based Visualizations** for:
  - The Bernoulli effect and flapping film interactions.
  - Electric potential distribution vs. contact area variation.
  - Output performance metrics: Voltage, Current, Power, and Efficiency.

## Equations Implemented
### **Energy Harvesting Calculations**
1. Volume Flow Rate: \( Q_v = v \times S \)
2. Mass Flow Rate: \( Q_m = Q_v \times \rho \)
3. Instantaneous Output Power: \( P_o = \frac{V^2}{R} \)
4. Instantaneous Input Power: \( P_i = \frac{Q_m \times v^2}{2} \)
5. Conversion Efficiency: \( \eta = \frac{P_o}{P_i} \times 100 \%

### **Fluid–Structure Interaction & Film Dynamics**
- Governing equations model **geometrically nonlinear deformation** of the flapping film.
- Includes **tensile stress, flexural rigidity, and immersed boundary effects**.

### **Electrical Output Analysis**
- Simulates **voltage, current, and power output** based on wind velocity.
- Computes **efficiency as a function of input power and harvested power**.

## Installation & Usage
### **Requirements**
- MATLAB (Recommended: R2021a or later)

### **Running the Simulation**
1. Clone the repository:
   ```sh
   git clone https://github.com/yourusername/bteng-simulation.git
   cd bteng-simulation
   ```
2. Open MATLAB and navigate to the project folder.
3. Run the main script:
   ```matlab
   run('BTENG.m')
   ```

## Results & Outputs
The simulation generates plots for:
- **Bernoulli effect and film motion**
- **Electrical output vs. wind velocity**
- **Frequency response analysis using FFT**
- **Energy harvesting efficiency**

## Contributions
Feel free to contribute by submitting pull requests or reporting issues!


---
**Author:** Sadra Saremi

