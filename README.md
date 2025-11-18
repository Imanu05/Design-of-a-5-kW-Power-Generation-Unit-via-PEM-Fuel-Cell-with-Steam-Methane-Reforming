
## Overview
This project implements a full MATLAB-based simulation of a **5 kW Proton Exchange Membrane (PEM) fuel cell system** with detailed electrochemical modeling.  
The script performs a multi-parameter search across temperature, active area, current density, and number of cells to identify operating conditions capable of producing **~5 kW** electrical power.

The model also computes voltage, hydrogen consumption, efficiency, and heat loss for each operating point, enabling systematic design and analysis of PEM fuel-cell systems.

## What the Code Does
The simulation evaluates thousands of operating combinations and calculates:

- Cell voltage (activation, ohmic, and concentration losses)
- Stack power output  
- Hydrogen flow rate (mol/s and converted to L/hr)  
- Efficiency based on hydrogen HHV  
- Heat loss from the stack  
- V–I characteristic curves  
- All parameter sets that achieve **5 kW ± 5 W**  


## Adjustable Parameters
The following ranges are swept by the model:

- **Current density (I):** 0.01 → 2 A/cm²  
- **Temperature (T):** 325 → 355 K  
- **Active area (A):** 400 → 500 cm²  
- **Number of cells (No):** 80 → 120  

Step sizes for each range are user-controlled via terminal input.

## Key Outputs
- Voltage vs. current (V–I) plot  
- Power output for each operating point  
- Efficiency and heat-loss profiles  
- Hydrogen flow-rate map  
  


## Tools Used
- **MATLAB**  
  - Multivariable parametric sweep  
  - Electrochemical performance modeling  
  - Data extraction & result filtering  
  - Visualization (V–I curves)  

