# Chebyshev Filter Design Tool

A comprehensive Python-based tool for designing Chebyshev filters with multiple synthesis methods. This tool provides various circuit realization techniques for implementing Chebyshev filters in both passive (LC) and active configurations.

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Parameter Guide](#parameter-guide)
- [Synthesis Methods](#synthesis-methods)
- [Usage Examples](#usage-examples)
- [Output Interpretation](#output-interpretation)
- [Module Structure](#module-structure)
- [Mathematical Background](#mathematical-background)
- [Contributing](#contributing)
- [License](#license)

## Overview

This tool implements various classical and modern synthesis methods for Chebyshev filter design, originally based on MATLAB implementations. It's particularly useful for:

- RF/Microwave filter design
- Analog circuit design education
- Research in filter synthesis methods
- Comparative analysis of different realization techniques

## Features

### Filter Design Capabilities
- **Chebyshev Type I filters** with specified ripple characteristics
- **Dual-section filter design** (F1 and F2 sections)
- **DC-drop and DC-drop-free configurations**
- **Arbitrary filter orders** (N ≥ 2)
- **Flexible impedance matching** (arbitrary source and load resistances)

### Synthesis Methods
- **Transmission Matrix Synthesis** - ABCD parameter decomposition
- **Yanagisawa Method** - Partial fraction decomposition
- **Mathews-Seifert Method** - Alternative decomposition approach
- **Lovering Method** - Simplified synthesis technique
- **Mitra Method** - Active filter realization
- **Kuh Method** - Active circuit synthesis
- **Cauer Synthesis** - LC ladder network realization
- **State Variable Method** - Active filter implementation

### Visualization
- **Pole-zero plots** for both filter sections
- **Group delay characteristics**
- **Magnitude response plots**
- **Interactive plotting** with matplotlib

## Installation

### Prerequisites
```bash
# Required Python version
Python 3.7+
```

### Install Dependencies
```bash
pip install -r requirements.txt
```

### Required Packages
```
numpy>=1.20.0
scipy>=1.7.0
matplotlib>=3.3.0
sympy>=1.8.0
```

### Setup
```bash
git clone <repository-url>
cd chebyshev-filter-design
python chebyshev_filter.py
```

## Quick Start

### Basic Usage
```bash
python chebyshev_filter.py
```

The program will prompt you for the following parameters:

```
Enter N: 4
Enter e1: 0.5
Enter e2: 0.5
Do you want the design to be DC-drop-free? Press 1 for yes, 0 for no: 0
Enter Rg: 50
Enter RL: 50
```

This creates a 4th-order Chebyshev filter with moderate ripple in both sections.

## Parameter Guide

### N (Filter Order)
- **Range**: Integer ≥ 2
- **Description**: Determines the steepness of the filter transition
- **Common values**: 3, 4, 5, 6, 7, 8
- **Impact**: 
  - Higher N → Steeper rolloff, more poles
  - Lower N → Gentler rolloff, simpler circuit

### e1 (First Section Ripple Factor)
- **Range**: Positive real number
- **Description**: Controls passband ripple in the first filter section
- **Relationship**: Ripple (dB) = 10 × log₁₀(1 + e1²)
- **Common values**:
  - `e1 = 0.1` → 0.043 dB ripple (very low ripple)
  - `e1 = 0.5` → 0.97 dB ripple (moderate ripple)
  - `e1 = 1.0` → 3.01 dB ripple (high ripple)

### e2 (Second Section Ripple Factor)
- **Range**: Positive real number
- **Description**: Controls passband ripple in the second filter section
- **Usage**: Often set equal to e1 for symmetric design

### DC-drop Configuration
- **0**: Allow DC drop (traditional Chebyshev response)
- **1**: DC-drop-free design (available only for even N)
- **Impact**: Affects the low-frequency behavior of the filter

### Rg (Source Resistance)
- **Range**: Positive real number (Ohms)
- **Description**: Generator/source resistance
- **Common values**: 50Ω, 75Ω, 300Ω

### RL (Load Resistance)
- **Range**: Positive real number (Ohms)
- **Description**: Load/termination resistance
- **Note**: Can be different from Rg for impedance transformation

## Synthesis Methods

### 1. Transmission Matrix Synthesis
Decomposes the filter into cascade of unit elements with ABCD parameters.
```python
synthesizer = TransmissionMatrixSynthesis(core)
T1, T2 = synthesizer.compute_transmission_matrices()
L1_vals, C1_vals, L2_vals, C2_vals = synthesizer.decompose_matrices(T1, T2)
```

### 2. Cauer Synthesis (LC Ladders)
Realizes the filter as a ladder network of inductors and capacitors.
```python
cauer = CauerSynthesis()
L_values, C_values = cauer.cauer_synthesis(impedance_function)
```

### 3. Active Synthesis Methods
Various methods for op-amp based implementations:
- **Kuh Method**: Admittance matrix approach
- **State Variable**: Using integrators and summers
- **Mitra Method**: Alternative active realization

## Usage Examples

### Example 1: Low-Ripple Bandpass Filter
```python
# Design parameters for low-ripple filter
N = 5           # 5th order
e1 = 0.1        # Low ripple in first section
e2 = 0.1        # Low ripple in second section  
dc = 0          # Allow DC drop
Rg = 50         # 50Ω source
RL = 50         # 50Ω load
```

### Example 2: High-Selectivity Filter
```python
# Design parameters for sharp cutoff
N = 8           # 8th order for sharp transition
e1 = 1.0        # Higher ripple acceptable for sharper response
e2 = 1.0        
dc = 1          # DC-drop-free (if N is even)
Rg = 75         # 75Ω source
RL = 50         # 50Ω load (impedance transformation)
```

### Example 3: Educational/Analysis
```python
# Compare different synthesis methods
N = 4
e1 = 0.5
e2 = 0.5

# This will generate:
# - Pole-zero plots
# - Group delay plots  
# - Magnitude response
# - Component values for multiple synthesis methods
```

## Output Interpretation

### Plots Generated
1. **Pole Locations**: Shows filter stability and frequency response characteristics
2. **Group Delay**: Indicates phase linearity (important for signal integrity)
3. **Magnitude Response**: Shows passband ripple and stopband attenuation

### Component Values
The tool outputs component values for various realizations:
```
Transmission matrix decomposition:
L1 values: [0.123, 0.456, 0.789, ...]  # Inductor values in Henries
C1 values: [1.23e-9, 4.56e-9, ...]     # Capacitor values in Farads
L2 values: [0.234, 0.567, ...]
C2 values: [2.34e-9, 5.67e-9, ...]

Cauer synthesis F1 - L values: [0.145, 0.234, ...]
Cauer synthesis F1 - C values: [1.45e-9, 2.34e-9, ...]
```

### Understanding the Results
- **L values**: Inductance in Henries (for practical circuits, often mH or µH)
- **C values**: Capacitance in Farads (for practical circuits, often nF or pF)
- **Multiple methods**: Different synthesis approaches may give different component values for the same response

## Module Structure

```
chebyshev-filter-design/
│
├── chebyshev_filter.py      # Main application entry point
├── filter_core.py           # Core filter design algorithms
├── synthesis_methods.py     # Various synthesis implementations
├── cauer_synthesis.py       # LC ladder synthesis
├── utils.py                 # Utility functions and plotting
├── requirements.txt         # Python dependencies
└── README.md               # This file
```

### Key Classes
- **ChebyshevCore**: Main filter design engine
- **TransmissionMatrixSynthesis**: ABCD parameter method
- **CauerSynthesis**: LC ladder realization
- **PlotUtils**: Visualization functions

## Mathematical Background

### Chebyshev Polynomials
The tool uses Chebyshev polynomials Tₙ(x) defined by the recurrence:
- T₀(x) = 1
- T₁(x) = x  
- Tₙ(x) = 2x·Tₙ₋₁(x) - Tₙ₋₂(x)

### Filter Transfer Function
For a Chebyshev filter of order N:
```
|H(jω)|² = 1 / (1 + ε²Tₙ²(ω/ωc))
```
Where ε is the ripple factor and ωc is the cutoff frequency.

### Pole Locations
Chebyshev filter poles are located on an ellipse in the s-plane:
```
poles = -σₖ ± jωₖ
σₖ = sinh(β) sin(θₖ)
ωₖ = cosh(β) cos(θₖ)
β = (1/N) asinh(1/ε)
θₖ = (2k-1)π/(2N)
```

## Troubleshooting

### Common Issues

#### "N must be equal or greater than 2"
- Ensure filter order N ≥ 2
- Single-pole filters are not supported

#### "DC-drop-free design only available for even N"
- DC-drop-free option requires even filter order
- Set dc=0 for odd filter orders

#### "Function not suitable for Cauer synthesis"
- Some transfer functions cannot be realized as LC ladders
- Try different synthesis methods
- Check for negative component values

#### Poor numerical accuracy
- Reduce filter order for very high N
- Check input parameter ranges
- Ensure e1, e2 are reasonable values

### Performance Tips
- For high-order filters (N > 10), expect longer computation times
- Very small ripple factors (e < 0.01) may cause numerical issues
- Large impedance ratios (Rg/RL >> 10) may affect synthesis accuracy

## Contributing

Contributions are welcome. Areas for improvement:
- Additional synthesis methods
- GUI interface
- Optimization routines
- More filter types (Butterworth, Elliptic)
- Improved numerical stability

### Development Setup
```bash
git clone <repository-url>
cd chebyshev-filter-design
pip install -r requirements.txt
# Make your changes
# Test with various parameter combinations
# Submit pull request
```

## Contact

[parham.soltany@gmail.com]

---

*This tool is designed for educational and research purposes. For production applications, verify all component values through simulation and measurement.*
