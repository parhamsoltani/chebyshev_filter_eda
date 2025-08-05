#!/usr/bin/env python3
"""
Chebyshev Filter Design Tool
A comprehensive tool for designing Chebyshev filters with multiple synthesis methods.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import sympy as sp
from sympy import symbols, sinh, cosh, asinh, sin, cos, pi, I, solve, limit, simplify
from sympy import poly, Poly, expand, collect, cancel, apart

from filter_core import ChebyshevCore
from synthesis_methods import (
    TransmissionMatrixSynthesis, 
    YanagisawaSynthesis,
    MathewsSeifertSynthesis,
    LoveringSynthesis,
    MitraSynthesis,
    KuhSynthesis,
    StateVariableSynthesis
)
from utils import FilterUtils, PlotUtils
from cauer_synthesis import CauerSynthesis

def main():
    """Main function for Chebyshev filter design."""
    
    # Input parameters
    N = int(input("Enter N: "))
    if N <= 1:
        print("N must be equal or greater than 2")
        return
        
    e1 = float(input("Enter e1: "))
    e2 = float(input("Enter e2: "))
    dc = int(input("Do you want the design to be DC-drop-free? Press 1 for yes, 0 for no: "))
    
    Rg_denorm = float(input("Enter Rg: "))
    RL_denorm = float(input("Enter RL: "))
    
    # Initialize core filter design
    core = ChebyshevCore(N, e1, e2, Rg_denorm, RL_denorm)
    
    # Design filters based on DC-drop preference
    if dc == 0 or N % 2 != 0:
        F1_tf, F2_tf = core.design_with_dc_drop()
    else:
        F1_tf, F2_tf = core.design_without_dc_drop()
    
    # Plot results
    plotter = PlotUtils()
    plotter.plot_poles(core)
    plotter.plot_group_delay(core)
    plotter.plot_magnitude_response(core)
    
    # Apply synthesis methods
    synthesizer = TransmissionMatrixSynthesis(core)
    T1, T2 = synthesizer.compute_transmission_matrices()
    L1_vals, C1_vals, L2_vals, C2_vals = synthesizer.decompose_matrices(T1, T2)
    
    # Apply other synthesis methods
    yanagisawa = YanagisawaSynthesis(core)
    ya1_y, yb1_y, yc1_y, yd1_y, ya2_y, yb2_y, yc2_y, yd2_y = yanagisawa.synthesize()
    
    # Cauer synthesis
    cauer = CauerSynthesis()
    if hasattr(core, 'Zin1'):
        L_vals1, C_vals1 = cauer.cauer_synthesis(core.Zin1)
        print(f"Cauer synthesis F1 - L values: {L_vals1}")
        print(f"Cauer synthesis F1 - C values: {C_vals1}")
    
    # Display results
    print(f"Transmission matrix decomposition:")
    print(f"L1 values: {L1_vals}")
    print(f"C1 values: {C1_vals}")
    print(f"L2 values: {L2_vals}")
    print(f"C2 values: {C2_vals}")
    
    plt.show()

if __name__ == "__main__":
    main()