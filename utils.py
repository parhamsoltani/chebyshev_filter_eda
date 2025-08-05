"""
Utility functions for filter design and visualization.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import sympy as sp

class FilterUtils:
    """Utility functions for filter operations."""
    
    @staticmethod
    def clean_complex_array(arr, tolerance=1e-10):
        """Clean complex array by removing small imaginary parts."""
        if np.isscalar(arr):
            return arr.real if abs(arr.imag) < tolerance else arr
        
        cleaned = []
        for val in arr:
            if isinstance(val, complex):
                if abs(val.imag) < tolerance:
                    cleaned.append(val.real)
                else:
                    cleaned.append(val)
            else:
                cleaned.append(val)
        return np.array(cleaned)
    
    @staticmethod
    def poly_from_roots(roots):
        """Create polynomial coefficients from roots."""
        poly = [1]
        for root in roots:
            # Multiply by (s - root)
            new_poly = [0] * (len(poly) + 1)
            for i, coeff in enumerate(poly):
                new_poly[i] += coeff
                new_poly[i+1] -= coeff * root
            poly = new_poly
        return np.array(poly)

class PlotUtils:
    """Plotting utilities for filter visualization."""
    
    def __init__(self):
        self.fig_counter = 1
    
    def plot_poles(self, core):
        """Plot pole locations."""
        plt.figure(self.fig_counter)
        self.fig_counter += 1
        
        plt.subplot(2, 1, 1)
        if hasattr(core, 'pk11') and core.pk11 is not None:
            plt.plot(np.real(core.pk11), np.imag(core.pk11), 'x', label='F11 poles')
        if hasattr(core, 'pk12') and core.pk12 is not None:
            plt.plot(np.real(core.pk12), np.imag(core.pk12), 'x', label='F12 poles')
        plt.title('Poles of F1')
        plt.xlabel('Real')
        plt.ylabel('Imaginary')
        plt.grid(True)
        plt.legend()
        
        plt.subplot(2, 1, 2)
        if hasattr(core, 'pk2') and core.pk2 is not None:
            plt.plot(np.real(core.pk2), np.imag(core.pk2), 'x', label='F2 poles')
        plt.title('Poles of F2')
        plt.xlabel('Real')
        plt.ylabel('Imaginary')
        plt.grid(True)
        plt.legend()
        
        plt.tight_layout()
    
    def plot_group_delay(self, core):
        """Plot group delay characteristics."""
        plt.figure(self.fig_counter)
        self.fig_counter += 1
        
        if core.TF_F1 is not None:
            plt.subplot(2, 1, 1)
            w, h = signal.freqs(core.DEN1, core.NUM1, worN=1000)
            group_delay = -np.diff(np.unwrap(np.angle(h))) / np.diff(w)
            plt.semilogx(w[1:], group_delay)
            plt.xlabel('Frequency (rad/s)')
            plt.ylabel('Group delay (s)')
            plt.title('F1 Group Delay')
            plt.grid(True)
        
        if core.TF_F2 is not None:
            plt.subplot(2, 1, 2)
            w, h = signal.freqs(core.DEN2, core.NUM2, worN=1000)
            group_delay = -np.diff(np.unwrap(np.angle(h))) / np.diff(w)
            plt.semilogx(w[1:], group_delay)
            plt.xlabel('Frequency (rad/s)')
            plt.ylabel('Group delay (s)')
            plt.title('F2 Group Delay')
            plt.grid(True)
        
        plt.tight_layout()
    
    def plot_magnitude_response(self, core):
        """Plot magnitude response."""
        plt.figure(self.fig_counter)
        self.fig_counter += 1
        
        if core.TF_F1 is not None:
            plt.subplot(2, 1, 1)
            w = np.logspace(-2, 2, 1000)
            w_eval, h = signal.freqs(core.NUM1, core.DEN1, worN=w)
            plt.semilogx(w_eval, 20*np.log10(np.abs(h)))
            plt.xlabel('Frequency (rad/s)')
            plt.ylabel('Magnitude (dB)')
            plt.title('F1 Magnitude Response')
            plt.grid(True)
        
        if core.TF_F2 is not None:
            plt.subplot(2, 1, 2)
            w = np.logspace(-2, 2, 1000)
            w_eval, h = signal.freqs(core.NUM2, core.DEN2, worN=w)
            plt.semilogx(w_eval, 20*np.log10(np.abs(h)))
            plt.xlabel('Frequency (rad/s)')
            plt.ylabel('Magnitude (dB)')
            plt.title('F2 Magnitude Response')
            plt.grid(True)
        
        plt.tight_layout()