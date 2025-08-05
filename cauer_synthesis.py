"""
Cauer synthesis for LC ladder networks.
"""

import numpy as np
from scipy import signal

class CauerSynthesis:
    """Cauer synthesis implementation."""
    
    def cauer_synthesis(self, impedance_func):
        """
        Perform Cauer synthesis on an impedance function.
        
        Args:
            impedance_func: Impedance function as transfer function
            
        Returns:
            L_values, C_values: Lists of inductor and capacitor values
        """
        if hasattr(impedance_func, 'num') and hasattr(impedance_func, 'den'):
            num = np.array(impedance_func.num).flatten()
            den = np.array(impedance_func.den).flatten()
        else:
            # Assume it's already in num, den format
            num, den = impedance_func
        
        # Check if suitable for Cauer synthesis
        if not self._is_suitable_for_cauer(num, den):
            print("Function not suitable for Cauer synthesis")
            return [], []
        
        kn_values = []
        num_current = num.copy()
        den_current = den.copy()
        
        # Initial step
        den_s_current = np.concatenate([[0], den_current])  # Multiply by s
        
        if len(num_current) == len(den_s_current):
            # Polynomial long division
            k = num_current[0] / den_s_current[0] if den_s_current[0] != 0 else 0
            if k < 0:
                print("Negative value encountered, cannot continue Cauer synthesis")
                return [], []
            kn_values.append(k)
            
            # Update for next iteration
            num_current = den_current.copy()
            den_current = num_current - k * den_s_current[1:]  # Remove leading zero
            den_current = self._remove_leading_zeros(den_current)
        
        # Continue Cauer synthesis
        iteration = 1
        while len(den_current) > 1:
            iteration += 1
            
            if len(num_current) == len(den_current) and iteration % 2 == 0:
                # Polynomial division step
                if den_current[0] != 0:
                    q = num_current[0] / den_current[0]
                    if q < 0:
                        print("Negative value encountered, cannot continue")
                        return [], []
                    kn_values.append(q)
                    
                    remainder = num_current - q * den_current
                    num_current = den_current.copy()
                    den_current = self._remove_leading_zeros(remainder)
                
            elif len(num_current) == len(den_current) + 1 and iteration % 2 == 1:
                # Cauer step
                den_s_current = np.concatenate([[0], den_current])
                
                if len(den_s_current) > 0 and den_s_current[0] != 0:
                    k = num_current[0] / den_s_current[0]
                    if k < 0:
                        print("Negative value encountered, cannot continue")
                        return [], []
                    kn_values.append(k)
                    
                    num_current = den_current.copy()
                    den_current = num_current - k * den_s_current[1:]
                    den_current = self._remove_leading_zeros(den_current)
            else:
                break
        
        # Final steps
        if len(den_current) == 1 and len(num_current) >= 1:
            if den_current[0] != 0:
                kn_values.append(num_current[0] / den_current[0])
            if len(num_current) > 1:
                kn_values.append(num_current[1] / den_current[0])
        
        # Separate into L and C values
        L_values = []
        C_values = []
        
        for i, k in enumerate(kn_values):
            if k < 0:
                print(f"Warning: Negative value at position {i}")
                continue
                
            if i % 2 == 0:  # Even index -> Inductor
                L_values.append(k)
            else:  # Odd index -> Capacitor
                C_values.append(k)
        
        return L_values, C_values
    
    def _is_suitable_for_cauer(self, num, den):
        """Check if the function is suitable for Cauer synthesis."""
        # Check if all coefficients are real
        if not np.all(np.isreal(num)) or not np.all(np.isreal(den)):
            return False
        
        # Additional checks can be added here
        return True
    
    def _remove_leading_zeros(self, poly):
        """Remove leading zeros from polynomial coefficients."""
        if len(poly) == 0:
            return poly
        
        # Find first non-zero coefficient
        start_idx = 0
        for i, coeff in enumerate(poly):
            if abs(coeff) > 1e-10:
                start_idx = i
                break
        
        return poly[start_idx:] if start_idx < len(poly) else [0]