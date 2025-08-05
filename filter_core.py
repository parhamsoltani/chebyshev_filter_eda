"""
Core Chebyshev filter design functionality.
"""

import numpy as np
import sympy as sp
from sympy import symbols, sinh, cosh, asinh, sin, cos, pi, I, solve, limit, simplify
from scipy import signal

class ChebyshevCore:
    """Core class for Chebyshev filter design."""
    
    def __init__(self, N, e1, e2, Rg_denorm, RL_denorm):
        self.N = N
        self.e1 = e1
        self.e2 = e2
        self.e12 = e1 * e2
        self.Rg_denorm = Rg_denorm
        self.RL_denorm = RL_denorm
        self.Rg = Rg_denorm / RL_denorm
        self.RL = 1
        
        # Symbolic variables
        self.w = symbols('w', real=True)
        self.s = symbols('s')
        
        # Filter transfer functions
        self.TF_F1 = None
        self.TF_F2 = None
        self.NUM1 = None
        self.DEN1 = None
        self.NUM2 = None
        self.DEN2 = None
        
        # Poles
        self.pk11 = None
        self.pk12 = None
        self.pk2 = None
        
    def design_with_dc_drop(self):
        """Design filters with DC drop."""
        # Calculate poles for F11 and F12
        a11 = sinh(1/self.N * asinh(1/self.e1))
        b11 = cosh(1/self.N * asinh(1/self.e1))
        a12 = sinh(1/self.N * asinh(1/self.e2))
        b12 = cosh(1/self.N * asinh(1/self.e2))
        
        self.pk11 = np.zeros(self.N, dtype=complex)
        self.pk12 = np.zeros(self.N, dtype=complex)
        
        for i in range(self.N):
            angle = (2*i + 1) * pi / (2*self.N)
            self.pk11[i] = -a11*sin(angle) + I*b11*cos(angle)
            self.pk12[i] = -a12*sin(angle) + I*b12*cos(angle)
        
        # Calculate transfer functions F11 and F12
        NUM11 = 1/(self.e1 * 2**(self.N-1))
        NUM12 = 1/(self.e2 * 2**(self.N-1))
        
        DEN11 = self._poles_to_poly(self.pk11)
        DEN12 = self._poles_to_poly(self.pk12)
        
        # F1 = F11 * F12
        self.NUM1 = np.convolve(NUM11 if np.isscalar(NUM11) else [NUM11], 
                               NUM12 if np.isscalar(NUM12) else [NUM12])
        self.DEN1 = np.convolve(DEN11, DEN12)
        
        # Calculate F2
        a22 = sinh(1/(2*self.N) * asinh(1/self.e12))
        b22 = cosh(1/(2*self.N) * asinh(1/self.e12))
        
        self.pk2 = np.zeros(2*self.N, dtype=complex)
        for i in range(2*self.N):
            angle = (2*i + 1) * pi / (4*self.N)
            self.pk2[i] = -a22*sin(angle) + I*b22*cos(angle)
        
        self.NUM2 = 1/(self.e12 * 2**(2*self.N-1))
        self.DEN2 = self._poles_to_poly(self.pk2)
        
        # Create transfer functions
        self.TF_F1 = signal.TransferFunction(self._clean_coeffs(self.NUM1), 
                                           self._clean_coeffs(self.DEN1))
        self.TF_F2 = signal.TransferFunction(self._clean_coeffs(self.NUM2), 
                                           self._clean_coeffs(self.DEN2))
        
        return self.TF_F1, self.TF_F2
    
    def design_without_dc_drop(self):
        """Design filters without DC drop (for even N)."""
        if self.N % 2 != 0:
            raise ValueError("DC-drop-free design only available for even N")
        
        # Generate Chebyshev polynomials
        TNw_N = self._generate_chebyshev_poly(self.N)
        TNw_2N = self._generate_chebyshev_poly(2*self.N)
        
        # Calculate rotated Chebyshev polynomial
        rotation_arg = sp.sqrt((sp.cos(pi/(2*self.N))*self.w)**2 + (sp.sin(pi/(2*self.N)))**2)
        LNw_N = TNw_N.subs(self.w, rotation_arg)
        
        # Calculate poles for F11
        F11_squared_den = 1 + (self.e1*LNw_N)**2
        F11_squared_den_s = F11_squared_den.subs(self.w, self.s/I)
        poles_F11_squared = solve(F11_squared_den_s, self.s)
        
        # Select left half-plane poles
        poles_F11 = [p for p in poles_F11_squared if sp.re(p) < 0][:self.N]
        self.pk11 = np.array([complex(p.evalf()) for p in poles_F11])
        
        # Similar process for F12
        F12_squared_den = 1 + (self.e2*LNw_N)**2
        F12_squared_den_s = F12_squared_den.subs(self.w, self.s/I)
        poles_F12_squared = solve(F12_squared_den_s, self.s)
        poles_F12 = [p for p in poles_F12_squared if sp.re(p) < 0][:self.N]
        self.pk12 = np.array([complex(p.evalf()) for p in poles_F12])
        
        # Calculate numerator constants
        NUM11 = np.prod([0 - p for p in poles_F11])
        NUM12 = np.prod([0 - p for p in poles_F12])
        
        self.NUM1 = np.convolve([NUM11] if np.isscalar(NUM11) else NUM11,
                               [NUM12] if np.isscalar(NUM12) else NUM12)
        self.DEN1 = np.convolve(self._poles_to_poly(self.pk11), 
                               self._poles_to_poly(self.pk12))
        
        # F2 calculation
        LNw_2N = TNw_2N.subs(self.w, sp.sqrt((sp.cos(pi/(4*self.N))*self.w)**2 + 
                                             (sp.sin(pi/(4*self.N)))**2))
        F2_squared_den = 1 + (self.e12*LNw_2N)**2
        F2_squared_den_s = F2_squared_den.subs(self.w, self.s/I)
        poles_F2_squared = solve(F2_squared_den_s, self.s)
        poles_F2 = [p for p in poles_F2_squared if sp.re(p) < 0][:2*self.N]
        self.pk2 = np.array([complex(p.evalf()) for p in poles_F2])
        
        NUM2 = np.prod([0 - p for p in poles_F2])
        self.NUM2 = NUM2
        self.DEN2 = self._poles_to_poly(self.pk2)
        
        # Create transfer functions
        self.TF_F1 = signal.TransferFunction(self._clean_coeffs(self.NUM1), 
                                           self._clean_coeffs(self.DEN1))
        self.TF_F2 = signal.TransferFunction(self._clean_coeffs(self.NUM2), 
                                           self._clean_coeffs(self.DEN2))
        
        return self.TF_F1, self.TF_F2
    
    def _generate_chebyshev_poly(self, n):
        """Generate Chebyshev polynomial of order n."""
        if n == 0:
            return sp.Integer(1)
        elif n == 1:
            return self.w
        elif n == 2:
            return 2*self.w**2 - 1
        elif n == 3:
            return 4*self.w**3 - 3*self.w
        else:
            # Use recurrence relation: T_n(x) = 2*x*T_{n-1}(x) - T_{n-2}(x)
            T = [sp.Integer(1), self.w]
            
            # Calculate powers of 2 first
            k = 2
            while 2**k <= n:
                T.append(2*T[2**(k-1)]**2 - 1)
                k += 1
            
            # Fill remaining terms using recurrence
            for i in range(len(T), n + 1):
                if i > 3:
                    T.append(2*self.w*T[i-1] - T[i-2])
            
            return T[n]
    
    def _poles_to_poly(self, poles):
        """Convert poles to polynomial coefficients."""
        poly_coeffs = [1]
        for pole in poles:
            # Multiply by (s - pole)
            new_coeffs = [0] * (len(poly_coeffs) + 1)
            for i, coeff in enumerate(poly_coeffs):
                new_coeffs[i] += coeff
                new_coeffs[i+1] -= coeff * pole
            poly_coeffs = new_coeffs
        return np.array([complex(c).real if abs(complex(c).imag) < 1e-10 
                        else complex(c) for c in poly_coeffs])
    
    def _clean_coeffs(self, coeffs):
        """Clean up coefficients by removing small imaginary parts."""
        if np.isscalar(coeffs):
            return coeffs
        cleaned = []
        for c in coeffs:
            if isinstance(c, complex):
                if abs(c.imag) < 1e-10:
                    cleaned.append(c.real)
                else:
                    cleaned.append(c)
            else:
                cleaned.append(c)
        return np.array(cleaned)