# Forced Vibration Analysis via Fourier Series
## 📄 Project Overview

In this project, the motion of a single-degree-of-freedom system under a piecewise-polynomial forcing function is examined:

1. The Fourier series of the forcing function \(f(t)\) was derived by hand, splitting it into two intervals:
   - \(0 < t < 1\): \(f_1(t) = 100\,(a_1 t^5 + a_2 t^4 + a_3 t^3 + a_4 t^2 + a_5 t)\)  
   - \(1 < t < 1.5\): \(f_2(t) = -500\)

2. A general-purpose Python module (`fourier_calculator.py`) was implemented using Sympy to:
   - Compute the Fourier coefficients \(a_0, a_n, b_n\).  
   - Assemble and evaluate the truncated Fourier series for any chosen \(k\).  
   - Generate plots of \(f(t)\) for various \(k\), compute amplitude–frequency spectra, and visualize the response.

3. The trade-off between accuracy and computation time was analyzed, resulting in a truncation order of **k = 20**.

4. The system response was obtained by:
   - Computing the small-angle displacement \(x(t)\).  
   - Determining the support reaction force \(F_A(t)\) via static equilibrium.  


