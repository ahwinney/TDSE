# TDSE
TDSE

This is a numerical solution for the Time Dependent Schroedinger Equation (TDSE).
It propagates a wave function inside of a cylindrical geometry under the influence of an electric field characterized by a Gaussian wave form.
Under inense electromagnetic field, a phenomenon called High Harmonic Generation (HHG) takes place which generates odd (2n+1) harmonics of the source photon frequency.
Different boundary conditions and distances were tested to investigate the effects of limiting the space for HHG to take place.

Physical Results    
    We indeed observed reflections from the boundaries and as a result we observed increased HHG and higher maximum photon energy.
    We observed a similar effect when implementing absorbing boundary conditions. The absorbing boundary behaved like a very distant boundary, dampening reflection of the wave function and therefore dampening HHG.
    This result is significant because if it is possible to restrict the radial distance that the wave function can travel and reflect it, then HHG can be amplified including, at the highest energy harmonics.
    Increased efficiency of HHG means the posibiltiy for higher energy attosecond pulse generation and improved efficiency for electron dynamics studies.

Computational Results  
    Originally the code was written in MATLAB which can handle some linear algebra effiently, but for the scale of calculation we wanted, it was ineffective.
    After translating the code to C++, we saw a small simulation grid improve its calculation speed by 3x. This was originally a 24 hour simulation, and after translating it took about 8 hours.
    Further optimization decreased simulation time. But implementing the CUDA library with Arrayfire made the largest impact.
    Running the calculation on the GPU reduced the same simulation time to 30 min, an improvement of 16x from C++ and 48x from MATLAB.
    This increased efficiency allowed us to calculate much larger simulation sizes, thus approaching a result with practical implications.

The simulated electric field generated by an 800nm laser pulse.
A similar pulse was generated for other wavelengths. All are 6 periods wtih a gaussian distribution.
![Laser Pulse](https://github.com/ahwinney/TDSE/blob/f20af9e2a83011681069f7784f4477846ff16afc/Laser%20Pulse.JPG)

Heat mappped image of the wave function after interacting with the 800nm laser pulse.
Small difraction patterns form at the boudary.
![800nm](https://github.com/ahwinney/TDSE/blob/6a5b6c7fd89c74da7f8022e946b58c4db142acfb/800nm.JPG)

Same heat mapped image after interacting with a 2400nm laser pulse.
Much greater difraction patterns form at the boundary due to the longer wavelength and propagation times.
![Alt text](https://github.com/ahwinney/TDSE/blob/6a5b6c7fd89c74da7f8022e946b58c4db142acfb/2400nm.JPG)

The HHG yield for different boundary distances and a 800nm laser pulse.
HHG yield is seen to increases significantly with shorter boundary conditions.
![Alt text](https://github.com/ahwinney/TDSE/blob/6a5b6c7fd89c74da7f8022e946b58c4db142acfb/HHG.JPG)

The high energy cut-off harmonics at different boundary distances and laser wavelengths.
The yield increases by orders of magnitude as the boundary distance is reduced. This occurs at further boundaries for longer wavelengths.
![Alt text](https://github.com/ahwinney/TDSE/blob/6a5b6c7fd89c74da7f8022e946b58c4db142acfb/Cut-off%20Yield%20v%20Boundary.JPG)

A short video of the wave function for those who are interested. As the laser pulse propagates, the electric field drives the wave function back and forth, generating difraction patterns with itself.
[![Watch a short video of the wave function propagation](https://img.youtube.com/vi/0afxQRICnuQ/0.jpg)](https://youtu.be/0afxQRICnuQ)
