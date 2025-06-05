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

Practical Results  
    Originally the code was written in MATLAB which can handle some linear algebra effiently, but for the scale of calculation we wanted, it was ineffective.
    After translating the code to C++, we saw a small simulation grid improve its calculation speed by 3x. This was originally a 24 hour simulation, and after translating it took about 8 hours.
    Further optimization decreased simulation time. But implementing the CUDA library with Arrayfire made the largest impact.
    Running the calculation on the GPU reduced the same simulation time to 30 min, an improvement of 16x from C++ and 48x from MATLAB.
    This increased efficiency allowed us to calculate much larger simulation sizes, thus approaching a result with practical implications.

![Laser Pulse](https://github.com/ahwinney/TDSE/blob/f20af9e2a83011681069f7784f4477846ff16afc/Laser%20Pulse.JPG)]

![Alt text](image-url)

![Alt text](image-url)

![Alt text](image-url)

![Alt text](image-url)

[![Watch a short video of the wave function propagation](https://img.youtube.com/vi/0afxQRICnuQ/0.jpg)](https://youtu.be/0afxQRICnuQ)
