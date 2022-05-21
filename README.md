# Spectral-Method-Capillarity
Spectral methods for capillary surfaces

GPL-3.0-or-later

As of May 5th, 2022, three prototype problems have been treated, with corresponding codes.

Each of these programs needs Chebfun.  Please install Chebfun from chebfun.org.

# P1
symmetric_capillary_disc.m

This code solves for capillary surfaces that are the image of a disk.  The inputs are the radius $b$ and the inlination angle $\psi_b$ there.  We require $-\pi \leq\psi_b\leq\pi$.  The output is the surface.

# P2
symmetric_capillary_annular.m 

This code solves for capillary surfaces that are the image of an annulus.  The inputs are the radii $a$ and $b$ with $0<a<b<\infry$ and the inlination angles $\psi_a$ and $\psi_b$ there. We require $-\pi \leq\psi_a,\psi_b\leq\pi$.  The output is the surface.

# P3
symmetric_capillary_2D.m 

This code solves for lower simensional capillary surfaces that are the analogue of either of the above problems.  The inputs are the radii $a$ and $b$ with $-\infty<a<b<\infry$ and the inlination angles $\psi_a$ and $\psi_b$ there. We require $-\pi \leq\psi_a,\psi_b\leq\pi$.  The output is the surface.

# Comments

For modest problems, the convergence to solutions is extremely fast and uses very little memory.

For more challenging problems the adaptive algorithm automatically increases the number of Chebyshev points to achieve the prescribed tolerances.  This precess has worked very well in almost all of the cases we have found.  

It is possible sometimes to break these codes for extremes of the inclination angles near $\pm\pi$ and either radii too close to each other with $\psi_a\psi_b > 0$, inner radius $a$ too close to 0, $b$ too large, or the  radii too far apart .  Typically this happens when $|\psi_a|, |\psi_b| > \pi/2$, and the closer to $\pm\pi$, the more challenging.   For these multi-scale problems customizing the initial guess and carefully tuning the increase of the Chebyshev nodes inside the adaptive part of the algorithm can lead to success when the base code does not converge.  Still, some rather extreme problems may  not work even then, and for those problems a multi-scale approach is needed.  This is the subject of a future paper.
