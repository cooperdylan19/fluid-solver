# Using Riemann Solvers to emulate fluid behaviour

**Aim**: Produce shock tubes in terms of density, velocity, pressure, and specific internal energy. The Euler equations (governing fluid dynamics) were numerically solved using the Lax–Friedrichs scheme on both Cartesian and spherical grids. An attempt was made to implement a more accurate Riemann solver, the HLL scheme, although this was not successful.

### Shock tubes (LF scheme)


![Figure 1](/LFScheme_Plot1_CC.png)
<br />
**Figure 1**: Comparison of the exact shock tube (black line) to approximations determined using the Lax–Friedrichs scheme. The scheme was run once for 100 cells (light green) and again for 1000 cells (green). The initial conditions were: Left state (x < 0.3 at t = 0), ρ = 1, p = 1, v = 0.75. Right state: ρ = 0.125, p = 0.1, v = 0. Snapshots are taken at t = 0.2.
<br />

![Figure 2](/LFScheme_Plot2_CC.png)
<br />
**Figure 2**: Comparison of the exact shock tube (black line) to approximations determined using the Lax–Friedrichs scheme. The scheme was run once for 100 cells (light green) and again for 1000 cells (green). The initial conditions were: Left state (x < 0.5 at t = 0), ρ = 1, p = 0.4, v = −2. Right state: ρ = 1, p = 0.4, v = 2. Snapshots are taken at t = 0.15.
<br />

![Figure 3](/LFScheme_Plot1_SPC.png)
<br />
**Figure 3**: Comparison of the exact shock tube (black line) to approximations determined using the Lax–Friedrichs scheme. The scheme was run once for 100 cells (light green) and again for 1000 cells (green). Here, a spherical grid was used instead of a Cartesian grid. The initial conditions were: Left state (r < 0.4 at t = 0), ρ = 1, p = 1, v = 0. Right state: ρ = 0.125, p = 0.1, v = 0. Snapshots are taken at t = 0.25.
<br />

The attempt to implement the HLL scheme was unsuccessful (plots not shown). As the framework of the function was the same as that of the Lax–Friedrichs scheme, the error likely arose from a misunderstanding of how to compute the flux. Had the implementation been successful, we would expect the HLL scheme to provide a more accurate representation of the shock, as it is a more sophisticated Riemann solver. The plots show that increasing the number of cells yields results closer to the true shock profile.




