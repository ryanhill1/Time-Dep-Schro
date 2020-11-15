### Solving the time dependent Schrödinger equation using the Crank-Nicolson method

The time dependent Schrodinger equation for the wavefunction $\psi$ of a particle of mass $m$ moving in a potential energy $V(x,t)$ is:

$ i\hbar\frac{\partial}{\partial t}\psi = \frac{-\hbar^2}{2m}\frac{\partial^2}{\partial x^2} \psi(x,t) + V(x,t)\psi(x,t)$

In this program calculated the time dependent propogation of an electron wavepacket through a potential barrier. I performed the calculation in a region of $L=500$ Angstroms. I started with an initial (complex valued) Gaussian wave function (for an electron) of:

$$ \psi(x,t=0) = exp \bigg[-\bigg(\frac{x-0.3L}{s}\bigg)^2 + ixk_0 \bigg] $$

with a width of $s=10$ Angstroms and average wavenumber $k_0=1$ Angstroms$^{-1}$. The potential energy $V(x) models a one dimensional crystal surface with periodic peaks to mimic atomic layers in the crystal:

$$ V(x) = V_1 \bigg[ 0.75 - \cos \bigg(\frac{x-x_0}{\omega_x}\bigg)\bigg] \quad \text{for} \quad x > x_0$$

$$ = 0 \quad \text{otherwise}$$

where $V_1=2.0$ eV, $x_0=0.5L$, and $\omega_x=5$ Angstroms.

I used the Crank-Nicolson method for my calcuations. The Crank-Nicolson method is a finite difference method that can be used to solve partial differential equations. It is stable and accurate to $\mathcal{O}(\Delta t^2)$ globally.

A finite difference form of the Schrodinger equation for use in the Crank-Nicolson method is:

$$\psi(x-\Delta x, t+ \Delta t) + \bigg[ \frac{2 m \omega i}{\hbar} - 2 - \frac{2 m \Delta x^2}{\hbar^2} V(x) \bigg] \psi(x,t+\Delta t) + \psi(x+\Delta x, t + \Delta t) $$

where $\omega = 2\Delta x^2/\Delta t$, $\Delta x$ is the sampling size in space and $\Delta t$ is the sampling size in time. For this problem the Crank-Nicolson equation has the form:

$$a_j \psi(x_{j-1}, t_{n+1}) + b_j\psi(x_j,t_{n+1})+c_j(\psi(x_{j+1},t_{n+1})=d_j $$

which can be written as a tri-diagonal matrix equation. As shown in the equation above, the Crank-Nicolson involves solving a set of simultaneous equations. $\psi_i = \psi(x_i, t_{n+1})$ are the unknowns to be found. The two endpoints $\psi_{-1}$ and $\psi_{N_x}$ are fixed at 0.
