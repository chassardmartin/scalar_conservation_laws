
This repo is the synthesis of an internship work done at the CERMICS lab of Ecole des Ponts ParisTech with [V. Ehrlacher](https://team.inria.fr/matherials/team-members/virginie-ehrlacher-galland/) devoted to the development of codes for non-local (and non-linear) *scalar conservation laws*. 


## Scalar conservation laws

A (non-linear) scalar conservation law with flux $f: \mathbb{R} \rightarrow \mathbb{R}$ describes the evolution of a concentration $u: \mathbb{R} \times [0, +\infty) \rightarrow \mathbb{R}$ according to 

$$
    \partial_t u(x,t) + \partial_x f(u(x,t)) = 0
$$

The taxonomy "conservation law" is explained when ones conducts the following computation 

$$
    \frac{d}{dt} \left( \int_{(a,b)} u(x,t) \, dx\right) = \int_{(a,b)} \partial_t u(x,t) \, dx = - \int_{(a,b)} \partial_x f(u(x,t)) \, dx
    = f(u(a, t)) - f(u(b, t)),
$$

so that the total concentration in $(a,b)$ is evolving exactly with respect to the flux on the boundaries, without  "loss" and the total concentration is conserved.

Here are some examples of well-known scalar conservation laws :

- **Transport equation** : $f(u) = b\,u$, with $b$ constant, the velocity. (linear)
- **Fick's Law (or Heat equation)** : $f(u) = - D \,\partial_x u $, with $D>0$ the diffusion coefficient. (linear)
- **Burger's equation** : $f(u) = u^2/2$. (non-linear)

## Non-local scalar conservation laws

A non-local scalar conservation law is a conservation law that has a slightly different form than the one above 

$$ 
\partial_t u(x,t) + \partial_x f(V(u(\cdot,t) \ast \omega)(x)) = 0,
$$

where $\ast$ is the convolution product 

$$
    (u(\cdot,t) \ast \omega) (x) = \int_{\mathbb{R}} u(x-y, t) \, \omega(y) \, dy,
$$
,$\omega$ is an integrable kernel and $V$ is a velocity function. The convolution is precisely what makes the conservation law "non-local" as it requires the knowledge of the function $u$ on the whole space (more precisely on the support of a translation of the kernel $\omega$) whereas a derivation is computed point-wise, that is, locally. 

## Schemes 

We present in this repo two different schemes to solve the non-local scalar conservation law : 

- a semi-lagrangian scheme, using the characteristics to derive an analytic expression of the solution, alongside a Fourier interpolation, which suffers from the Gibbs phenomenon, but requires no $\text{CFL}$-type condition. The associated notebook is ```fourier_int.ipynb```.

- a "lagrangian-eulerian" scheme directly extracted from [E. Abreu et al. 2022] that lets the space discretization evolve with the characteristic curves for a better approximation of discontinuities, but which suffers from a restrictive $\text{CFL}$ condition. The associated notebook is ```lag_euler.ipynb```. We added a translation of the code in Julia for performance issues : for identical discretization parameters, a simulation of time $T=0.1$ s took 23min for the python code against 12s on the Julia version. The julia code is ```lag_euler.jl```.

## Numerical Examples

- We have used those schemes to study (under the label "burgers") a non-local extension of Burger's equation with velocity $V(x) = 1-x$ as a simple example of non-linearity, using a parametrized kernel with form 
$$
    \omega_{\eta}(x) = \mathbb{1}_{[0,\eta)}(x) (1/\eta),
$$
which implies that $\text{supp}(\omega) = [0, \eta]$, and a gaussian initial condition $u_0(x) = \frac{4}{\sqrt{7\pi}} \exp(-x^2/7)$.

- We also studied (under the label "rectangle") the evolution of a characteristic function as initial condition $u_0(x) = \mathbb{1}_{[0,1)}(x)$ and $V(x)=x$, with kernel 
$$
\omega(x) = \mathbb{1}_{[0, 1)}(x),
$$ 
and flux $f(u) = u$.

Plots can be seen in the graphs/ dir and source codes can be found in the code/ dir. 


# References 

1) E. Abreu, R. De la Cruz, JC. Juajibioy, W. Lambert, *Lagrangian-eulerian approach for nonlocal conservation laws* in Journal of Dynamics and Differential Equations, Springer, 2022.

