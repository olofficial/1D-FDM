# Heat conduction in 1D
## Physical system
Assumption: A one-dimensional temperature field $T(x, t)$ dependent on the position $x$ and time $t$.
The heat equation in 1D is 
$$
\begin{align}
\frac{ \partial T }{ \partial t } = \lambda \frac{ \partial^{2} T }{ \partial x^{2} },
\end{align}$$
where $\lambda$ is the thermal conductivity of the material.

The system is defined by a 1D interval of length $L$, and time starting at $t = t_0$.

## Discrete approximation
Instead of continuous coordinates, the interval $A: x \in [0, L]$ and $B: t \in [t_0, t_1]$ are divided into $N_1$ and $N_2$ steps, respectively:
$$
\begin{align}
&\Delta x = \frac{L}{N_{1}}\\
&\Delta t = \frac{t_{1}-t_{0}}{N_{2}}
\end{align}
$$

### First order derivative
The first order derivative is approximated by the Euler forward method:
