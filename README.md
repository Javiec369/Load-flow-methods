# Numerical Methods and Power Systems Analysis in Julia Repository

> Universidad Tecnológica de Pereira\
> Electrical enginnering program 

This repository is a collection of exercises and projects related to the programming of numerical methods and power systems analysis in the programming language Julia.

This repository uses knowledge acquired in the learning of different numerical methods, such as:

> PLease keep in mind that excel file mustb be dowloaded in oder to analize the load flow

 - **End Point Method**: Iterative calculation method that allows finding the end point of a dynamic system.
	 $$x_{k+1} = f(x_{k},t_{k},\Delta t) $$
 - **Fixed Point Method**: Numerical method for finding the fixed point of a function.
	 $$x = f(x)$$
 - **Ybus Formulation**: Method for calculating the Ybus admittance matrix of an electrical system.
	$$Y_{bus} = \sum_{i=1}^{n} Y_{ij}$$
 - **Load Flow - Quasy Dynamic Flow**: Power system analysis method for calculating power flows $[2]$.
	$$P_{i} - P_{i}^{D} - \sum_{j=1}^{n} V_{i}V_{j}(G_{ij}cos\theta_{ij} + B_{ij}sin\theta_{ij}) = 0$$
	$$Q_{i} - Q_{i}^{D} - \sum_{j=1}^{n} V_{i}V_{j}(G_{ij}sin\theta_{ij} - B_{ij}cos\theta_{ij}) = 0$$
 - **Newton's Method**: Method for solving non-linear equations.
	$$f(x_{k}) + J(x_{k})(x_{k+1}-x_{k}) = 0$$

 - **Linearized Load Flow**: A linear approximation is developed on the complex numbers $[1]$ and not on the reals as in the conventional load flow formulations.
$$\frac{1}{1-\Delta V}=\sum_{n=0}^{+\infty}(\Delta V)^{n},\hspace{10pt}|\Delta V|<1$$
 
Where:    
    - $x_k$ is the current value of the variable in the iterative method.\
    - $f$ is the function describing the dynamical system.\
    - $t_k$ is the current time in the dynamic system.\
    - $\Delta t$ is the time increment in the dynamical system.\
    - $x$ is the fixed point of the function $f(x)$.\
    - $Y_{bus}$ is the admittance matrix of the electrical system.\
    - $Y_{ij}$ is the admittance between node $i$ and node $j$.\
    - $P_{i}$ is the active power at node $i$.\
    - $P_{i}^{D}$ is the active power demanded at node $i$.\
    - $Q_{i}$ is the reactive power at node $i$.\
    - $Q_{i}^{D}$ is the reactive power demanded at node $i$.
    - $V_{i}$ is the voltage at node $i$.\
    - $G_{ij}$ is the conductance between node $i$ and node $j$.\
    - $B_{ij}$ is the susceptance between node $i$ and node $j$.\
    - $\theta_{ij}$ is the phase angle between node $i$ and node $j$.\
    - $f(x_k)$ is the function evaluated at $x_k$.\
    - $J(x_{k})$ is the Jacobian matrix evaluated on $x_k$.\
    - $x_{k+1}$ is the next value of the variable in the iterative method.\
    - $\Delta P$ is the variation of the active power in the power system.\
    - $\Delta Q$ is the variation of the reactive power in the electrical system.

## Power Systems Analysis

In this repository, the analysis of power systems is performed using different numerical methods. Knowledge in the field of electrical engineering and, in particular, in the programming of load flow methods is used.

## Julia Programming Language

Julia is a high-performance programming language designed to perform numerical and scientific calculations. Julia is especially useful for programming numerical and parallel algorithms, as well as for data analysis and visualisation.

For more information about Julia, please visit its [homepage](https://julialang.org/).

## Contributions

If you wish to contribute to this repository, please send your pull requests with new implementations of numerical methods and power system analysis in Julia. You can also send proposals for improvements or bug fixes  all contributions are welcome!


## References

*[1]: A. Garces, "A Linear Three-Phase Load Flow for Power Distribution Systems," IEEE Transactions on Power Systems, vol. 31, no. 1, pp. 480-487, Jan. 2016.*

*[2]: A. Garces-Ruiz, "Power Flow in Unbalanced Three-Phase Power Distribution Networks Using Matlab: Theory, analysis, and quasi-dynamic simulation," Revista Ingeniería, received: 2nd-April-2022, modified: 28th-April-2022, accepted: 23th-May-2022.*