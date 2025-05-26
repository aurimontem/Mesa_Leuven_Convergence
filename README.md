
# MESA@Leuven Friday Morning Lab

# Best Practices: Convergence Testing 

In this brief morning lab session, we will go over best practices for solving partial differential equations numericall, with applications to the MESA Stellar Evolution Code.

The lecture slides can be found at [this repository, to be updated soon](https://broken-url.com). 

# Overview: 
When solving a (partial) differential equation numerically, you are essentially approximating a _derivative_ as a _difference_. 
For a quantity $y$ which varies with some coordinate $x$, this resembles the following:

$$ \frac{\partial y}{\partial x} \approx \frac{\Delta y}{\Delta x} = \frac{y[i+1] - y[i]}{x[i+1] - x[i]}$$

This is a "forward-difference," since it is the difference between y at a step forward, $i+1$, and $y$ at a given step, $i$, and using that to approximate the derivative at step $i$. In practice, modern numerical techniques are just _slightly_ fancier ways of approximating a derivative as a difference between zones or between times, and then solving the set of equations corresponding to the values at each point. 

How do we know this is a good approximation? Without diving deep into the formal theory of numerical errors, let's note that a derivative 

$$
f'(x) =\lim _{h \rightarrow 0} \frac{f(a+h)-f(a)}{(a+h) - (a)}
$$

So, the essential question is: "Are we in the limit of small $h$ ?"

We now turn to exploring this in the context of the MESA stellar evolution code. 

# Mini-mini lab 1: Spatial Resolution

In MESA, the fundamental spatial coordinate is the "mesh", which is broken up into "zones" (sometimes referred to as "cells" or "shells") of varying mass $dm$, such that $\Sigma_i(dm_i)=M_* - m_\mathrm{IB}$ where $M_*$ is the star mass and $m_\mathrm{IB}$ is the mass inside the model inner boundary, which is 0 for most uses of MESA. 

To help ensure that the zones are small enough such that we are in fact "in the limit of small $h$", at each timestep, MESA can "adaptively" split and merge zones in order to achieve some tolerances in how various quantities vary from zone to zone. 


