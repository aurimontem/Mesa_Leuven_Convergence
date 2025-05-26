
# MESA@Leuven Friday Morning Lab

# Best Practices: Convergence Testing 

In this brief morning lab session, we will go over best practices for solving partial differential equations numerically, specifically in the context of the MESA Stellar Evolution Code.

The lecture slides can be found at [this repository, to be updated soon](https://broken-url.com). 

Solutions can be found at [this repository, to be updated soon](https://broken-url.com). 

## Lab Overview: 
When solving a (partial) differential equation numerically, you are essentially approximating a _derivative_ as a _difference_. 
For a quantity $y$ which varies with some coordinate $x$, this resembles the following:

$$ \frac{\partial y}{\partial x} \approx \frac{\Delta y}{\Delta x} = \frac{y[i+1] - y[i]}{x[i+1] - x[i]}$$

This is a "forward-difference," since it evaluates the difference between y at a step forward, $i+1$, and $y$ at a given step, $i$, and using that to approximate the derivative at step $i$. In practice, modern numerical techniques are just _slightly_ fancier ways of approximating a derivative as a difference between zones or between times, and then solving the set of equations corresponding to the values at each point. 

How do we know this is a good approximation? Without diving deep into the formal theory of numerical errors, let's note that a derivative is defined by taking the limit of the slope, $\left(f(x + h) - f(x)\right)/h$, as $h$ approaches 0: 

$$
f'(x) =\lim _{h \rightarrow 0} \frac{f(a+h)-f(a)}{(a+h) - (a)}
$$

So, the essential question is: "Are we in the limit of small $h$ ?"

When we solve the same equations on a finer and finer grid, and find that the answers do not change for the quantities we care about, we call this "converged" or "numerically converged." 

We now turn to exploring this in the context of the MESA stellar evolution code. This will be broken up in 3 mini-mini-labs: 

In **Mini-mini lab 1**, we will explore changing resolution in space. In **Mini-mini lab 2**, we will explore changing resolution in time. In **Mini-mini lab 3**, we will explore changing physical approximations within reasonable model uncertainties. While this third task is not exactly the same as testing the "limit of small $h$", it is still **testing the numerical assumptions we are making in modeling the star as a sphere with finite shells**, so it is still an important and relevant aspect of convergence testing for astrophysical simulations. 

# Mini-mini lab 1: Spatial Resolution 

## Overview
In MESA, the fundamental spatial coordinate is the "mesh", which is broken up into "zones" (sometimes referred to as "cells" or "shells" or "mesh points") of varying mass $dm$, such that $\Sigma_i(dm_i)=M_* - m_\mathrm{IB}$ where $M_*$ is the star mass and $m_\mathrm{IB}$ is the mass inside the model inner boundary, which is 0 for most uses of MESA. The indexing is such that zone `1` corresponds to the surface of the star, and zone `nz` corresponds to the center of the star (or inner boundary). 

To help enforce that the zones are small enough that we are in fact "in the limit of small $h$", at each timestep, MESA can "adaptively" split and merge zones in order to achieve some tolerances in how various quantities vary from zone to zone. 

However, in choosing a mesh, MESA is guessing at what consitutes "small $h$". 
We can make it make better guesses, and we must always check it for errors. 

To change how MESA discretizes its mesh, we can do 3 things: 

1) We can tell MESA to increase or decrease the number of zones, e.g. take whatever it thinks the mesh should be and (double? triple? 10x? halve?) the number of zones. This is controlled by setting `mesh_delta_coeff` (`=1` by default). A smaller value means more grid points, with less delta (difference) between them. A larger value means fewer grid points, with larger allowed "delta" between them. 
   
2) We can also tell MESA to increase or decrease the tolerance for various physical targets directly. For example, perhaps MESA wants to have at most a relative change of 50% in density from zone `i` to zone `i+1`, and perhaps we think that's not good enough; we can specify that we want only 10% variations (Though, note that in this specific example you may end up with a TON of mesh points, because the density varies by tens of orders of magnitude between the core and the surface). There are _many_ controls for this; see `$MESA_DIR/star/defaults/controls.defaults` under the header
   ```fortran
   ! mesh adjustment   
   ! ===============
   ```
   Or the equivalent [on the MESA documentation website](https://docs.mesastar.org/en/latest/reference/controls.html#mesh-adjustment)

3) We can create our own custom mesh scheme in `src/run_star_extras.f90`. We may turn to this as a bonus task, time permitting.


## Mini-mini Lab 1 Instructions: 

Let's set up an example that illustrates (1) the importance of testing spatial resolution and (2) how _bad_ the spatial resolution of a lot of default MESA setups are. In general, we cannot emphasize enough that these labs, the `test_suite`, and the basic `$MESA_DIR/star/work` directory are NOT converged numerically. 

To see this, copy a clean work directory and enter:

```bash
cp -r $MESA_DIR/star/work ./work_space
cd work_space
```

The default work directory takes a $15M_\odot$ star and evolves it until ZAMS. In the massive star community, a lot of attention recently has been given to binarity, stellar winds, and other physics which changes the properties of the H-rich envelope (much of which we've discussed this week). So, for our purposes, let's engineer a partially-stripped star by cranking up the stellar wind, and evolve it until a slightly more advanced phase, namely, core He depletion. 

To do this, add the following to the `&controls` section of `inlist_project` under the header `! winds`: 

```fortran
hot_wind_scheme = 'Dutch'
cool_wind_RGB_scheme = 'Dutch'
cool_wind_AGB_scheme = 'Dutch'
Dutch_scaling_factor = dutchwindscale
```


Each member of your table should pick a different `mesh_delta_coeff`. For those with slower computers, you should choose larger values of `mesh_delta_coeff`. 


