
# MESA@Leuven Best Practices Lab — Convergence Testing

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

In **Mini-mini lab 1**, we will explore changing resolution in space and time. In **Mini-mini lab 2**, we will briefly discuss what to do when a resolution test _fails_. 
In **Mini-mini lab 3**, we will explore changing physical approximations within reasonable model uncertainties. Though not explicitly about numerical resolution testing, the third task is likewise **testing the numerical assumptions we are making in modeling the star as a sphere with finite shells**, so it is still a relevant aspect of convergence testing for astrophysical simulations. 

## MESA-specific background

### Lagrangian Mesh 
In the equations MESA solves, the fundamental spatial coordinate is made up of concentric shells each with a given mass. This is often referred to as the "mesh", which is broken up into "zones" (sometimes referred to as "cells" or "shells" or "mesh points"). The mass per zone $dm$ can vary, under the constraint that $\Sigma_i(dm_i)=M_* - m_\mathrm{IB}$ where $M_*$ is the star mass and $m_\mathrm{IB}$ is the mass inside the model inner boundary, which is 0 for most uses of MESA. The indexing is such that zone `1` corresponds to the surface of the star, and zone `nz` corresponds to the center of the star (or inner boundary). 

To help enforce that the zones are small enough that we are in fact "in the limit of small $h$", at each timestep, MESA can "adaptively" split and merge zones in order to achieve some tolerances in how various quantities vary from zone to zone. However, in choosing a mesh, MESA is guessing at what consitutes "small $h$". 

We can make MESA make better guesses, and we must always check it for errors. To change how MESA discretizes its mesh, we can do 3 things: 

S1) We can tell MESA to increase or decrease the number of zones, e.g. take whatever it thinks the mesh should be and (double? triple? 10x? halve?) the number of zones. This is controlled by setting `mesh_delta_coeff` (`=1` by default). A smaller value means more grid points, with less delta (difference) between them. A larger value means fewer grid points, with larger allowed "delta" between them. 
   
S2) We can also tell MESA to increase or decrease the tolerance for various physical targets directly. For example, perhaps MESA wants to have at most a relative change of 50% in density from zone `i` to zone `i+1`, and perhaps we think that's not good enough; we can specify that we want only 10% variations (Though, note that in this specific example you may end up with a TON of mesh points, because the density varies by tens of orders of magnitude between the core and the surface). There are _many_ controls for this; see `$MESA_DIR/star/defaults/controls.defaults` under the header
   ```fortran
   ! mesh adjustment   
   ! ===============
   ```
   Or the equivalent [on the MESA documentation website](https://docs.mesastar.org/en/latest/reference/controls.html#mesh-adjustment)

S3) We can create our own custom mesh scheme in `src/run_star_extras.f90`. We may turn to this as a bonus task, time permitting.

### Adaptive Timesteps

MESA is an implicit code, meaning it chooses its timestep adaptively and iterates until it achieves a solution within specified tolerances (i.e. specified differences between the right-hand-side and left-hand-side of the equations it's solving, as other controls on how much one model can deviate from one timestep to another). If the errors are too large in a given timestep, then MESA will cut the timestep in an attempt to get closer to "the limit of small $h$" (where now $h$ represents an increment in time $dt$). 

However, like choosing a mesh, MESA is guessing at what consitutes "small $h$". To change how MESA selects its timestep, we can likewise have a few options: 

T1) We can tell MESA to multiply the timestep it originally selects by a `time_delta_coeff`, analogous to `mesh_delta_coeff` (`=1` by default). Like with mesh, a smaller value means finer sampling in time, with less delta (difference) between them. A larger value means coarser sampling in time, with larger allowed "delta" between them. 
   
T2) We can tell MESA to increase or decrease the tolerance for various physical targets directly, in the `&controls` section of the inlist. Note that at different evolutionary phases, different tolerances may be setting the timestep. 
```fortran
   ! timestep controls
   ! =================
```
A broad control that's often used is `varcontrol_target` which specifies how much the model should deviate (defined with a broad metric encompassing a handful of physical quantities) from timestep to timestep. 
A summary of these may also be found [on the MESA documentation website](https://docs.mesastar.org/en/latest/reference/controls.html#timestep-controls)

T3) We can also tell MESA to cut or increase the timestep in `src/run_star_extras.f90` via user-defined criteria, especially in the `extras_check_model` routine, by directly manipulating the `s% dt` in the star info structure. 

# Mini-mini Lab 1: Spatial and Temporal Resolution testing

Let's set up an example that illustrates (1) the importance of testing resolution and (2) how _bad_ the default resolution in MESA is for certain regimes. 
In general, we cannot emphasize enough that these labs, the `test_suite`, and the basic `$MESA_DIR/star/work` directory are NOT converged numerically. 

DO NOT USE DEFAULTS AS A STARTING POINT FOR SCIENCE RUNS UNLESS YOU HAVE DONE ROBUST RESOLUTION TESTING! 

Let's start as close to MESA defaults as possible. Copy a clean work directory and enter it:

```bash
cp -r $MESA_DIR/star/work ./work_res1
cd work_res1
```

The default work directory takes a $15M_\odot$ star and evolves it until ZAMS. In the massive star community, a lot of attention recently has been given to stellar winds, binarity, and other physics which may impact the properties of the H-rich envelope (much of which we've discussed this week). Let's evolve that model until it's later along its core He burning phase. 

First, so that everyone can easily share their last HR diagram with eachother at the table, add the following to the `&star_job` section of `inlist_project`: 

```fortran
pause_before_terminate = .true.
```

Let's also start on the main sequence, to save us a minute or two of runtime. Add the following to the `&star_job` section of `inlist_project`: 
```fortran
create_pre_main_sequence_model = .false. ! previously .true. 
```

To stop during core He burning, change the following in the  `&controls` section of `inlist_project`: 
```fortran
stop_near_zams = .false. ! previously .true.
```
and 
```fortran
xa_central_lower_limit_species(1) = 'he4' ! previously 'h1'
xa_central_lower_limit(1) = 2d-1          ! previously 1d-3
```

It is often good practice to change your time and mesh resoltion together, though in principle these can be varied independently. Today we will use methods S1/T1. We mention the other methods above in order to remind the user that there are lots of good ways to do things in MESA depending on your problem. 

Each member of your table should pick a unique `mesh_delta_coeff` from the set [0.3, 0.5, 1, 2]. Make sure everyone at your table is exploring a different value. 

**NOTE: For those with slower computers, you should choose larger values of `mesh_delta_coeff`. ** If you have a very fast computer, feel free to try other values, but it's recommended not to go below 0.2 for the sake of time in this lab block (or, if you do, be prepared to kill the run). 

To change your resolution, add the following controls to your inlist, replacing `VALUE` with the appropriate value: 

```fortran
! timesteps
time_delta_coeff = VALUE ! 1 by default
max_model_number = 2000 ! off by default -- putting a cap here in case things get too crazy

! mesh
mesh_delta_coeff = VALUE ! 1 by default
max_allowed_nz = 16000 ! default 8000
```

Great, you're ready to run! In the terminal, from your working directory, clean make and run! 

```bash
./clean && ./mk && ./rn 
```

Watch the run evolve, and watch the runs of others at your table. Compare the HR diagram that pops up with those produced by people at your table with a different mesh_delta_coeff / time_delta_coeff. Do your diagrams agree? Disagree? Which agree better? 
For comparison, record the final **Mass**, **Radus**, **$T_\mathrm{eff}$**, **Luminosity**, and **star age**. 




# Mini-mini lab 2: Temporal Resolution 

## Overview
MESA takes timesteps that it chooses based on various criteria. To help enforce that the time steps are small enough that we are in fact "in the limit of small $h$" where now $h$ represents an increment in time $dt$. 

Like choosing a mesh, MESA is guessing at what consitutes "small $h$" (though MESA is doing this adaptively, picking a timestep to reach tolerance targets). GIVE MORE EXPLANATION OF THIS. 

To change how MESA selects its timestep, we can likewise do 3 things: 

1) We can tell MESA to multiply the timestep it originally selects by a `time_delta_coeff`, analogous to `mesh_delta_coeff` (`=1` by default). Like with mesh, a smaller value means finer sampling in time, with less delta (difference) between them. A larger value means coarser sampling in time, with larger allowed "delta" between them. 
   
2) We can also tell MESA to increase or decrease the tolerance for various physical targets directly.
  
FIX BELOW: 
5) For example, perhaps MESA wants to have at most a relative change of 50% in density from zone `i` to zone `i+1`, and perhaps we think that's not good enough; we can specify that we want only 10% variations (Though, note that in this specific example you may end up with a TON of mesh points, because the density varies by tens of orders of magnitude between the core and the surface). There are _many_ controls for this; see `$MESA_DIR/star/defaults/controls.defaults` under the header
   ```fortran
   ! mesh adjustment   
   ! ===============
   ```
   Or the equivalent [on the MESA documentation website](https://docs.mesastar.org/en/latest/reference/controls.html#mesh-adjustment)

6) We can create our own custom timestepping in `src/run_star_extras.f90`, but in this case instead of a `use_other_mesh` routine, we will choose timesteps in `extras_finish_step` and `extras_check_model`. We may turn to this as a bonus task, time permitting.

## Mini-mini Lab 2 Instructions: 

Let's set up an example that illustrates (1) the importance of testing spatial resolution and (2) how _bad_ the spatial resolution of a lot of default MESA setups are. In general, we cannot emphasize enough that these labs, the `test_suite`, and the basic `$MESA_DIR/star/work` directory are NOT converged numerically. 

To see this, copy a clean work directory and enter:

```bash
cp -r $MESA_DIR/star/work ./work_time
cd work_time
```

Let's do the same thing as before, adding the following to the `&star_job` section of `inlist_project`: 

```fortran
pause_before_terminate = .true.
create_pre_main_sequence_model = .false. ! previously .true. 
```

Just like before, change the winds and stopping condition in the `&controls` section of `inlist_project`: 

```fortran
hot_wind_scheme = 'Dutch'
cool_wind_RGB_scheme = 'Dutch'
cool_wind_AGB_scheme = 'Dutch'
Dutch_scaling_factor = 1

stop_near_zams = .false. ! previously .true.
xa_central_lower_limit_species(1) = 'he4' ! previously 'h1'
xa_central_lower_limit(1) = 2d-1          ! previously 1d-3
```

Next, each member of your table should pick a different `time_delta_coeff` from the set [0.3, 0.5, 1, 2].
For those with slower computers, you should choose larger values of `time_delta_coeff`. 
If you have a very fast computer, you can try other values, but it's recommended not to go below 0.2 for the sake of this exercise 
(or, if you do, be prepared to kill the run). 

Set this by adding the following to your inlist and replacing `VALUE` with the appropriate value: 
```fortran
! timesteps
time_delta_coeff = VALUE 
```

Now you're ready! In the terminal, from your working directory, clean make and run! 

```bash
./clean && ./mk && ./rn 
```

Compare the HR diagram that pops up with those produced by people at your table with different timestepping. Do your diagrams agree? Disagree? Which agree better? 
For comparison, record the final **Mass**, **Radus**, **$T_\mathrm{eff}$**, **Luminosity**, and **star age**. 



One more note: Often it is good practice to run the resolution test by varying the spatial and temporal resolution together. 
