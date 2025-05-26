
# MESA@Leuven Friday Morning Lab

# Best Practices: Convergence Testing 

In this brief morning lab session, we will go over best practices for solving partial differential equations numericall, with applications to the MESA Stellar Evolution Code.

The lecture slides can be found at [this repository, to be updated soon](https://broken-url.com). 

# Overview: 
When solving a (partial) differential equation numerically, you are essentially approximating a _derivative_ as a _difference_. For a quantity $y$ which varies with some coordinate $x$, this looks like

$$ \frac{\partial y}{\partial x} \approx \frac{\Delta y}{\Delta x} = \frac{y[i+1] - y[i]}{x[i+1] - x[i]}$$

Shown above is a ``forward-difference," since you are taking the difference between y at step $i+1$ and $y$ at step $i$, and using that to approximate the derivative at step $i$; in practice, modern numerical techniques have _slightly_ fancier ways of doing this. 

# Mini-mini lab 1: 
