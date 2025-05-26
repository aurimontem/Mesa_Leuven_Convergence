
# MESA@Leuven Friday Morning Lab

# Best Practices: Convergence Testing 

In this brief morning lab session, we will go over best practices for solving partial differential equations numericall, with applications to the MESA Stellar Evolution Code.

The lecture slides can be found at [this repository, to be updated soon](https://broken-url.com). 

# Overview: 
When solving a (partial) differential equation numerically, you are essentially approximating a _derivative_ as a _difference_. For a quantity $y$ which varies with some coordinate $x$, this looks like

$$ \frac{\partial y}{\partial x} \approx \frac{\Delta y}{\Delta x} = \frac{y[i+1] - y[i]}{x[i+1] - x[i]}$$

Shown above is a ``forward-difference," since you are taking the difference between y at step $i+1$ and $y$ at step $i$, and using that to approximate the derivative at step $i$; in practice, modern numerical techniques have _slightly_ fancier ways of doing this. 

How do we know this is a good approximation? Without diving deep into the formal theory of numerical errors, let's note that a derivative 

$$
f'(x) =\lim _{h \rightarrow 0} \frac{f(a+h)-f(a)}{(a+h) - (a)}
$$

So, the essential question is: "Are we in the limit of small $h$ ?"

# Mini-mini lab 1: Spatial Resolution

In MESA, 
