# Lsode

This package enables you to solve systems of ODEs (Ordinary Differential Equations) by calling the DLSODE function from the ODEPACK package (FORTRAN).
To solve differential equations in Julia I highly recommend using the [DifferentialEquations](https://github.com/JuliaDiffEq/DifferentialEquations.jl) package.
I created this package mainly to try out Julia's ability to call FORTRAN code and to learn how to do this.

That being said, to use this you need to have the library installed.
The things you need for the library can be found at: <http://www.netlib.org/odepack/> but similar files are provided here, in the odepack directory.
After downloading them, you can run:
```
gfortran -shared -O2 odepack*.f -o libodepack.so -fPIC 
```
to create a shared library.
This needs to be placed where Julia's `ccall()` can find it (e.g `/usr/lib/`).

## Usage
Use `@diff_eq` and `@diff_eq_jac` to construct functions needed as inputs for the solver.
To solve `ode()` function is used. 
This can be called with or without the Jacobian.

For further details look at the help from Julia REPL for the function/macros above or have a look at the tests or the code itself.

If you have any ideas or comments feel free to make an issue or contact me.
