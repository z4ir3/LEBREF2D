# LEBREF2D 1.0

## A MATLAB 2D mesh refinement package

LEBREF2D is a MATLAB package for the mesh refinement of 2D domains based on an efficient usage of MATLAB built-in functions and vectorization. 

<img src="/docs/lebref-pic.png" height="400" width="600">

It includes a uniform mesh refinement routine as well as routines for local adaptive mesh refinements based on the longest edge bisection (LEB) algorithm as particular case of the newest vertex bisection (NVB). 

Functions can be easily adapted to be used in an adaptive Finite Element Code for solving Partial Differential Equations (PDEs) in the setting of the *a posteriori* error estimation. 

Run the adaptive demo [lebdemo.m](src/lebdemo.m) for a demonstration of adaptive mesh refinements of some built-in domains (square, L-shaped, crack) and please read the short [documentation](docs/doc.pdf) for information.

### Compatibility 

Tested with **MATLAB 8.4.0.150421 (R2014b)**

### License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
