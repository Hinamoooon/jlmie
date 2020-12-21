# jlmie

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Hinamoooon.github.io/jlmie/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Hinamoooon.github.io/jlmie/dev)
[![Build Status](https://travis-ci.com/Hinamoooon/jlmie.jl.svg?branch=master)](https://travis-ci.com/Hinamoooon/jlmie)

Julia implementation of Mie theory for nanophotonics.[[1]](#reference)  
This project is under development.  
Author is checking the validity of the implementation as much as possible, but any warranty is not provided.  

# Required packages
- `SpecialFunctions` for computation of Bessel functions
- `Plots` and `PyPlot` for visualization

Use Julia Package Manager (]) and run
```
(@v1.5) pkg> add SpecialFunctions Plots PyPlot
```

# How to use
Important functions are included in `./code/src/jlmie.jl`. After downloading or cloning jlmie, above file need to be included, for example, by
```
include "jlmie\\code\\src\\jlmie.jl"
```
Note that this expression is valid for Windows users. Mac or Linux users need to rewrite it according to your environment.  

# Example
1. import packages
```
include "..\\src\\jlmie.jl"  # correct appropreately depending on where you are
using SpecialFunctions
using Plots
```
2. define problem
```
# structure
nmat = 4  # refractive index of the material constituting a sphere
radius = 75*1e-9

# vacuum wavelenth range to be calculated
lbd0 = (400:800)*1e-9;
```
3. calculation
```
# Conversion of parameters into relative refractive index (m) and size parameter (x)
m,x = jlmie_mx(nmat,radius,lbd0)  # Environment is vacuum (nenv = 1) by default

# Scattering efficiency
Qsca = jlmie_Qsca(m,x)
```

4. plot the calculated spectrum
```
plot(lbdp,Qsca,
    xlabel       = "Wavelength (nm)",
    ylabel       = "Scattering efficiency",
    legend       = false,
)
```

Some examples are included in `./code/examples`.
- `Qsca_n4r75.jl` calculates total scattering efficiency spectrum of a high-refractive-index nanosphere (n = 4) of 75 nm in radius.
- `QscaMultipoles_n4r100.jl` calculates total scattering efficiency of a sphere with n = 4 and r = 100 nm and each contribution from n-th order electric and magnetic Mie resonance.
![Output example in jlmie](https://github.com/Hinamoooon/jlmie/blob/master/result/QscaMultipoles_n4r100.png?raw=true)

# License
jlmie is distributed under MIT license.

# Reference
1. Absorption and Scattering of Light by Small Particles; Bohren, C. F., Huffman, D. R., Eds.; Wiley-VCH Verlag GmbH: Weinheim, Germany, 1998.