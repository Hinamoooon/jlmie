![jlmie](https://repository-images.githubusercontent.com/323283354/7572ec00-4553-11eb-8f5d-2978b93331aa)

# jlmie
<!--
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Hinamoooon.github.io/jlmie/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Hinamoooon.github.io/jlmie/dev)
-->
[![CI](https://github.com/Hinamoooon/jlmie/workflows/CI/badge.svg)](https://github.com/Hinamoooon/jlmie/actions?query=workflow%3ACI)
[![Build Status](https://travis-ci.com/Hinamoooon/jlmie.svg?branch=master)](https://travis-ci.com/Hinamoooon/jlmie)

Julia implementation of Mie theory for nanophotonics.[[1]](#reference)  
Author is checking the validity of the implementation as much as possible, but any warranty is not provided.  
This project is under development. Because jlmie is just a transported program from author's homebuilt package written in different language, its implementation might not be optimized for Julia language...

# Installation
Use Julia Package Manager (]) and run
```
(@v1.5) pkg> add https://github.com/Hinamoooon/jlmie.git
```

# Dependency
jlmie depends on a following package.
- `SpecialFunctions` for computation of Bessel functions

# Recommendation
Installation of following packages are recommended to visualize calculated results. Example codes in `./example` require prior installation of them.
- `Plots` for visualization
- `PyPlot` for visulization of 3D graphs

Use Julia Package Manager (]) and run
```
(@v1.5) pkg> add Plots PyPlot
```

# How to use
[To be updated once it is published in the Julia's public repositry]  
Install jlmie following the [Installation section](#installation).
1. import packages
```
using jlmie
using Plots
```
2. define problem
```
# structure
nmat = 4  # refractive index of the material constituting a sphere
radius = 75*1e-9

# vacuum wavelenth range to be calculated
lbd0 = (400:800)*1e-9;
lbdp = lbdp*1e9;  # used for plot
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

# Examples
Some example codes are included in `./example`.
- `ex01_Qsca_n4r100.jl` calculates total scattering efficiency spectrum of a high-refractive-index nanosphere (n = 4) of 100 nm in radius.  
![Output example in jlmie](https://github.com/Hinamoooon/jlmie/blob/master/result/ex01_Qsca_n4r100.png?raw=true)
- `ex02_QscaMultipoles_n4r100.jl` calculates total scattering efficiency of a sphere with n = 4 and r = 100 nm and each contribution from n-th order electric and magnetic Mie resonance.  
![Output example in jlmie](https://github.com/Hinamoooon/jlmie/blob/master/result/ex02_QscaMultipoles_n4r100.png?raw=true)
- `ex03_ForwardBackward_n4r100.jl` calculates forward (θ = 0°) and backward (θ = 180°) scattering intensities.  
![Output example in jlmie](https://github.com/Hinamoooon/jlmie/blob/master/result/ex03_ForwardBackward_n4r100.png?raw=true)
- `ex04_2DRadiationPattern_n4r100.jl` calculates far-field radiation patterns at a fixed wavelength.  
![Output example in jlmie](https://github.com/Hinamoooon/jlmie/blob/master/result/ex04_2DRadiationPattern_n4r100.png?raw=true)
- `ex05_3DRadiationPattern_n4r100.jl` is a 3D version of ex04 which calculates radiation pattern.  
![Output example in jlmie](https://github.com/Hinamoooon/jlmie/blob/master/result/ex05_3DRadiationPattern_n4r100.png?raw=true)

# Directory structure [To be updated]
Important functions are included in `./src`. If you want to see source codes of any functions, see `./src/jlmie.jl`. 

# Author
Tatsuki Hinamoto@Kobe University, Japan

# License
jlmie is distributed under MIT license.

# Reference
1. Absorption and Scattering of Light by Small Particles; Bohren, C. F., Huffman, D. R., Eds.; Wiley-VCH Verlag GmbH: Weinheim, Germany, 1998.

# To do list
- read wavelength dependent refractive indices
- Computation of near-field profiles
- Detailed documents describing theoretical backgrounds