# import needed packages
# for jlmie.jl
include("..\\src\\jlmie.jl")
using SpecialFunctions
# for this script
using Plots

# settings
# structure
nmat = 4  # refractive index of a sphere
radius = 100*1e-9

# environment
nenv = 1.00

# wavelenth 
lbd0 = (400:1000)*1e-9; lbdp = lbd0.*1e9;

# xy plot range
x = (-150:2:150)*1e-9
y = 0.0
z = (-150:2:150)*1e-9
xp = x*1e9; yp = y*1e9; zp = z*1e9;

# other settings
nmax = -1  # -1: namx large enough (determined from x)

# main
Et,Er,Ep = jlmie_Enf(nmat,radius,lbd0,nenv,x,y,z,nmax)

# save
if true
    runningfilename = splitext(splitpath(@__FILE__)[end])[1]
    outfilename = ".\\result\\" * runningfilename
    savefig(plt,outfilename)
    println("Graph has been saved in result directory")
end