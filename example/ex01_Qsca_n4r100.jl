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

# other settings
nmax = -1  # -1: namx large enough (determined from x)

# main
m,x = jlmie_mx(nmat,radius,lbd0,nenv)
Qsca = jlmie_Qsca(m,x,nmax)
println("calculation completed")

# show results
# validation can be done by comparing to ITMO Mie Calculator
# https://physics.itmo.ru/en/mie
plt = plot(lbdp,Qsca,
    xlabel       = "Wavelength (nm)",
    ylabel       = "Scattering efficiency",
    legend       = false,
    # label        ="LABEL",
    # xlims        =(-3,3),
    # ylims        =(-3,3),
    # aspect_ratio =0.5,
    # title        ="TITLE",
    # linecolor    =:blue,
    # linewidth    =5,
    # linestyle    =:dot,
    # size         =(400,300),
) # https://qiita.com/I_ppp/items/dca3552affa6a672e4bd
display(plt)

# save
if true
    runningfilename = splitext(splitpath(@__FILE__)[end])[1]
    outfilename = ".\\result\\" * runningfilename
    savefig(plt,outfilename)
    println("Graph has been saved in result directory")
end