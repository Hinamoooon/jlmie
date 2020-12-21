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
lbd0 = (400:1:1000)*1e-9; lbdp = lbd0.*1e9;

# angular range
theta = [0 pi]
phi = [0 0]

# other settings
nmax = -1  # -1: namx large enough (determined from x)

Isff_F,_,_ = jlmie_Isff(nmat,radius,lbd0,nenv,theta[1],phi[1],nmax)
Isff_B,_,_ = jlmie_Isff(nmat,radius,lbd0,nenv,theta[2],phi[2],nmax)

# show results
# validation can be done by comparing to ITMO Mie Calculator
# https://physics.itmo.ru/en/mie
plt = plot(lbdp,Isff_F,
    xlabel       = "Wavelength (nm)",
    ylabel       = "Scattering efficiency",
    # legend       = false,
    label        = "Forward",
    # xlims        =(-3,3),
    # ylims        =(-3,3),
    # aspect_ratio =0.5,
    # title        ="TITLE",
    # linecolor    =:blue,
    # linewidth    =5,
    # linestyle    =:dot,
    # size         =(400,300),
) # https://qiita.com/I_ppp/items/dca3552affa6a672e4bd
plot!(lbdp,Isff_B,
    label = "Backward",
)
display(plt)

# save
if true
    runningfilename = splitext(splitpath(@__FILE__)[end])[1]
    outfilename = ".\\result\\" * runningfilename
    savefig(plt,outfilename)
    println("Graph has been saved in result directory")
end

# find peaks
lbdp_peaks = zeros(2)
lbdp_peaks[1] = lbdp[findmax(Isff_F./Isff_B)[2]]
lbdp_peaks[2] = lbdp[200+findmax(Isff_B[201:end]./Isff_F[201:end])[2]]  # avoid high order modes

println("maximum of F/B ratio is at " * string(round(lbdp_peaks[1])))
println("maximum of B/F ratio is at " * string(round(lbdp_peaks[2])))