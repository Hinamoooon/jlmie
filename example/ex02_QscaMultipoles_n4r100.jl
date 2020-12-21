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
Qsca_ED = jlmie_Qsca_a_n(m,x,1)
Qsca_MD = jlmie_Qsca_b_n(m,x,1)
Qsca_EQ = jlmie_Qsca_a_n(m,x,2)
Qsca_MQ = jlmie_Qsca_b_n(m,x,2)
Qsca_EH = jlmie_Qsca_a_n(m,x,3)
Qsca_MH = jlmie_Qsca_b_n(m,x,3)
println("calculation completed")

# show results
# validation can be done by comparing to ITMO Mie Calculator
# https://physics.itmo.ru/en/mie
plt = plot(lbdp,Qsca,
    xlabel       = "Wavelength (nm)",
    ylabel       = "Scattering efficiency",
    # legend       = false,
    label        = "Total",
    # xlims        =(-3,3),
    # ylims        =(-3,3),
    # aspect_ratio =0.5,
    # title        ="TITLE",
    # linecolor    =:blue,
    # linewidth    =5,
    # linestyle    =:dot,
    # size         =(400,300),
) # https://qiita.com/I_ppp/items/dca3552affa6a672e4bd
plot!(lbdp,Qsca_ED,
    label = "ED",
    linestyle = :dash,
)
plot!(lbdp,Qsca_MD,
    label = "MD",
    linestyle = :dash,
)
plot!(lbdp,Qsca_EQ,
    label = "EQ",
    linestyle = :dash,
)
plot!(lbdp,Qsca_MQ,
    label = "MQ",
    linestyle = :dash,
)
plot!(lbdp,Qsca_EH,
    label = "EH",
    linestyle = :dash,
)
plot!(lbdp,Qsca_MH,
    label = "MH",
    linestyle = :dash,
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
lbdp_peaks = zeros(6)
lbdp_peaks[1] = lbdp[findmax(Qsca_ED)[2]]
lbdp_peaks[2] = lbdp[findmax(Qsca_MD)[2]]
lbdp_peaks[3] = lbdp[findmax(Qsca_EQ)[2]]
lbdp_peaks[4] = lbdp[findmax(Qsca_MQ)[2]]
lbdp_peaks[5] = lbdp[findmax(Qsca_EH)[2]]
lbdp_peaks[6] = lbdp[findmax(Qsca_MH)[2]]

println("maximum of ED is at " * string(round(lbdp_peaks[1])))
println("maximum of MD is at " * string(round(lbdp_peaks[2])))
println("maximum of EQ is at " * string(round(lbdp_peaks[3])))
println("maximum of MQ is at " * string(round(lbdp_peaks[4])))
println("maximum of EH is at " * string(round(lbdp_peaks[5])))
println("maximum of MH is at " * string(round(lbdp_peaks[6])))