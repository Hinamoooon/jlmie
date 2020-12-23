using jlmie
using Plots

n = 3
# angular range
theta = range(0,2π,length=73) # per 5 degrees with length 73
mu = cos.(theta)
pin,taun = jlmie_pt(mu,n)
# show results
# validation can be done by comparing to ITMO Mie Calculator
# https://physics.itmo.ru/en/mie
plt = plot(theta, pin, proj=:polar,
    # xlabel       = "Wavelength (nm)",
    # ylabel       = "Scattering efficiency",
    # legend       = false,
    label        = "pi",
    # xlims        =(-3,3),
    # ylims        =(-3,3),
    # aspect_ratio =0.5,
    # title        ="TITLE",
    # linecolor    =:blue,
    # linewidth    =5,
    # linestyle    =:dot,
    # size         =(400,300),
) # https://qiita.com/I_ppp/items/dca3552affa6a672e4bd
plot!(theta, taun, proj=:polar,label="tau")
display(plt)

# save
if true
    runningfilename = splitext(splitpath(@__FILE__)[end])[1]
    outfilename = ".\\result\\" * runningfilename
    savefig(plt,outfilename)
    println("Graph has been saved in result directory")
end