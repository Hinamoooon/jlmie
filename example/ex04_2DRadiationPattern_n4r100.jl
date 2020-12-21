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
# lbd0 = 618*1e-9  # ED
# lbd0 = 831*1e-9  # MD
lbd0 = 918*1e-9  # 1nd Kerker
# lbd0 = 756*1e-9  # 2nd Kerker
lbdp = lbd0.*1e9

# angular range
theta = range(0,2π,length=73) # per 5 degrees with length 37
phi = 0.0
# phi = pi/2

# other settings
nmax = -1  # -1: namx large enough (determined from x)

Isff = zeros(length(theta))
for i = 1:length(theta)
    Isff[i],_,_ = jlmie_Isff(nmat,radius,lbd0,nenv,theta[i],phi,nmax)
end
Isff = Isff/maximum(Isff)

# show results
plt = plot(theta, Isff,
           proj=:polar,
           legend=false
)
display(plt)

# save
if true
    runningfilename = splitext(splitpath(@__FILE__)[end])[1]
    outfilename = ".\\result\\" * runningfilename
    savefig(plt,outfilename)
    println("Graph has been saved in result directory")
end