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
theta = range(0,π,length=37) # per 5 degrees with length 37
phi = range(0,2π,length=73) # per 5 degrees with length 73

# other settings
nmax = -1  # -1: namx large enough (determined from x)

Isff = zeros(length(theta),length(phi))
for i = 1:length(theta)
    for j = 1:length(phi)
        Isff[i,j],_,_ = jlmie_Isff(nmat,radius,lbd0,nenv,theta[i],phi[j],nmax)
    end
end
Isff = Isff/maximum(Isff)
mat_x = Isff.* (sin.(theta) * cos.(phi)')
mat_y = Isff.* (sin.(theta) * sin.(phi)')
mat_z = Isff.* (cos.(theta) * cos.(phi*0)')

# Just a sphere
# mat_x = 1 .* (sin.(theta) * cos.(phi)')
# mat_y = 1 .* (sin.(theta) * sin.(phi)')
# mat_z = 1 .* (cos.(theta) * cos.(phi*0)')

# show results
pyplot()
plt = plot(mat_x,mat_y,mat_z,
           st=:surface,
           camera=(45,30), # azimuth, elevate
           xlabel="x",
           ylabel="y",
           zlabel="z",
           xlims=(-1,1),
           ylims=(-1,1),
           zlims=(-1,1),
)
display(plt)

## Followings are plots with Plots
# plt = plot(theta,phi,mat_z', st=:surface,
#            xlims = (-5,5),
#            ylims = (-5,5),
#            zlims = (-5,5),
# )
# plt = plot(theta,phi,mat_z', marker=:circle )
# plt = plot(theta,phi,mat_z')
# plt = plot(theta,phi,Isff')
# plt = plot(theta,phi,mat_z', st=:heatmap)

# save
if true
    runningfilename = splitext(splitpath(@__FILE__)[end])[1]
    outfilename = ".\\result\\" * runningfilename
    savefig(plt,outfilename)
    println("Graph has been saved in result directory")
end