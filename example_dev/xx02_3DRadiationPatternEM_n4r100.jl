# This code separately visualize the 3D radiation pattern
# of electric and magnetic resonances.
# Because the amplitude of scattering intensity is normalized, 
# the setting of wavelength is not very important...

using jlmie
using Plots

# settings
# structure
nmat = 3.5  # refractive index of a sphere
radius = 100*1e-9

# environment
nenv = 1.00

# wavelenth (not very important)
lbd0ED = 618*1e-9  # ED
lbd0MD = 831*1e-9  # MD
lbd0 = [lbd0ED lbd0MD]

# angular range
theta = range(0,π,length=37) # per 5 degrees with length 37
phi = range(0,2π,length=73) # per 5 degrees with length 73

# other settings
n = 1
EorH = 2  # 1: Eelectric, 2: Magnemtic
Isff = zeros(length(theta),length(phi))
for i = 1:length(theta)
    for j = 1:length(phi)
        Isff[i,j] = jlmie_Isff_EM_n(nmat,radius,lbd0[EorH],nenv,theta[i],phi[j],n)[EorH]
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
           camera=(30,30), # azimuth, elevate
           xlabel="x",
           ylabel="y",
           zlabel="z",
           xlims=(-1,1),
           ylims=(-1,1),
           zlims=(-1,1),
)
display(plt)

# save
if true
    runningfilename = splitext(splitpath(@__FILE__)[end])[1]
    outfilename = ".\\result\\" * runningfilename
    savefig(plt,outfilename)
    println("Graph has been saved in result directory")
end