# import needed packages
# for jlmie.jl
include("..\\src\\jlmie.jl")
using SpecialFunctions
# for this script
using Plots

# # xz plot range
# x = (-150:5:150)
# y = 0.0
# z = (-150:5:150)
# # calculation using broadcast
# r, theta, phi = cart2sph(x',y,z)
# # heatmap(x,z,r)
# # heatmap(x,z,theta)
# # heatmap(x,z,phi)

# meshgrid type calculation
# mat_x = x' .* ones(length(z))
# mat_z = ones(length(x))' .* z
# mat_r,mat_theta,mat_phi = cart2sph(mat_x,y,mat_z)
# calculation using broadcast
# r, theta, phi = cart2sph(x',y,z)
# heatmap(x,z,r)
# heatmap(x,z,theta)
# heatmap(x,z,phi)

x = (-150:1:150)
y = (-150:1:150)
z = (-150:1:150)

xyz = 2  # 1:x, 2:y, 3:z
if xyz == 1
    # xz plot range
    mat_x = x' .* ones(length(z))
    mat_z = ones(length(x))' .* z
    mat_r,mat_theta,mat_phi = cart2sph(mat_x,0,mat_z)
    p1 = heatmap(x,z,mat_r)
    p2 = heatmap(x,z,mat_theta)
    p3 = heatmap(x,z,mat_phi)
elseif xyz == 2
    # yz plot range
    mat_y = y' .* ones(length(y))
    mat_z = ones(length(y))' .* z
    mat_r,mat_theta,mat_phi = cart2sph(0,mat_y,mat_z)
    p1 = heatmap(y,z,mat_r)
    p2 = heatmap(y,z,mat_theta)
    p3 = heatmap(y,z,mat_phi)
elseif xyz == 3
    # xy plot range
    mat_x = x' .* ones(length(y))
    mat_y = ones(length(x))' .* y
    mat_r,mat_theta,mat_phi = cart2sph(mat_x,mat_y,0)
    p1 = heatmap(x,y,mat_r)
    p2 = heatmap(x,y,mat_theta)
    p3 = heatmap(x,y,mat_phi)
end
plot(p1,p2,p3)