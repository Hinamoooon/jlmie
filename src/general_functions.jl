# spherical bessel functions
sbesselj(nu, x) = √(π./2x)*besselj.(nu+1/2, x)
sbessely(nu, x) = √(π./2x)*bessely.(nu+1/2, x)
shankelh1(nu,x) = sbesselj.(nu,x) + 1im.*sbessely.(nu,x)

function cart2d2sph2d(xyz,ax_horz,ax_vert)
# generate 2D meshgrid in spherical coordinates system
# corresponding to Cartesian 2D profile (xz, yz, xy)
# for calculation of field profiles
# Input:
# xyz: 1->xz, 2->yz, 3->xy
# ax_horz, ax_vert: 1D array of axes
    mat1 = ax_horz' .* ones(length(ax_vert))
    mat2 = ones(length(ax_horz))' .* ax_vert

    if xyz == 1
        # xz plot range
        mat_r,mat_theta,mat_phi = cart2sph(mat1,0,mat2)
    elseif xyz == 2
        # yz plot range
        mat_r,mat_theta,mat_phi = cart2sph(0,mat1,mat2)
    elseif xyz == 3
        # xy plot range
        mat_r,mat_theta,mat_phi = cart2sph(mat1,mat2,0)
    end
    return mat_r,mat_theta,mat_phi
end


function cart2sph(x,y,z)
# convert Cartesian coordinates(x,y,z) to spherical coordinates(r,θ,φ)
# θ is a polar angle (θ=0 is z-axis), φ is an azimuthal angle (φ=0 is x-axis)
    r = sqrt.(x.^2 .+ y.^2 .+ z.^2)
    φ = mod.(atan.(y,x) .+ 2*π, 2*π)
    θ = acos.(z./r)
    return r,θ,φ
end