# copy following lines to import required packages
## for jlmie.jl
# include("jlmie\\code\\src\\jlmie.jl")
# using SpecialFunctions

function jlmie_mx(nmat,radius,lbd0,nenv=1.00)
# conversion of parameters for calculation of light scattering
# by a small sphere based on Mie theory
# Input:
# nmat: n_material; refractive index of the material of the sphere
# radius: radius of the sphere
# lbd0: target vacuum wavelength (range)
# nenv: n_environement; refractive index of the environment
#       around the sphere (default: air)
# Output:
# m: relative refractive index (n_material / n_environment)
# x: size parameter (wavenumber * radius)
    m = nmat./nenv
    lbd = lbd0./nenv  # effective wavelength in the environment
    k = 2*pi./lbd  # free space wavenumber in the environment
    x = k.*radius
    return m,x
end

function jlmie_Qsca(m,x,nmax::Int=-1)
# Calculate scattering efficiency of a sphere
# Input:
# m: relative refractive index (n_material / n_environment)
# x: size parameter (wavenumber * radius)
# namx: maximun order of resonances being considered
#       (nmax = 1: dipole, 2: dipole+quadrupole, ...,
#        -1: sufficiently large value determined by x)
# Output:
# Qsca: Scattering efficiency (total up to n-th order resonances)
    if nmax == -1
        nmax = largenmax(x)
    end

    Qsca = zeros(length(x))
    for i = 1:nmax
        Qsca = Qsca .+ jlmie_Qsca_n(m,x,i)
    end
    return Qsca
end

function jlmie_Qsca_n(m,x,n::Int)
    an,bn,_,_ = jlmie_abcd(m,x,n)
    nth_Qsca = 2 ./ x.^2 .* (2*n+1).*(abs.(an).^2 .+ abs.(bn).^2)
    return nth_Qsca
end

function jlmie_Qsca_a_n(m,x,n::Int)
# calculate scattering efficiency due to the n-th order electric multipole
# (n = 1: dipole, 2: quadrupole, ...)
    an,_,_,_ = jlmie_abcd(m,x,n)
    nth_Qsca_a = 2 ./ x.^2 .* (2*n+1).*(abs.(an).^2)
    return nth_Qsca_a
end

function jlmie_Qsca_b_n(m,x,n::Int)
# calculate scattering efficiency due to the n-th order magnetic multipole
# (n = 1: dipole, 2: quadrupole, ...)
    _,bn,_,_ = jlmie_abcd(m,x,n)
    nth_Qsca_b = 2 ./ x.^2 .* (2*n+1).*(abs.(bn).^2)
    return nth_Qsca_b
end

function jlmie_Enf(nmat,radius,lbd0,nenv,x,y,z,nmax::Int=-1)
    if nmax == -1
        nmax = largenmax(jlmie_mx(nmat,radius,lbd0,nenv)[2])
    end

end

function jlmie_Enf1_n(nmat,radius,lbd0,nenv,r,theta,phi,n::Int=1)
# implementation of eq.(4.40) of B&H
    m = nmat./nenv
    lbd = lbd0./nenv  # effective wavelength in the environment
    lbd1 = lbd0./nmat  # effective wavelength in the sphere
    k = 2*pi./lbd  # free space wavenumber in the environment
    k1 = 2*pi./lbd1  # free space wavenumber in the sphere
    x = k.*radius
    x1 = k1.*radius
    mu = cos.(theta)
    _,_,cn,dn = jlmie_abcd(m,x,n)
    pin,taun = jlmie_pt(cos.(theta),n)
    
    jx1 = sbesselj.(n,x1)
    d1x1jx1 = x1.*sbesselj.(n-1,x1) .- n.*sbesselj.(n,x1)
    # d1mxjmx = m.*x.*sbesselj.(n-1,m.*x) .- n.*sbesselj.(n,m.*x)

    E0 = 1.0;
    En = (1im)^n *E0*(2*n+1)/(n*(n+1))
    
    M1o1n_r = 0
    M1o1n_t = cos.(phi).*pin.*jx1
    M1o1n_p = -1*sin.(phi).*taun.*jx1
    # M1e1n_r = 0
    # M1e1n_t = -1*sin.(phi).*pin.*jx1
    # M1e1n_p = -1*cos.(phi).*taun.*jx1
    # N1o1n_r = sin.(phi).*n.*(n+1).*sin.(theta).*pin.*jx1./x1
    # N1o1n_t = sin.(phi).*tau.*d1x1jx1./x1
    # N1o1n_p = cos.(phi).*pin.*d1x1jx1./x1
    N1e1n_r = cos.(phi).*n.*(n+1).*sin.(theta).*pin.*jx1./x1
    N1e1n_t = cos.(phi).*tau.*d1x1jx1./x1
    N1e1n_p = -1*sin.(phi).*pin.*d1x1jx1./x1

    E1r_n = En.*(cn.*M1o1n_r .- 1im.*dn.*N1e1n_r)
    E1t_n = En.*(cn.*M1o1n_t .- 1im.*dn.*N1e1n_t)
    E1p_n = En.*(cn.*M1o1n_p .- 1im.*dn.*N1e1n_p)
    return E1r_n,E1t_n,E1p_n
end

function jlmie_Enfs_n(nmat,radius,lbd0,nenv,r,theta,phi,n::Int=1)
    # implementation of eq.(4.45) of B&H
    m = nmat./nenv
    lbd = lbd0./nenv  # effective wavelength in the environment
    k = 2*pi./lbd  # free space wavenumber in the environment
    x = k.*radius
    mu = cos.(theta)
    an,bn,_,_ = jlmie_abcd(m,x,n)
    pin,taun = jlmie_pt(cos.(theta),n)

    hx = shankelh1.(n,x)
    d1xhx = x.*shankelh1.(n-1,x) .- n.*shankelh1.(n,x)

    E0 = 1.0;
    En = (1im)^n *E0*(2*n+1)/(n*(n+1))
    
    M3o1n_r = 0
    M3o1n_t = cos.(phi).*pin.*hx
    M3o1n_p = -1*sin.(phi).*taun.*hx
    # M3e1n_r = 0
    # M3e1n_t = -1*sin.(phi).*pin.*hx
    # M3e1n_p = -1*cos.(phi).*taun.*hx
    # N3o1n_r = sin.(phi).*n.*(n+1).*sin.(theta).*pin.*hx./x
    # N3o1n_t = sin.(phi).*tau.*d1xhx./x
    # N3o1n_p = cos.(phi).*pin.*d1xhx./x
    N3e1n_r = cos.(phi).*n.*(n+1).*sin.(theta).*pin.*hx./x
    N3e1n_t = cos.(phi).*tau.*d1xhx./x
    N3e1n_p = -1*sin.(phi).*pin.*d1xhx./x

    Esr_n = En.*(1im.*an.*N3e1n_r .- bn.*M3o1n_r)
    Est_n = En.*(1im.*an.*N3e1n_t .- bn.*M3o1n_t)
    Esp_n = En.*(1im.*an.*N3e1n_p .- bn.*M3o1n_p)
    return Esr_n,Est_n,Esp_n
end

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

function jlmie_Isff_EM_n(nmat,radius,lbd0,nenv,theta,phi,n)
# calculate far-feild intensity of scattered field
# from n-th order Electric and Magnetic multipoles separately
    m = nmat./nenv
    lbd = lbd0./nenv  # effective wavelength in the environment
    k = 2*pi./lbd  # free space wavenumber in the environment
    x = k.*radius
    mu = cos.(theta)
    an,bn,_,_ = jlmie_abcd(m,x,n)
    pin,taun = jlmie_pt(mu,n)
    S1EDn = (2*n+1)/(n*(n+1)).*(an.*pin)
    S2EDn = (2*n+1)/(n*(n+1)).*(an.*taun)
    S1MDn = (2*n+1)/(n*(n+1)).*(bn.*taun)
    S2MDn = (2*n+1)/(n*(n+1)).*(bn.*pin)
    # settings
    E0 = 1  # amplitude of incident electric field
    r = 1  # distance from the sphere to an observation point (1 m)
    # main
    Est_EDn = E0 * exp.(1im.*k*r)./(-1im.*k*r) .* cos.(phi).*S2EDn
    Esp_EDn = -1*E0 * exp.(1im.*k*r)./(-1im.*k*r) .* sin.(phi).*S1EDn
    Est_MDn = E0 * exp.(1im.*k*r)./(-1im.*k*r) .* cos.(phi).*S2MDn
    Esp_MDn = -1*E0 * exp.(1im.*k*r)./(-1im.*k*r) .* sin.(phi).*S1MDn
    Is_EDn = abs.(Est_EDn).^2 .+ abs.(Esp_EDn).^2
    Is_MDn = abs.(Est_MDn).^2 .+ abs.(Esp_MDn).^2
    return Is_EDn, Is_MDn
end

function jlmie_Isff(nmat,radius,lbd0,nenv,theta,phi,nmax::Int=-1)
# calculate far-field intensity of scattered light
# (as a total from n = 1 to nmax)
    if nmax == -1
        nmax = largenmax(jlmie_mx(nmat,radius,lbd0,nenv)[2])
    end

    Est,Esp = jlmie_Esff_n(nmat,radius,lbd0,nenv,theta,phi,1)
    if nmax > 1
        for n = 2:nmax
            Est_n,Esp_n = jlmie_Esff_n(nmat,radius,lbd0,nenv,theta,phi,n)
            Est = Est .+ Est_n
            Esp = Esp .+ Esp_n
        end
    end
    Ist = abs.(Est).^2
    Isp = abs.(Esp).^2
    Is = Ist .+ Isp
    return Is,Ist,Isp
end

function jlmie_Esff_n(nmat,radius,lbd0,nenv,theta,phi,n::Int)
# calculate scattered electric far-fields E_s,theta and E_s,phi
# resulting from n-th order resonances
# See Sec. 4.4.4 of B&H for derivation
# Input:
# nmat: n_material; refractive index of the material of the sphere
# radius: radius of the sphere
# lbd0: target vacuum wavelength (range)
# nenv: n_environement; refractive index of the environment
#       around the sphere
    # preparation
    m = nmat./nenv
    lbd = lbd0./nenv  # effective wavelength in the environment
    k = 2*pi./lbd  # free space wavenumber in the environment
    x = k.*radius
    mu = cos.(theta)
    S1n,S2n = jlmie_S12(m,x,mu,n)

    # settings
    E0 = 1  # amplitude of incident electric field
    r = 1  # distance from the sphere to an observation point (1 m)

    # main
    Est_n = E0 * exp.(1im.*k*r)./(-1im.*k*r) .* cos.(phi).*S2n
    Esp_n = -1*E0 * exp.(1im.*k*r)./(-1im.*k*r) .* sin.(phi).*S1n
    return Est_n,Esp_n
end

function jlmie_S12(m,x,mu,n::Int)
# calculate n-th order S1,S2 component of scattering matrix.
    an,bn,_,_ = jlmie_abcd(m,x,n)
    pin,taun = jlmie_pt(mu,n)

    S1n = (2*n+1)/(n*(n+1)).*(an.*pin.+bn.*taun)
    S2n = (2*n+1)/(n*(n+1)).*(an.*taun.+bn.*pin)
    return S1n,S2n
end

function jlmie_pt(mu,n::Int)
# See Sec. 4.3.1 of B&H
# ref of pi,tau: [D. Deirmendjian, “Electromagnetic Scattering on Spherical Polydispersions"]
    ck = sum(mu .< -1) .+ sum(mu .> 1)
    if ck > 0
        error("mu is cos(theta), -1 <= mu <= 1")
    end

    if n == 0
        pin = 0*mu  # pi(0,mu) = 0
        taun = []  # not used but tau(0,mu) = 0
    elseif n == 1
        pin = 0*mu.+1  # pi(1,mu) = 1
        taun = n.*mu.*pin  # note that pi(0,mu) = 0
    elseif n > 1
        # recurrence, eq.(4.47) of B&H
        pin = (2*n-1)/(n-1).*mu.*jlmie_pt(mu,n-1)[1] .-
              (n)/(n-1).*jlmie_pt(mu,n-2)[1]
        taun = n.*mu.*pin .- (n+1).*jlmie_pt(mu,n-1)[1] 
    else
        error("n should be larger than 0")
    end
    return pin,taun
end

function jlmie_abcd(m,x,n::Int)
# calculate n-th order Mie coefficients a,b,c,d
# input:
# m: relative refractive index (n_material / n_environment)
# x: size parameter (wavenumber * radius)

    # Bessel/Hankel functions
    jx = sbesselj.(n,x)
    jmx = sbesselj.(n,m.*x)
    hx = shankelh1.(n,x)

    # derivative
    d1xjx = x.*sbesselj.(n-1,x) .- n.*sbesselj.(n,x)
    d1mxjmx = m.*x.*sbesselj.(n-1,m.*x) .- n.*sbesselj.(n,m.*x)
    d1xhx = x.*shankelh1.(n-1,x) .- n.*shankelh1.(n,x)
    
    an = (m.^2 .*jmx.*d1xjx .- jx.*d1mxjmx) ./
         (m.^2 .*jmx.*d1xhx .- hx.*d1mxjmx)
    bn = (jmx.*d1xjx .- jx.*d1mxjmx) ./
         (jmx.*d1xhx .- hx.*d1mxjmx)
    cn = (jx.*d1xhx .- hx.*d1xjx) ./
         (jmx.*d1xhx .- hx.*d1mxjmx)
    dn = (m.*jx.*d1xhx .- m.*hx.*d1xjx) ./
         (m.^2 .*jmx.*d1xhx .- hx.*d1mxjmx)
    return an,bn,cn,dn
end

function largenmax(x)
    nmax = round(maximum(x.+x.^(1/3).+2))
    nmax = Int(nmax)
    # println("nmax = "*string(nmax))
    return nmax
end

# spherical bessel functions
sbesselj(nu, x) = √(π./2x)*besselj.(nu+1/2, x)
sbessely(nu, x) = √(π./2x)*bessely.(nu+1/2, x)
shankelh1(nu,x) = sbesselj.(nu,x) + 1im.*sbessely.(nu,x)