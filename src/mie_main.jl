# export jlmie_mx, jlmie_abcd, jlmie_pt, jlmie_S12

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


function jlmie_S12(m,x,mu,n::Int)
# calculate n-th order S1,S2 component of scattering matrix.
    an,bn,_,_ = jlmie_abcd(m,x,n)
    pin,taun = jlmie_pt(mu,n)

    S1n = (2*n+1)/(n*(n+1)).*(an.*pin.+bn.*taun)
    S2n = (2*n+1)/(n*(n+1)).*(an.*taun.+bn.*pin)
    return S1n,S2n
end


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


function largenmax(x)
    nmax = round(maximum(x.+x.^(1/3).+2))
    nmax = Int(nmax)
    # println("nmax = "*string(nmax))
    return nmax
end