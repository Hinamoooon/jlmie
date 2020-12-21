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