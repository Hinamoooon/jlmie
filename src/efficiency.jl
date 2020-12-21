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