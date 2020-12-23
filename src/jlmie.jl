module jlmie
    using SpecialFunctions

    include("general_functions.jl")
    include("mie_main.jl")
    export jlmie_mx, jlmie_abcd, jlmie_pt, jlmie_S12
    include("efficiency.jl")
    export jlmie_Qsca, jlmie_Qsca_a_n, jlmie_Qsca_b_n
    include("farfield.jl")
    export jlmie_Isff, jlmie_Isff_EM_n
end