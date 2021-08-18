module InelasticScattering

    using DocStringExtensions
    using Parameters
    using OffsetArrays
    #import PhysicalConstants.CODATA2018:c_0, h, k_B

    include("raman_constants.jl")
    include("molecular_constructors.jl")
    include("inelastic_cross_section.jl")

    "Speed of light in vacuum in `[cm/s]`"
    global const c = 2.99792458e10
    "Planck constant in `[erg ̇s]`"
    global const h = 6.62607015e-27
    "Boltzmann constant in `[erg/K]`"
    global const k_B = 1.380649e-16

    export compute_effective_coefficents!
    export compute_σ_Rayl_coeff!
    export compute_σ_Rayl_VibRaman_coeff_hires!
    export compute_energy_levels!

    #export compute_σ_Raman_coeff!
    export compute_σ_VibRaman_coeff!, compute_σ_RoVibRaman_coeff!
end # module
