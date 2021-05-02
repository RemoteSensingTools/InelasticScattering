abstract type AbstractMolecule end

struct N₂ <: AbstractMolecule end
struct O₂ <: AbstractMolecule end
struct H₂ <: AbstractMolecule end

"""
    struct PolTensor{FT}

A struct which provides all polarizability tensor elements 
    relevant to Rayleigh, and rotational+vibrational Raman cross-sections  
    Ref: Buldakov et al.,(1999), Molecular Spectroscopy

    The polarizability at arbitrary wavenumbers [m⁻¹] and temperatures [K] can be computed as 
    
    α̅(2πcν, T) = α̅₀₀(1 + α_b⋅T + α_c⋅T²)/(1-(2πcν/ω₀)²

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PolarizationTensor{FT}
    "Average Polarizability `[m³]` at an angular frequency of ω₀ `[Hz]` and T=0K"
    α̅₀₀::FT        
    "First derivative of `α̅₀₀` wrt the internuclear distance rₑ `[m²]`"
    α_prime::FT  
    "Reference frequency `[Hz]`"
    ω₀::FT     
    "Linear T-dependence coefficient `[K⁻¹]`"
    α_b::FT      
    "Quadratic T-dependence coefficient `[K⁻²]`"
    α_c::FT
    "Polarizability anisotropy `[m³]` " 
    γ̅ ::FT
    "First derivative of `γ̅` with respect to the internuclear distance `rₑ` `[m²]`"
    γ_prime::FT       
end

"""
    struct DunhamCoefficients{FT}

A struct which provides 
the Dunham expansion coefficients to compute the molecular energy levels 
corresponding to different rotational (J) and vibrational (v) states. The overall state 
of the molecule is given by its electronic ground state. #Constants from "Molecular 
Spectra and Molecular Structure IV. Constants of diatomic molecules" by K.P. Huber 
and G. Herzberg, 1978 (https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Mask=1000)

Needs description of Y elements

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct DunhamCoefficients{FT}
    "Dunham Expansion Coefficients (5x5 Matrix)"
    Y::AbstractArray{FT}    
end

Base.@kwdef struct MolecularConstants{FT}
    vmr::FT
    PolTensor::PolarizationTensor{FT}
    DunhamCoeffs::DunhamCoefficients{FT}
end







