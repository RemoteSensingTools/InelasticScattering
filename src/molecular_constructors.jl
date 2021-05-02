"Construct constants for N₂:"
function getMolecularConstants(::N₂, vmr::FT) where {FT}
    @info "Constructing N₂ constants, using VMR (unitless) of " vmr
    @assert 0 ≤ vmr ≤ 1 "VMR has to be between [0,1]"
    p = PolarizationTensor{FT}(
            α̅₀₀     = 1.7406e-30,
            α_prime = 1.86e-30,
            ω₀      = 2.6049e16, 
            α_b     = 1.8e-6,
            α_c     = 0.0,
            γ̅       = 0.691e-30,  #Ref: Asawaroengchai & Rosenblatt (1980), J. Chem. Phys. (gamma is written as beta_e in the paper) );
            γ_prime = 1.40e-30)
    
    Y = zeros(FT, 5,5)
    # Fill in all non-zero elements:
    Y[1,2] =  1.99824  # rotational constant in equilibrium position (cm^{-1})
    Y[1,3] = -5.76E-6  # centrifugal distortion constant (cm^{-1})
    Y[2,1] =  2358.57  # vibrational constant - first term (cm^{-1}) 
    Y[2,2] = -0.017318 # rotational constant - first term (cm^{-1})
    Y[3,1] = -14.324
    Y[4,1] = -2.26e-3  #vibrational constant - third term (cm^{-1})
    MolecularConstants(vmr=vmr,PolTensor=p, DunhamCoeffs=DunhamCoefficients(Y) )
end