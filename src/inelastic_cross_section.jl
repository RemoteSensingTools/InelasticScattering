function σ_Rayleigh(ν, molecules::Array{MolecularConstants{FT}}) where {FT}
    Q=FT(0) 
    for mol in molecules
        @unpack α̅₀₀, γ̅ = mol.PolTensor
        @unpack vmr = mol
        F_King = 1 + 2(γ̅ / 3α̅₀₀)^2
        Q_Rayl =  (128/3)π^5 * α̅₀₀^2 * F_King * ν^4
        Q += vmr * Q_Rayl
    end
    return Q
end