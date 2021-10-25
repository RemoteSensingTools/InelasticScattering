using InelasticScattering
using Plots
λ = 440*1.e-7 #cm
ν̃ = 1/λ
T = 250 #k

n2 = InelasticScattering.getMolecularConstants(InelasticScattering.N₂(), (0.8));
compute_effective_coefficents!(ν̃, T, n2);
compute_energy_levels!(n2);
compute_σ_Rayl_coeff!(n2);
compute_σ_Rayl_VibRaman_coeff_hires!(T, n2);
compute_σ_VibRaman_coeff!(T, n2);
compute_σ_RoVibRaman_coeff!(T, n2);

# Create grid:
grid_out = -250:0.002:250
σ_out = similar(collect(grid_out));
σ_out .= 0;

apply_lineshape!(0, n2.effCoeff.σ_Rayl_coeff, collect(grid_out), σ_out, 1, 300.0, 28);
apply_lineshape!(n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2, n2.effCoeff.σ_RoRaman_coeff_JtoJp2, collect(grid_out), σ_out, 1, 300.0, 28);
apply_lineshape!(n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2, n2.effCoeff.σ_RoRaman_coeff_JtoJm2, collect(grid_out), σ_out, 1, 300.0, 40);

plotly()
plot(grid_out, σ_out*1e50, yscale=:log10)