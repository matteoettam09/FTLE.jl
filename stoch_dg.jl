
using DifferentialEquations

using Plots;
using Statistics;
using Distributions;
using LinearAlgebra;

#use GR module
gr();

A = 0.25
ω = 2π/10

piA = π*A
ppA = π^2*A

α = 0.0;
ϵ = 0.25;

function doublegyre(dz, z, p, t)
  a = ϵ * sin(ω * t)
  b = 1 - 2 * a
  f = a * z[1]^2 + b * z[1]
  spf = sin(π * f)
  cpf = cos(π * f)
  spy = sin(π * z[2])
  cpy = cos(π * z[2])

  df_dx = (2 * a * z[1] + b)

  dz[1] = -piA * spf * cpy
  dz[2] = piA * cpf * spy * df_dx
end

function σ_DJ(du, u, p, t)
  du[1] = α # * u[1]
  du[2] = α # * u[2]
end

Δt = 0.05

u0 = [0.1, 0.6];
T = 40.0;
tspan = (0.0, T);

prob_sde_DJ = SDEProblem(doublegyre, σ_DJ, u0, tspan)
sol = solve(prob_sde_DJ, EM(), dt=Δt)

plot(sol, vars=(1, 2))
savefig("S_DJ.png")