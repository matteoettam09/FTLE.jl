#= 
Author(s): Matteo Manzi, Mattia Petrini
=#

using DifferentialEquations
using LinearAlgebra
using Plots; 
gr()
using GLMakie
using RecipesBase

using Distributions

function pendulum!(dz,z,α,t)
    dz[1] = z[2]
    dz[2] = (α*cos(5t)-1)*sin(z[1])
    dz_dz0 = [z[3] z[4]; z[5] z[6]]
    dg_dz = [0 1; cos(z[1])*(α*cos(5t) - 1) 0]
    ddt_dz_dz0 = dg_dz * dz_dz0
    dz[3] = ddt_dz_dz0[1, 1]
    dz[4] = ddt_dz_dz0[1, 2]
    dz[5] = ddt_dz_dz0[2, 1]
    dz[6] = ddt_dz_dz0[2, 2]
end

α = 2.5 
T = 10.0 
t_span = (0.0,T)

len = 250
FTLE = zeros((len, len))

ran = pi
x0_list = range(-ran, ran,length=len)
v0_list = range(-ran, ran,length=len)

i = 0
for x_0 in x0_list
    global i = i + 1
    j = 0
    for v_0 in v0_list
        j = j + 1
        z_0 = [x_0, v_0, 1, 0, 0, 1]
        prob = ODEProblem(pendulum!,z_0,t_span,α)
        sol = solve(prob)
        dzt_dz0_array = last(sol)[3:6]
        dzt_dz0 = [dzt_dz0_array[1] dzt_dz0_array[2]; dzt_dz0_array[3] dzt_dz0_array[4]]
        Delta = transpose(dzt_dz0)*dzt_dz0
        eig_val = eigen(Delta).values
        FTLE[j, i] = 1/T * log(sqrt(maximum(eig_val)))
    end
end

Plots.heatmap(x0_list, v0_list, FTLE)