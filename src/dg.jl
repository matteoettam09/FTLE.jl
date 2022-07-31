#= 
Author(s): Matteo Manzi
=#

using DifferentialEquations
using LinearAlgebra
using Plots; gr()
using PlotlyJS
using GLMakie
using RecipesBase

A = 0.1
ω = 2π/10
ϵ = 0.1

piA = π*A
ppA = π^2*A

function doublegyre!(dz,z,ϵ,t)
    a = ϵ * sin(ω*t)
    b = 1 - 2*a
    f = a*z[1]^2 + b*z[1]
    spf = sin(π*f)
    cpf = cos(π*f)
    spy = sin(π*z[2])
    cpy = cos(π*z[2])

    df_dx = (2*a*z[1] + b)

    dz[1] = - piA*spf*cpy
    dz[2] = piA*cpf*spy*df_dx

    dz_dz0 = [z[3] z[4]; z[5] z[6]]

    dg_dz_11 = - ppA * cpy * cpf * df_dx
    dg_dz_12 = ppA* spf * spy
    dg_dz_21 = ppA * ( - spf * df_dx^2 + cpf * 2 * a)
    dg_dz_22 = ppA * cpf * df_dx * cpy

    dg_dz = [dg_dz_11  dg_dz_12; dg_dz_21 dg_dz_22] 

    ddt_dz_dz0 = dg_dz * dz_dz0
    dz[3] = ddt_dz_dz0[1, 1]
    dz[4] = ddt_dz_dz0[1, 2]
    dz[5] = ddt_dz_dz0[2, 1]
    dz[6] = ddt_dz_dz0[2, 2]
end

T = 20.0
t_span = (0.0,T) 

len = 150

FTLE = zeros((len, len))

x0_list = range(10, 11,length=len)
# y sicuramente non oltre 2
v0_list = range(0, 1,length=len)


i = 0
for x_0 in x0_list
    global i = i + 1
    j = 0
    for v_0 in v0_list
        j = j + 1
        z_0 = [x_0, v_0, 1, 0, 0, 1]
        prob = ODEProblem(doublegyre!,z_0,t_span,ϵ)
        sol = solve(prob)
        dzt_dz0_array = last(sol)[3:6]
        dzt_dz0 = [dzt_dz0_array[1] dzt_dz0_array[2]; dzt_dz0_array[3] dzt_dz0_array[4]]
        Delta = transpose(dzt_dz0)*dzt_dz0
        eig_val = eigen(Delta).values
        FTLE_ij = 1/T * log(sqrt(maximum(eig_val)))
        FTLE[i, j] = FTLE_ij 
    end
end

Plots.heatmap(x0_list, v0_list, transpose(FTLE))