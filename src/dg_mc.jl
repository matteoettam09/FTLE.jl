#= 
Author(s): Matteo Manzi
=#

using DifferentialEquations
using LinearAlgebra
using Plots; gr()
using GLMakie
using RecipesBase

A = 0.1
ω = 2π/10
ϵ = 0.1

piA = π*A
ppA = π^2*A

T = 50.0
t_span = (0.0,T) 

# function doublegyre!(dz,z,ϵ,t)
#     a = ϵ * sin(ω*t)
#     b = 1 - 2*a
#     f = a*z[1]^2 + b*z[1]
#     spf = sin(π*f)
#     cpf = cos(π*f)
#     spy = sin(π*z[2])
#     cpy = cos(π*z[2])

#     df_dx = (2*a*z[1] + b)

#     dz[1] = - piA*spf*cpy
#     dz[2] = piA*cpf*spy*df_dx

#     dz_dz0 = [z[3] z[4]; z[5] z[6]]

#     dg_dz_11 = - ppA * cpy * cpf * df_dx
#     dg_dz_12 = ppA* spf * spy
#     dg_dz_21 = ppA * ( - spf * df_dx^2 + cpf * 2 * a)
#     dg_dz_22 = ppA * cpf * df_dx * cpy

#     dg_dz = [dg_dz_11  dg_dz_12; dg_dz_21 dg_dz_22] 

#     ddt_dz_dz0 = dg_dz * dz_dz0
#     dz[3] = ddt_dz_dz0[1, 1]
#     dz[4] = ddt_dz_dz0[1, 2]
#     dz[5] = ddt_dz_dz0[2, 1]
#     dz[6] = ddt_dz_dz0[2, 2]
# end

# len = 100

# FTLE = zeros((len, len))

# x0_list = range(0, 2,length=len)
# # y sicuramente non oltre 2
# v0_list = range(0, 1,length=len)


# i = 0
# for x_0 in x0_list
#     global i = i + 1
#     j = 0
#     for v_0 in v0_list
#         j = j + 1
#         z_0 = [x_0, v_0, 1, 0, 0, 1]
#         prob = ODEProblem(doublegyre!,z_0,t_span,ϵ)
#         sol = solve(prob)
#         dzt_dz0_array = last(sol)[3:6]
#         dzt_dz0 = [dzt_dz0_array[1] dzt_dz0_array[2]; dzt_dz0_array[3] dzt_dz0_array[4]]
#         Delta = transpose(dzt_dz0)*dzt_dz0
#         eig_val = eigen(Delta).values
#         FTLE_ij = 1/T * log(sqrt(maximum(eig_val)))
#         FTLE[j, i] = FTLE_ij 
#     end
# end

# Plots.heatmap(x0_list, v0_list, FTLE)

# println(x0_list[argmax(FTLE)[1]])
# println(v0_list[argmax(FTLE)[2]])

function dg2d!(dz,z,ϵ,t)
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
end

nt = 200
times = range(0, T, length=nt)

Plots.plot()
x0 = 2
v0 = 0.5353535353535354

mc_samples = 1000

fact = 0.01
x0_list = range(x0-fact, x0+fact,length=mc_samples)
v0_list = range(v0-fact, v0+fact,length=mc_samples)

sol_mc = zeros((2, nt, mc_samples^2))
i = 0 
for x_0 in x0_list, v_0 in v0_list
    global i = i + 1
    println(i/1000000)
    z_0 = [x_0, v_0]
    prob = ODEProblem(dg2d!,z_0,t_span,ϵ)
    sol = solve(prob)
    # print(sol(times)[1, :])
    sol_mc[1, :, i] = sol(times)[1, :]
    sol_mc[2, :, i] = sol(times)[2, :]
    # Plots.plot!(sol(times)[1, :], sol(times)[2, :], label="")
end
# Plots.savefig("Traj.png")

x0 = 0.48484848484848486 # println(x0_list[argmin(FTLE)[1]])
v0 = 0.5959595959595959 # println(v0_list[argmin(FTLE)[2]])

x0_list = range(x0-fact, x0+fact,length=mc_samples)
v0_list = range(v0-fact, v0+fact,length=mc_samples)

sol_mc2 = zeros((2, nt, mc_samples^2))
i = 0 
for x_0 in x0_list, v_0 in v0_list
    global i = i + 1
    # println(i)
    z_0 = [x_0, v_0]
    prob = ODEProblem(dg2d!,z_0,t_span,ϵ)
    sol = solve(prob)
    # print(sol(times)[1, :])
    sol_mc2[1, :, i] = sol(times)[1, :]
    sol_mc2[2, :, i] = sol(times)[2, :]
end

anim = @animate for l in range(1, nt, length=nt)
    println(l/nt)

    sol_t = sol_mc[:, Int(l), :]
    Plots.scatter(sol_t[1, :], sol_t[2, :], label="", xlims=(0, 4), ylims=(0, 1), markercolor=:red, markersize = 1)
    Plots.xlabel!("x")
    Plots.ylabel!("v")
end 
``
gif(anim, "dg_mc_max.gif", fps = 10)

anim = @animate for l in range(1, nt, length=nt)
    println(l/nt)

    sol_t2 = sol_mc2[:, Int(l), :]
    Plots.scatter(sol_t2[1, :], sol_t2[2, :], label="",  xlims=(0.2, 0.8), ylims=(0.2, 0.8), markercolor=:blue, markersize = 1)
    Plots.xlabel!("x")
    Plots.ylabel!("v")
end 
``
gif(anim, "dg_mc_min.gif", fps = 10)
