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

α = 2.5 
T = 10.0 
t_span = (0.0,T)


# function pendulum!(dz,z,α,t)
#     dz[1] = z[2]
#     dz[2] = (α*cos(5t)-1)*sin(z[1])
#     dz_dz0 = [z[3] z[4]; z[5] z[6]]
#     dg_dz = [0 1; cos(z[1])*(α*cos(5t) - 1) 0]
#     ddt_dz_dz0 = dg_dz * dz_dz0
#     dz[3] = ddt_dz_dz0[1, 1]
#     dz[4] = ddt_dz_dz0[1, 2]
#     dz[5] = ddt_dz_dz0[2, 1]
#     dz[6] = ddt_dz_dz0[2, 2]
# end

# len = 50
# FTLE = zeros((len, len))

# ran = pi
# x0_list = range(-ran, ran,length=len)
# v0_list = range(-ran, ran,length=len)

# i = 0
# for x_0 in x0_list
#     global i = i + 1
#     j = 0
#     for v_0 in v0_list
#         j = j + 1
#         z_0 = [x_0, v_0, 1, 0, 0, 1]
#         prob = ODEProblem(pendulum!,z_0,t_span,α)
#         sol = solve(prob)
#         dzt_dz0_array = last(sol)[3:6]
#         dzt_dz0 = [dzt_dz0_array[1] dzt_dz0_array[2]; dzt_dz0_array[3] dzt_dz0_array[4]]
#         Delta = transpose(dzt_dz0)*dzt_dz0
#         eig_val = eigen(Delta).values
#         FTLE[j, i] = 1/T * log(sqrt(maximum(eig_val)))
#     end
# end

# Plots.heatmap(x0_list, v0_list, FTLE)



function pendulum2d!(dz,z,α,t)
    dz[1] = z[2]
    dz[2] = (α*cos(5t)-1)*sin(z[1])
end

nt = 200
times = range(0, T, length=nt)

Plots.plot()
x0 = -2.885136110439606
v0 = -2.3722230241392315

mc_samples = 1000

fact = 0.1
x0_list = range(x0-fact*abs(x0), x0+fact*abs(x0),length=mc_samples)
v0_list = range(v0-fact*abs(v0), v0+fact*abs(v0),length=mc_samples)

sol_mc = zeros((2, nt, mc_samples^2))
i = 0 
for x_0 in x0_list, v_0 in v0_list
    global i = i + 1
    # println(i)
    z_0 = [x_0, v_0]
    prob = ODEProblem(pendulum2d!,z_0,t_span,α)
    sol = solve(prob)
    # print(sol(times)[1, :])
    sol_mc[1, :, i] = sol(times)[1, :]
    sol_mc[2, :, i] = sol(times)[2, :]
    # Plots.plot!(sol(times)[1, :], sol(times)[2, :], label="")
end
# Plots.savefig("Traj.png")

x0 = -0.320570678937734 # println(x0_list[argmin(FTLE)[1]])
v0 = -0.320570678937734 # println(v0_list[argmin(FTLE)[2]])

x0_list = range(x0-fact*abs(x0), x0+fact*abs(x0),length=mc_samples)
v0_list = range(v0-fact*abs(v0), v0+fact*abs(v0),length=mc_samples)

sol_mc2 = zeros((2, nt, mc_samples^2))
i = 0 
for x_0 in x0_list, v_0 in v0_list
    global i = i + 1
    # println(i)
    z_0 = [x_0, v_0]
    prob = ODEProblem(pendulum2d!,z_0,t_span,α)
    sol = solve(prob)
    # print(sol(times)[1, :])
    sol_mc2[1, :, i] = sol(times)[1, :]
    sol_mc2[2, :, i] = sol(times)[2, :]
end

anim = @animate for l in range(1, nt, length=nt)
    println(l/nt)

    sol_t = sol_mc[:, Int(l), :]
    Plots.scatter(sol_t[1, :], sol_t[2, :], label="", xlims=(-40, -2), ylims=(-5, -1), markercolor=:red, markersize = 1)
    Plots.xlabel!("x")
    Plots.ylabel!("v")
end 
``
gif(anim, "matteo_juliacon/ppp_mc_max.gif", fps = 10)

anim = @animate for l in range(1, nt, length=nt)
    println(l/nt)

    sol_t2 = sol_mc2[:, Int(l), :]
    Plots.scatter(sol_t2[1, :], sol_t2[2, :], label="", xlims=(-0.6, 0.6), ylims=(-0.6, 0.6), markercolor=:blue, markersize = 1)
    Plots.xlabel!("x")
    Plots.ylabel!("v")
end 
``
gif(anim, "matteo_juliacon/ppp_mc_min.gif", fps = 10)
