using DifferentialEquations

using Plots;
using Statistics;
using Distributions;
using LinearAlgebra;

#use GR module
gr();

α = 1.0;
β = 1.0;
γ =-1.0;
ϵ = 0.025;

T = 30.0;
tspan = (0.0, T);

t_steps = 80
Δt = 1 / t_steps

# function Duffing(du, u, p, t)
#   du[1] = α * u[2]
#   du[2] = (β * u[1] + γ * u[1]^3)
#   du_du0 = [u[3] u[4]; u[5] u[6]]
#   dg_dz = [0 α; β+3*γ*u[1]^2 0]
#   ddt_dz_dz0 = dg_dz * du_du0
#   du[3] = ddt_dz_dz0[1, 1]
#   du[4] = ddt_dz_dz0[1, 2]
#   du[5] = ddt_dz_dz0[2, 1]
#   du[6] = ddt_dz_dz0[2, 2]
# end

# function σ_Duffing(du, u, p, t)
#   du[1] = 0.0
#   du[2] = ϵ * u[1]
#   du_du0 = [u[3] u[4]; u[5] u[6]]
#   db_dz = [0 0; ϵ 0]
#   ddt_dz_dz0 = db_dz * du_du0
#   du[3] = ddt_dz_dz0[1, 1]
#   du[4] = ddt_dz_dz0[1, 2]
#   du[5] = ddt_dz_dz0[2, 1]
#   du[6] = ddt_dz_dz0[2, 2]
# end

# len = 100
# FTLE = zeros((len, len))

# x0_list = range(-1.5, 1.5,length=len)
# v0_list = range(-1, 1,length=len)

# i = 0
# for x_0 in x0_list
#     global i = i + 1
#     println(i)
#     j = 0
#     Threads.@threads for v_0 in v0_list
#         j = j + 1
#         u0 = [x_0, v_0, 1, 0, 0, 1]
#         prob_sde_Duff = SDEProblem(Duffing, σ_Duffing, u0, tspan)
#         sol = solve(prob_sde_Duff, EM(), dt=Δt)
#         dzt_dz0_array = last(sol)[3:6]
#         dzt_dz0 = [dzt_dz0_array[1] dzt_dz0_array[2]; dzt_dz0_array[3] dzt_dz0_array[4]]
#         Delta = transpose(dzt_dz0)*dzt_dz0
#         eig_val = eigen(Delta).values
#         FTLE[j, i] = 1/T * log(sqrt(maximum(eig_val)))
#     end
# end

# Plots.heatmap(x0_list, v0_list, FTLE)
# Plots.savefig("matteo_juliacon/S_Duffing.png")

# # dx = a dt + b dW
# #d/dx0 dx = d/dx0 (a dt + b dW)
# #d dx/dx0 = da/dx dx/dx0 dt + db/dx dx/dx0 dW
# #....

function Duffing2d(du, u, p, t)
  du[1] = α * u[2]
  du[2] = (β * u[1] + γ * u[1]^3)
end

function σ_Duffing2d(du, u, p, t)
  du[1] = 0.0
  du[2] = ϵ * u[1]
end

nt = 200
times = range(0, T, length=nt)

Plots.plot()

x0 = 0.0; # x0_list[argmax(FTLE)[1]]
v0 = 0.0; # v0_list[argmax(FTLE)[2]]
# println(x0_list[argmax(FTLE)[1]])
# println(v0_list[argmax(FTLE)[2]])

mc_samples = 100 # 1000

fact = 0.1
x0_list = range(x0-fact, x0+fact,length=mc_samples)
v0_list = range(v0-fact, v0+fact,length=mc_samples)

sol_mc = zeros((2, nt, mc_samples^2))
i = 0 
for x_0 in x0_list, v_0 in v0_list
    global i = i + 1
    println(i/mc_samples^2)
    z_0 = [x_0, v_0]
    prob_sde_Duff = SDEProblem(Duffing2d, σ_Duffing2d, z_0, tspan)
    sol = solve(prob_sde_Duff, EM(), dt=Δt)
    sol_mc[1, :, i] = sol(times)[1, :]
    sol_mc[2, :, i] = sol(times)[2, :]
    # Plots.plot!(sol(times)[1, :], sol(times)[2, :], label="")
end

anim = @animate for l in range(1, nt, length=nt)
  println(l/nt)

  sol_t = sol_mc[:, Int(l), :]
  Plots.scatter(sol_t[1, :], sol_t[2, :], label="", xlims=(-1.5, 1.5), ylims=(-1, 1), markercolor=:red, markersize = 1)
  Plots.xlabel!("x")
  Plots.ylabel!("v")
end 
``
gif(anim, "matteo_juliacon/stoch_duff_mc.gif", fps = 10)

x0_min = 1.0 # x0_list[argmin(FTLE)[1]]
v0_min = 0.0 # v0_list[argmin(FTLE)[2]]

println(x0_min)
println(v0_min)

x0_min_list = range(x0_min-fact, x0_min+fact,length=mc_samples)
v0_min_list = range(v0_min-fact, v0_min+fact,length=mc_samples)

sol_mc2 = zeros((2, nt, mc_samples^2))
i = 0 
for x_0 in x0_min_list, v_0 in v0_min_list
    global i = i + 1
    println(i/mc_samples^2)
    z_0 = [x_0, v_0]
    prob_sde_Duff = SDEProblem(Duffing2d, σ_Duffing2d, z_0, tspan)
    sol = solve(prob_sde_Duff, EM(), dt=Δt)
    sol_mc2[1, :, i] = sol(times)[1, :]
    sol_mc2[2, :, i] = sol(times)[2, :]
end

anim = @animate for l in range(1, nt, length=nt)
  println(l/nt)

  sol_t2 = sol_mc2[:, Int(l), :]
  Plots.scatter(sol_t2[1, :], sol_t2[2, :], label="", xlims=(-1.5, 1.5), ylims=(-1, 1), markercolor=:blue, markersize = 1)
  Plots.xlabel!("x")
  Plots.ylabel!("v")
end 
``
gif(anim, "matteo_juliacon/stoch_duff_mc_min.gif", fps = 10)
