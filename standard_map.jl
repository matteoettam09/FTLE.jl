#= 
Author(s): Matteo Manzi
=#

using DifferentialEquations
using LinearAlgebra
using Plots; gr()
# using Colors

using PlotlyJS
using GLMakie
# using CairoMakie

using Distributed

function dg_dz(z, ϵ)
    dg_dz = [1 ϵ*cos(z[2]); 1 1 + ϵ*cos(z[2])]
    return dg_dz
end 


T = 5
t_span = 0:T

len = 300

FTLE = zeros((len, len))

function my_f(ϵ)
    global y0_list = range(-π, π,length=len)  # -3, 3
    global x0_list = range(0, 2π,length=len)
    i = 0
    @distributed for y0 in y0_list
        i = i + 1
        j = 0
        for x0 in x0_list

            j = j + 1
            
            dz_dz0 = [1 0; 0 1]
            yk = y0
            xk = x0
            for k in t_span
                x_old = xk
                xk = mod(xk + yk + ϵ * sin(xk), 2π)
                yk = yk + ϵ * sin(x_old)
                z = [yk xk]
                # print(z)
                dz_dz0 = dg_dz(z, ϵ) * dz_dz0
            end
            Delta = transpose(dz_dz0)*dz_dz0
            eig_val = eigen(Delta).values
            FTLE_ij = 1/T * log(sqrt(maximum(eig_val)))
            # print(FTLE_ij)

            FTLE[i, j] = FTLE_ij 
        end
    end
    return FTLE
end


# f = my_f(0.7)
# # print(f)
# Plots.heatmap(y0_list, x0_list, f, axis=([], false), legend =:none)

n = 10
N = 5
gif_range = range(0.01, stop = N, length = n)
gif_range_contrario = range(N, stop = 0.01, length = n)
gra = vcat(gif_range, gif_range_contrario) 

anim = @animate for l in gra
    print(l/N)
    print(' ')

    f = my_f(l)

    Plots.heatmap(y0_list, x0_list, f, axis=([], false), legend =:none)
end
``
gif(anim, "standard_map_sito.gif", fps = 10)