#= 
Author(s): Matteo Manzi
=#

using DifferentialEquations
using LinearAlgebra
using Plots; gr()
using PlotlyJS
using GLMakie


function dg_dz(z, ϵ)
    dg_dz = [1 ϵ*cos(z[2]); 1 1 + ϵ*cos(z[2])]
    return dg_dz
end 

ϵ = 10

T = 100
t_span = 0:T

len = 1000

FTLE = zeros((len, len, length(t_span)))

y0_list = range(-3, 3,length=len)
x0_list = range(0, 2π,length=len)
i = 0
for y0 in y0_list
    global i = i + 1
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
            dz_dz0 = dg_dz(z, ϵ) * dz_dz0

            Delta = transpose(dz_dz0)*dz_dz0
            eig_val = eigen(Delta).values
            FTLE_ijk = 1/T * log(sqrt(maximum(eig_val)))
            FTLE[i, j, k+1] = FTLE_ijk
        end
        
    end
end



anim = @animate for l in t_span
    print(l/T)
    print(' ')

    FTLE_t = FTLE[:, :, l+1]

    Plots.heatmap(x0_list, y0_list, FTLE_t)
end 
``
gif(anim, "standardmap_time.gif", fps = 10)