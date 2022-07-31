#= 
Author(s): Stefano Guidolotti, Matteo Manzi
=#

using Images
using Plots; gr()

N = 41 # n. of pixels

function my_f(l)
    im = Array(load("/home/matteo/Desktop/codes/StardustGroup/chaos/Arnold's Cat/recipes/Vladimir_Arnold.jpg"))

    x_s = 1:N;
    y_s = 1:N;
    y_0 = y_s' .* ones(N);
    x_0 = transpose(y_0);

    x = x_0
    y = y_0

    im1 = Array(load("/home/matteo/Desktop/codes/StardustGroup/chaos/Arnold's Cat/recipes/Vladimir_Arnold.jpg"))
    for k in 1:l
        x_old = x
        x = (2 * x + y) .% N
        y = (x_old + y) .% N;

        x.+=1
        y.+=1

        for i in x_s, j in y_s
            xcomp = Int(x[i, j])
            ycomp = Int(y[i, j])
            im[i, j] = im1[xcomp, ycomp]
        end

    end
    return im
end

n = 0:N/2

anim = @animate for l in n
    println(l)
    f = my_f(l)
    plot(f, axis=([], false))
end
``
gif(anim, "arnold.gif", fps = 3)