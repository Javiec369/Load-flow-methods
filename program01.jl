#=
    End point method
=#
using LinearAlgebra

function End_point_method(x_ini::Number)
    T(x) = cos(x)
    n_iter = 30
    x = zeros(n_iter, 1)
    err = zeros(n_iter)
    x[1] = x_ini
    err[1] = Inf
    for k in 2:n_iter
        x[k] = T(x[k-1])
        err[k] = norm(x[k]-x[k-1])
    end
    return x, err
end 

x_sol, errd = Met_punto_fijo(0.0)

#Plot err using log scale
using Plots
errd
x_sol
plot(errd)
plot(errd, yaxis=:log, title="Log scale graph")
