#=
    Fixed point method
=#
using LinearAlgebra
function punto_fijo(A::Matrix, b::Vector)
    x = [0;0]
    n_iter = 30
    err = zeros(n_iter)
    An = I - A
    L = norm(An)
    if L > 1
        println("No cumple con el criterio de convergencia")
    else
        for k = 1:n_iter
            x = b + An*x
            err[k] = norm(A*x-b)
        end
    return x,err
    end
end 

A = [1 0.5;0.3 1]
b = [3; 5]

x, err = punto_fijo(A, b)
using Plots
plot(err, yaxis=:log, title="Gráfico en escala logarítmica", 
    xlabel="Eje x", ylabel="Eje y log10")