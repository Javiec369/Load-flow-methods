#=
Jacobian and newton Rapshon's method
    1. Find the values of theta_2 and theta_3 for the system of equations 
    shown using Newton's method.
    2. Show a convergence graph on a logarithmic scale.
    3. Solve the problem with the ForwardDiff differentiation package 
    and directly formulating the derivative. Compare computation times.  

    Given a function f and variables x, the Jacobian can be calculated as follows:

        jac = x -> ForwardDiff.jacobian(f,x)

    this function is executed as

        Jx = jac(x) 

    4. Display a heatmap graph where the x-axis is the initial condition for theta_2, 
    the y-axis is the initial condition for theta_3 and the z-axis is the number of iterations.
=#

using LinearAlgebra
using SparseArrays
using Plots
using ForwardDiff

# Definir las funciones
function F(theta)
    theta2, theta3 = theta
    h1 = 0.4923 - 0.3978*cos(theta2) - 1.888*sin(theta2) - 0.2489*cos(theta2-theta3) - 2.0154*sin(theta2-theta3)
    h2 = 3.9234 - 1.6390*cos(theta3) - 7.269*sin(theta3) - 0.2439*cos(theta2-theta3) + 2.0154*sin(theta2-theta3)
    return [h1, h2]
end

# Definir la matriz jacobiana
function Djac_F(theta)
    theta2, theta3 = theta
    J11 = 0.3978*sin(theta2) - 1.888*cos(theta2) + 0.2489*sin(theta2-theta3) - 2.0154*cos(theta2-theta3)
    J12 = -0.2489*sin(theta2-theta3) + 2.0154*cos(theta2-theta3)
    J21 = 0.2439*sin(theta2-theta3) + 2.0154*cos(theta2-theta3)
    J22 = 1.6390*sin(theta3) - 7.269*cos(theta3) - 0.2439*sin(theta2-theta3) - 2.0154*cos(theta2-theta3)
    return [J11 J12; J21 J22]
end

# Calcular el Jacobiano a partir de la librería ForwardDiff
function jac_fowardDiff(θ)
    return ForwardDiff.jacobian(F, θ)
end

# Definir método de Newton-Raphson
function newton_raphson(F, DF, x0, tol=1e-6, maxiter=100)
    x = copy(x0)
    res = norm(F(x))
    iter = 0
    res_history = [res]
    while res > tol && iter < maxiter
        Δx = - DF(x) \ F(x)
        x += Δx
        res = norm(F(x))
        iter += 1
        push!(res_history, res)
    end
    return x, iter, res_history
end
# Encontrar los valores de theta_2 y theta_3 y graficar la convergencia 
@time begin
theta, iter, res_history = newton_raphson(F, jac_fowardDiff, [0.0, 0.0])
end
# plot(res_history, yaxis=:log10);

# Crear matriz de datos que almacena las iteraciones según X0
n = 300;
data = zeros(n, n);
theta2_vals = range(-0.8, stop=0.8, length=n);
theta3_vals = range(-0.8, stop=0.8, length=n);
for i in 1:n
    for j in 1:n
        x0 = [theta2_vals[i], theta3_vals[j]]
        _, iter, _ = newton_raphson(F, jac_fowardDiff, x0)
        data[i,j] = iter
    end
end
# Generar gráfico de tipo heatmap
heatmap(theta2_vals, theta3_vals, data, 
    xlabel="Condición inicial para theta_2", 
    ylabel="Condición inicial para theta_3", 
    zlabel="Número de iteraciones", 
    title="Convergencia del método de Newton-Raphson")