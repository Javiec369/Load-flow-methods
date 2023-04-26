#= 
    Exercise fixed point numerical method
=#
using LinearAlgebra
using Plots

function punto_fijo(v1, P, Q, x_eq)
    #Datos para el método
    n_iter = 40
    err = zeros(n_iter)
    #Valor inicial de V2 para iniciar con las iteraciones
    v2 = 0.5 + 0.5im
    #Con 'v2_iter' se crea una lista vacía de tipo compleja con valores flotantes
    v2_iter = Complex{Float64}[]
    for k in 1:n_iter
        #Ecuación de la tensión V2 determinada a partir de un previo análisis
        v2 = v1 - (x_eq*im) * (conj(P + Q*im)  / conj(v2))
        #Corriente en la carga
        I  = (v1 - v2) / (x_eq*im) 
        #Potencia compleja en la carga
        S  = (v2 * conj(I)) 
        err[k] = norm(P + Q*im - S)
        #'push!' agregar elementos a un vector o matriz. Sintaxis: push!(vector, elemento)
        push!(v2_iter, v2) 
    end
    return v2, err, v2_iter
end

#Datos del cto
v1   = 1        #Tensión en el primero nodo
P    = 10.0     #Potencia activa de la carga (pu)
Q    = 3.0      #Potencia reactiva de la carga (pu)
xt1  = 0.008    #Reactancia del primer transformador (pu)
xt2  = 0.009    #Reactancia del segundo transformador (pu)
xl   = 0.001    #Reactancia de la linea de transmisión (pu)
x_eq = xt1 + xt2 + xl #Reactancia equivalente

v2_sol, errd, v2 = punto_fijo(v1, P, Q, x_eq)

#Podemos llamar 'v2' puesto que se almacenó sus iteraciones en un vector.
v2
println("Solucion V2: ",v2_sol,"\nMagnitud: ",abs(v2_sol),
        ", Ángulo (grados): ",rad2deg(angle(v2_sol)))

plot(errd, title="Grafica del error", xlabel="Eje x", ylabel="Eje y")
plot(errd[1:25], yaxis=:log, title="Grafica del error en escala logarítmica", 
    xlabel="Eje x", ylabel="Eje y (log10)")
