#=
Computational methods in power systems - formulation of Ybus
    1. Make a function that calculates Ybus directly.  The function must be general.
    2. Calculate Zbus, analyse the sparsity properties of both matrices plot(spy(abs.(y_bus))).
    3. Make a function that calculates the nodal powers and power flows. 
    The output must be a NamedTuple.
    4. Make a function that computes the node-branch incidence matrix.
    5. Calculate the Ybus from the incidence matrix.
    6. Make a function that performs Kron node elimination, from the step nodes.

=#

using Plots
using SparseArrays
using LinearAlgebra
using Base.Math
using Printf

# función para imprimir el diccionario
function imprimir(vec_d::Vector{NamedTuple})
    printstyled("Vector de Tuplas de tamaño ",length(vec_d),"\n";color=:green) 
    w = vec_d[1]
    nd = length(w)
    q = keys(w)
    for k = 1:nd
          print(q[k],'\t')
    end
    print("\n")
    for d in vec_d
          for k in q
                if typeof(d[k])==Float64
                      r = round(d[k],digits=4)
                else
                      r = d[k]
                end
                print("$r\t")
          end
          print("\n")
    end
end

# parametros del sistema de potencia tomado del libro de Anderson
num_lineas = 9
num_nodos  = 9
vec_lineas = Array{NamedTuple,1}(undef, num_lineas)
vec_lineas[1] = (n1=1, n2=4, r=0.0000, x=0.0576, b=0.0000) 
vec_lineas[2] = (n1=2, n2=7, r=0.0000, x=0.0625, b=0.0000)
vec_lineas[3] = (n1=3, n2=9, r=0.0000, x=0.0586, b=0.0000)
vec_lineas[4] = (n1=4, n2=5, r=0.0100, x=0.0850, b=0.0880)
vec_lineas[5] = (n1=4, n2=6, r=0.0170, x=0.0920, b=0.0790)
vec_lineas[6] = (n1=5, n2=7, r=0.0320, x=0.1610, b=0.1530)
vec_lineas[7] = (n1=6, n2=9, r=0.0390, x=0.1700, b=0.1790)
vec_lineas[8] = (n1=7, n2=8, r=0.0085, x=0.0720, b=0.0745)
vec_lineas[9] = (n1=9, n2=8, r=0.0119, x=0.1008, b=0.1045)

# Voltajes nodales (los ángulos están en grados)
vn = [1.0400; 1.0250; 1.0250; 1.0258; 0.9956; 1.0127; 1.0258; 1.0159; 1.0324]
an = [0; 9.2800; 4.6648; -2.2168; -3.9888; -3.6874; 3.7197; 0.7275; 1.9667]

#---------Construcción de la función Ybus_matrix
function Ybus_matrix(vec_lineas::Array{NamedTuple})
    num_nodos  = length(vec_lineas)
    Ybus = complex(zeros(num_nodos, num_nodos))
    for linea in vec_lineas
        y_lin = 1/(linea.r + linea.x*im)
        Ybus[linea.n1, linea.n1] = Ybus[linea.n1, linea.n1] + y_lin + linea.b*im
        Ybus[linea.n2, linea.n2] = Ybus[linea.n2, linea.n2] + y_lin + linea.b*im
        Ybus[linea.n2, linea.n1] = Ybus[linea.n2, linea.n1] - y_lin
        Ybus[linea.n1, linea.n2] = Ybus[linea.n1, linea.n2] - y_lin
    end
    return Ybus
end 
Ybus = Ybus_matrix(vec_lineas)
Ybus_sparse = sparse(Ybus);
Zbus = inv(Ybus);
#plot(spy(abs.(Ybus)))
#plot(spy(abs.(Zbus)))

#Función para calcular la tensión nodal en magnitud y ángulo radián
function calcular_Vbus(tension_nodal::Vector{Float64}, Theta::Vector{Float64})
    v_bus = tension_nodal.*(exp.(deg2rad.(Theta)*im))
    return v_bus
end
v_bus = calcular_Vbus(vn, an);

#Cálculo de la corriente 'I_bus'
I_bus = Ybus*v_bus;
#Cálculo de las potencias nodales
#S_k = v_bus.*conj(I_bus)

#Función para calcular las potencias nodales y el flujo de potencia por las lineas
function calcular_potencias(Ybus, tension_nodal, vec_lineas::Array{NamedTuple})
    #------Calcular las potencias nodales
    S = tension_nodal.*conj.(Ybus*tension_nodal)
    #P = real(S)
    #Q = imag(S)
    #------Calcular los flujos de potencia
    Ikm = complex(zeros(length(vec_lineas)))
    flujos = complex(zeros(length(vec_lineas)))
    for i in 1:length(vec_lineas)
        linea = vec_lineas[i]
        y_lin = 1/(linea.r + linea.x*im)
        Ikm[i] = (tension_nodal[linea.n1] - tension_nodal[linea.n2])*y_lin
        flujos = tension_nodal.*conj(Ikm)
    end
    # Retornar un NamedTuple con los resultados
    Results = (flujo_potencia = flujos, potencia_nodal = S)
    return Results
end
PyF = calcular_potencias(Ybus, v_bus, vec_lineas)

#Función para calcular la matriz de incidencia nodo-rama
function matriz_incidencia(vec_lineas, num_nodos)
    num_lineas = length(vec_lineas)
    incidencia = zeros(num_lineas, num_nodos)
    for i in 1:num_lineas
        n1 = vec_lineas[i].n1
        n2 = vec_lineas[i].n2
        incidencia[n1, i] = 1
        incidencia[n2, i] = -1
    end
    return incidencia
end 
A = matriz_incidencia(vec_lineas, num_nodos);
#Función para calcular la matriz de admitancias primitivas
function calculo_Yprim(vec_lineas::Array{NamedTuple})
    diag_Yprim = [ 1/(linea.r + linea.x*im) for linea in vec_lineas ]
    yprim = Diagonal(diag_Yprim)
    return yprim
end
Yprim = calculo_Yprim(vec_lineas);
#Matriz de admitancias Ybus a partir de la matriz de incidencia nodo-rama
Ybus_secondcalc = A * Yprim * A' 

#Eliminación de nodos de Kron, nodos de paso
# Y_kron = copy(Ybus);
#Obtener los nodos en los cuales la potencia activa es casi despreciable, es decir, no hay inyección de corriente
node = findall(abs.(real(PyF.potencia_nodal)) .< 0.005);
# funcion para calcular la matriz de reduccion de kron
function calc_Ykron(nodes, Ybus)
    Ykron = Ybus
    i_pos = 1
    node = nodes[i_pos]
    while node != 0
        # Crear un arreglo de indices sin el nodo
        ind_m = [1:node-1 ; node+1:size(Ykron)[1]]
        Y_NR = Ykron[node, ind_m] # Calculo de la Ynr
        Y_RN = Ykron[ind_m, node] # Calculo de la Yrn
        Y_RR = Ykron[node, node] # Calculo de la Yrr
        Y_NN = Ykron[ind_m, ind_m]
        Ykron = Y_NN - Y_NR*conj(Y_RN)' / Y_RR # Calculo de la Ykron
        i_pos += 1
        # Actualizar los indices de los nodos
        for i in i_pos:length(nodes)
            if(nodes[i] != 0)
            nodes[i] -= 1
        end
    end
    if(i_pos > length(nodes))
        break
    end
    # Actualizar el nodo
    node = nodes[i_pos]
    end
    return Ykron
end
Ykron_result = calc_Ykron(node, Ybus)

#Función para mostrar los resultados
function Resultados()
    printstyled("Potencias nodales:","\n";color=:green)
    for i in 1:length(PyF.potencia_nodal)
        imag_part = imag(PyF.potencia_nodal[i])
        sign = ifelse(signbit(imag_part), "-", "+")
        @printf("%.5f %s %.5fim\n", real(PyF.potencia_nodal[i]), sign, abs(imag_part))
    end
    printstyled("Flujos de potencia:","\n";color=:green)
    for i in 1:length(PyF.flujo_potencia)
        imag_part = imag(PyF.flujo_potencia[i])
        sign = ifelse(signbit(imag_part), "-", "+")
        @printf("%.5f %s %.5fim\n", real(PyF.flujo_potencia[i]), sign, abs(imag_part))
    end
    printstyled("Matriz Ybus sin los nodos de paso", "\n";color=:green)
    Ykron_result
end
Resultados()


# lines_page = @namedtuple(
#     bus1 = df_lines[!, :Bus1], 
#     bus2 = df_lines[!, :Bus2], 
#     length = df_lines[!, :Length], 
#     linecode = df_lines[!, :LineCode]
# );
# linecod_page = @namedtuple(
#     Name = df_linecod[!, :Name], 
#     R1 = df_linecod[!, :R1_Ohmkm], 
#     X1 = df_linecod[!, :X1_Ohmkm], 
#     R0 = df_linecod[!, :R0_Ohmkm], 
#     X0 = df_linecod[!, :X0_Ohmkm]
# );
# loads_page = @namedtuple(
#     Bus = df_loads[!, :Bus],
#     Phase = df_loads[!, :phase],
#     PF = df_loads[!, :PF],
#     Profile = df_loads[!, :Profile]
# );
# coord_page = @namedtuple(
#     Bus = df_coord[!, :Bus],
#     X_cord = df_coord[!, :X],
#     Y_cord = df_coord[!, :Y]
# );
# general_page = @namedtuple(
#     Vll = df_general[!, :Voltage_line_line_kV],
#     NomPower = df_general[!, :Nominal_Power_MW],
#     Vsubest = df_general[!, :Voltage_at_subestation_pu]
# );