#=
Computational methods for power systems. Consider the 900-node CIGRE distribution system:
    1. Identify the maximum demand scenario.
    2. Calculate the load flow using Newton's method for this scenario.
    3. Show the convergence graph
    4. Compare with the fixed point method in terms of calculation time,
    number of iterations and accuracy.
=#

using LinearAlgebra
using Plots
using SparseArrays
using ForwardDiff
using DataFrames
using XLSX

#Cargar la base de datos en data frames
df_lines = DataFrame(XLSX.readtable("FEEDER900.xlsx", "lines")...);
df_linecod = DataFrame(XLSX.readtable("FEEDER900.xlsx", "line_codes")...);
df_coord = DataFrame(XLSX.readtable("FEEDER900.xlsx", "coordinates")...);
df_general = DataFrame(XLSX.readtable("FEEDER900.xlsx", "general")...);
df_loads = DataFrame(XLSX.readtable("FEEDER900.xlsx", "loads")...);
df_profile = DataFrame(XLSX.readtable("FEEDER900.xlsx", "profiles")...);

#Función para calcular Zbase
function calc_Zbase(general_sheet)
    Vbase = general_sheet[1, 1]
    Sbase = general_sheet[1, 2]
    Zbase = Vbase^2 / Sbase
    return Zbase
end
Zbase = calc_Zbase(df_general);

#Función que añade datos de cada hoja del archivo xlsx a un NamedTuple
function convert2NamedTuple(lines_sheet)
    num_rows = nrow(lines_sheet)
    tuple_lines = Array{NamedTuple, 1}(undef, num_rows)
    for rows in 1:num_rows
        Bus1 = lines_sheet[rows, 1]
        Bus2 = lines_sheet[rows, 2]
        #Es necesario pasar la longitud entre nodos a km
        Length = lines_sheet[rows, 3]/1000
        LineCodes = lines_sheet[rows, 4]
        #cambio a pu de R1 y X1 multiplicando por la dist en km y la Zbase
        r1 = (df_linecod[LineCodes, 2])*(Length)/Zbase
        x1 = (df_linecod[LineCodes, 3])*(Length)/Zbase
        coord_x1 = df_coord[Bus1, 2]
        coord_x2 = df_coord[Bus2, 2]
        coord_y1 = df_coord[Bus1, 3]
        coord_y2 = df_coord[Bus2, 3]
        #Se crea el NamedTuple con los datos
        tuple_lines[rows] = (bus1 = Bus1, bus2 = Bus2, length = Length, lincode = LineCodes, R1 = r1, X1 = x1, 
        Cx1 = coord_x1, Cx2 = coord_x2, Cy1 = coord_y1, Cy2 = coord_y2)

    end
    return tuple_lines
end
# [t.() for t in tuple_lines] vector columna con cada valor asignado en el NamedTuple
tuple_lines = convert2NamedTuple(df_lines);

#Función para calcular la matriz de admitancias
function matrix_Ybus(vec_lines)
    num_nodos = last(max([t.bus1 for t in tuple_lines], [t.bus2 for t in tuple_lines]))
    Ybus = complex(zeros(num_nodos, num_nodos))
    for line in vec_lines
        y_lin = 1/(line.R1 + line.X1*im)
        Ybus[line.bus1, line.bus1] = Ybus[line.bus1, line.bus1] + y_lin 
        Ybus[line.bus2, line.bus2] = Ybus[line.bus2, line.bus2] + y_lin 
        Ybus[line.bus2, line.bus1] = Ybus[line.bus2, line.bus1] - y_lin
        Ybus[line.bus1, line.bus2] = Ybus[line.bus1, line.bus2] - y_lin
    end
    return Ybus
end
#Ybus en por unidad y como matriz dispersa
Ybus = matrix_Ybus(tuple_lines);
#Ybus = sparse(Ybus);

function loadsNamedTuple(loads_sheet, time)
    num_rows = nrow(loads_sheet)
    tuple_loads = Array{NamedTuple, 1}(undef, num_rows)
    for rows in 1:num_rows
        BusL = loads_sheet[rows, 1] 
        Phase = loads_sheet[rows, 2]
        PF = loads_sheet[rows, 3]
        profile = loads_sheet[rows, 4]
        Pactv = df_profile[time, profile]/(df_general[1, 2]*1000)
        S = Pactv/PF
        Qreac = sqrt(S^2 - Pactv^2)
        tuple_loads[rows] = (busload = BusL, phase = Phase, pf = PF, P = Pactv, Q = Qreac)
    end
    return tuple_loads
end
tuple_loads = loadsNamedTuple(df_loads, 1);

#Función que cálcula el vector de potencias
function Power_vector(Loads_sheet)
    num_nodos = last(max([t.bus1 for t in tuple_lines], [t.bus2 for t in tuple_lines]))
    power = complex(zeros(num_nodos))
    for i in 1:num_nodos
        for loads in Loads_sheet
            if loads.busload == i
                power[loads.busload] = loads.P + loads.Q*im
            end
        end
    end
    return power
end

# Función para calcular el flujo de carga para operación nominal usando Método de Punto fijo
function end_point_method(Ybus, power)
    num_nodos = last(max([t.bus1 for t in tuple_lines], [t.bus2 for t in tuple_lines]))
    Yss = Ybus[1, 1]
    Yns = Ybus[1, 2:end]
    Ynn = Ybus[2:end, 2:end]
    Ynn = sparse(Ynn)
    Vs = 1.05+0*im
    Vs = [Vs]     
    # Método punto fijo
    Vn = complex(ones(num_nodos - 1))
    Sn = power[2:end]*-1
    n_iter = 10
    plot()
    Vn_real = complex(zeros(num_nodos - 1))
    Vn_old = complex(zeros(num_nodos - 1))
    err = complex(zeros(n_iter))
    Vbus = complex(zeros(num_nodos))
    for k in n_iter
        In = conj(Sn./Vn)
        Vn = Ynn \ (In - Yns .* Vs)
        Vbus = append!(Vs, Vn)
        Vn_real = real(Vn)
        err[k] = norm(Vn - Vn_old)
        Vn_old = Vn
    end
    return Vbus
end

#Función que calcula el flujo cuasi-dinámico 
function quasi_dynamic_flow(df_profile, df_loads, tuple_lines, Ybus)
    num_nodos = last(max([t.bus1 for t in tuple_lines], [t.bus2 for t in tuple_lines]))
    rowsprofile = nrow(df_profile)
    Vbus = complex(zeros(num_nodos))
    Ibus = complex(zeros(num_nodos))
    Sbus = complex(zeros(num_nodos))
    Slost = complex(zeros(rowsprofile))
    St = complex(zeros(rowsprofile))
    for rows in 1:rowsprofile
        V_loads = loadsNamedTuple(df_loads, rows)
        Pt = Power_vector(V_loads)
        Vbus = end_point_method(Ybus, Pt)
        Ibus = Ybus*Vbus
        Sbus = Vbus.*conj(Ibus)
        S_lost = Vbus'*Ibus
        Slost[rows] = S_lost
        St[rows] = Sbus[1,1]
    end
    dfVbus = DataFrame(Mag=abs.(Vbus[:,1]),Angle=angle.(Vbus[:,1]))
    return dfVbus, Slost, St, Vbus
end
df_Vbus, Slost, St, Vbus  = quasi_dynamic_flow(df_profile, df_loads, tuple_lines, Ybus);

display(plot(real(Slost), legend=false; title="System Losses", xlabel="Time (min)", ylabel="Power(p.u)"))
display(plot(real(St), legend=false; title="Total active power", xlabel="Time (min)", ylabel="Power (p.u)"))
display(plot(imag(St), legend=false; title="Total reactive power", xlabel="Time (min)", ylabel="power (p.u)"))
display(plot(abs.(Vbus), legend=false; title="Voltage", xlabel="Time (min)", ylabel="Voltage (p.u)"))

# Identificar el escenario de demanda máxima de un DF
function max_sum(df)
    max = 0
    index = 0
    for i in 1:size(df, 1)
        suma = sum(df[i, :])
        if suma > max
            max = suma
            index = i
        end
    end
    return max, index
end
_, index = max_sum(df_profile);

# Vector de cargas con la información del escenario de máxima demanda
num_cargas = size(df_loads, 1)
vec_loads = Array{NamedTuple, 1}(undef, num_cargas)
for i in 1:num_cargas
    local n = df_loads[i, 1]
    local phase = df_loads[i, 2] 
    local pf = df_loads[i, 3] 
    local P = df_profile[index, df_loads[i, 4]]/(df_general[1, 2]*1000) 
    local S = P/pf 
    local Q = sqrt(S^2-P^2) 
    local Scomp = P + Q*im 
    vec_loads[i] = (n = n, phase = phase, pf = pf, Scomp = Scomp)
end

# Función para calcular la potencia aparente, dado el escenario de máxima demanda
function Calc_Sn(vec_loads)
    num_nodos = last(max([t.bus1 for t in tuple_lines], [t.bus2 for t in tuple_lines]))
    S = complex(zeros(num_nodos - 1))
    for load in vec_loads
        S[load.n-1] = -load.Scomp
    end
    return S
end
Sn = Calc_Sn(vec_loads);

# Función para calcular la potencia activa, dado el escenario de máxima demanda
function Calc_Pn(vec_loads)
    num_nodos = last(max([t.bus1 for t in tuple_lines], [t.bus2 for t in tuple_lines]))
    P = zeros(num_nodos - 1)
    for load in vec_loads
        P[load.n-1] = -real(load.Scomp)
    end
    return P
end

# Función para calcular la potencia reactiva, dado el escenario de máxima demanda
function Calc_Qn(Vec_loads)
    num_nodos = last(max([t.bus1 for t in tuple_lines], [t.bus2 for t in tuple_lines]))
    Q = zeros(num_nodos - 1)
    for load in Vec_loads
        Q[load.n-1] = -imag(load.Scomp)
    end
    return Q
end

# Método de punto fijo
function end_point_M(Ybus, Sn, Vs, num_iter)
    num_nodos = last(max([t.bus1 for t in tuple_lines], [t.bus2 for t in tuple_lines]))
    Vn = complex(ones(num_nodos - 1))
    Yns = Ybus[1,2:end]
    Ynn = Ybus[2:end,2:end]
    Ynn = sparse(Ynn)
    conv = []
    for i in 1:num_iter
        In = conj(Sn./Vn)
        append!(conv, norm(Vn-Ynn\(In-Yns*Vs)))
        Vn = Ynn\(In-Yns*Vs)
        if conv[i] < 1e-6
            break
        end
    end
    Voltages = zeros(num_nodos, 2)
    Voltages[1,1] = abs(Vs)
    Voltages[1,2] = angle(Vs)
    for i in 1:num_nodos-1
        Voltages[i+1,1] = abs(Vn[i])
        Voltages[i+1,2] = angle(Vn[i])
    end
    return Voltages, Vn, conv
end

# Ecuaciones de flujo
function calc_F(Voltages, Ybus)
    P = zeros(size(Voltages)[1] - 1)
    Q = zeros(size(Voltages)[1] - 1)
    for k in 2:size(Voltages)[1]
        for m in 1:size(Voltages)[1]
            P[k - 1] +=
            real(Ybus[k, m])*Voltages[k, 1]*Voltages[m, 1]*cos(Voltages[k, 2]-Voltages[m, 2])+
            imag(Ybus[k, m])*Voltages[k, 1]*Voltages[m, 1]*sin(Voltages[k, 2]-Voltages[m, 2])
            Q[k - 1] +=
            real(Ybus[k, m])*Voltages[k, 1]*Voltages[m, 1]*sin(Voltages[k, 2]-Voltages[m, 2])-
            imag(Ybus[k, m])*Voltages[k, 1]*Voltages[m, 1]*cos(Voltages[k, 2]-Voltages[m, 2])
        end
    end
    return [P, Q]   
end

# Ecuación ΔF donde se realiza la diferencia entre la F y la Función objetivo (Max_demand)
function calc_deltaF(Fob, Voltages, Ybus)
    F0 = calc_F(Voltages, Ybus)
    return [Fob[1]-F0[1], Fob[2]-F0[2]]
end

function calc_jaco(Voltages, Ybus)
    P, Q = calc_F(Voltages, Ybus)
    Ybus = Ybus[2:end, 2:end]       # Separar el nodo slack de la Ybus
    Voltages = Voltages[2:end, :]   # Quitamos el nodo slack
    nNodos = size(Voltages)[1]
    Jptheta = zeros(nNodos, nNodos)
    Jpu = zeros(nNodos, nNodos)
    Jqtheta = zeros(nNodos, nNodos)
    Jqu = zeros(nNodos, nNodos)
    for k in 1:nNodos
        for m in 1:nNodos
            b = imag(Ybus[k, m])
            g = real(Ybus[k, m])
            theta = Voltages[k, 2] - Voltages[m, 2]
            if k == m
            # Cálculo de la diagonal principal de las submatrices del jacobiano
                Jptheta[k, m] = -b*Voltages[k, 1]^2 - Q[k]
                Jpu[k, m] = g*Voltages[k, 1] + P[k]/Voltages[k, 1]
                Jqtheta[k, m] = -g*Voltages[k, 1]^2 + P[k]
                Jqu[k, m] = -b*Voltages[k, 1] - Q[k]/Voltages[k, 1]
            else
            # Cálculo del resto de las submatrices del jacobiano
                Jpu[k, m] = g*Voltages[k, 1]*cos(theta) + b*Voltages[k, 1]*sin(theta)
                Jqu[k, m] = g*Voltages[k, 1]*sin(theta) - b*Voltages[k, 1]*cos(theta)
                Jptheta[k, m] = Jqu[k, m]*Voltages[m, 1]
                Jqtheta[k, m] = -Jpu[k, m]*Voltages[m, 1]
            end
        end
    end
    # Union de las submatrices del jacobiano
    Jacobiano = [Jptheta Jpu; Jqtheta Jqu]
    return Jacobiano
end

function Newton_Raphson(Fob, Voltages, Ybus)
    iSize = size(Voltages)[1]
    deltaF = calc_deltaF(Fob, Voltages, Ybus)
    error = []
    i = 1
    append!(error, norm(deltaF))
    while error[i] > 1e-6
        Jacobiano = calc_jaco(Voltages, Ybus)
        deltaF = [deltaF[1]; deltaF[2]]
        delta_volts = Jacobiano \ deltaF
        delta_volts1 = [0 0]
        delta_volts2 = [delta_volts[iSize:end] delta_volts[1:iSize-1]]
        delta_volts = [delta_volts1; delta_volts2]
        Voltages = delta_volts + Voltages
        deltaF = calc_deltaF(Fob, Voltages, Ybus)
        append!(error, norm(deltaF))
        i = i + 1
        if i > 10
            break
        end
    end
    return Voltages, error, i
end
# Tiempo de cálculo para método de Newton_Raphson
@time begin
    num_nodos = last(max([t.bus1 for t in tuple_lines], [t.bus2 for t in tuple_lines]))
    Pn = Calc_Pn(vec_loads) # Potencias activas de los nodos
    Qn = Calc_Qn(vec_loads) # Potencias reactivas de los nodos
    Fob = [Pn, Qn]
    Volts = ones(num_nodos - 1) # Voltajes iniciales
    theta = zeros(num_nodos - 1) # Angulos iniciales
    volts1 = [1.05 0]
    volts2 = [Volts theta]
    Voltages = [volts1; volts2]
    # Calcular Voltajes
    Volts, error, i = Newton_Raphson(Fob, Voltages, Ybus);
end
# Tiempo de cálculo para método de punto fijo
@time begin
    Sn = Calc_Sn(vec_loads) # Potencias de los nodos
    Vbus = complex(ones(1))
    Vbus[1] = 1.05
    # Calcular Voltajes
    Voltages, Vn, conv = end_point_M(Ybus, Sn, Vbus[1], 10)
end

# Gráfica de convergencia del método de punto fijo
display(plot(conv, title = "Convergencia del método de punto fijo",
xlabel = "Iteraciones", ylabel = "Error", label = "Error", legend = :topleft))
# Gráfica de convergencia del método de Newton-Raphson
display(plot(error, title = "Convergencia del método de Newton-Raphson",
xlabel = "Iteraciones", ylabel = "Error", label = "Error", legend = :topleft))