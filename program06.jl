#=
Computational methods for power systems.
    1. Take the 900-node IEEE system.  Assume balanced three-phase operation. 
    Calculate the Ybus in per unit and model the loads as constant power.
    2. Calculate load flow for nominal operation using the fixed point method. 
    Show a graph of the magnitude of the nodal voltages. Show voltages and angles 
    in tables (data-frame).
    3. Calculate the quasi-dynamic load flow. Show a graph of the system losses 
    and the active and reactive power in the substation.
=#

using Plots
using Base.Math
using SparseArrays
using LinearAlgebra
using DataFrames
using XLSX

# Cargar la base de datos en un  DataFrame
df_lines = DataFrame(XLSX.readtable("FEEDER900.xlsx", "lines")...);
df_linecod = DataFrame(XLSX.readtable("FEEDER900.xlsx", "line_codes")...);
df_coord = DataFrame(XLSX.readtable("FEEDER900.xlsx", "coordinates")...);
df_general = DataFrame(XLSX.readtable("FEEDER900.xlsx", "general")...);
df_loads = DataFrame(XLSX.readtable("FEEDER900.xlsx", "loads")...);
df_profile = DataFrame(XLSX.readtable("FEEDER900.xlsx", "profiles")...);

# Función para calcular Zbase
function calc_Zbase(general_sheet)
    Vbase = general_sheet[1, 1]
    Sbase = general_sheet[1, 2]
    Zbase = Vbase^2 / Sbase
    return Zbase
end
Zbase = calc_Zbase(df_general);

# Función que añade datos de cada hoja del archivo xlsx a un NamedTuple
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

# Función para calcular la matriz de admitancias
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
power = Power_vector(tuple_loads);

# Función para calcular el flujo de carga para operación nominal usando Método de Punto fijo
function end_point_method(Ybus, power)
    num_nodos = last(max([t.bus1 for t in tuple_lines], [t.bus2 for t in tuple_lines]))
    Yss = Ybus[1, 1]
    Yns = Ybus[1, 2:end]
    Ynn = Ybus[2:end, 2:end]
    Ynn = sparse(Ynn)
    Vs = 1.05+1.0*im
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
Vbus = end_point_method(Ybus, power);

#Función que calcula el flujo cuasi-dinámico 
function quasi_dynamic_flow(df_profile, df_loads, tuple_lines, Ybus)
    num_nodos = last(max([t.bus1 for t in tuple_lines], [t.bus2 for t in tuple_lines]))
    Vs = 1.05+0*im
    Vs = [Vs]
    rowsprofile = nrow(df_profile)
    Vbus = complex(zeros(num_nodos))
    Ibus = complex(zeros(num_nodos))
    Sbus = complex(zeros(num_nodos))
    Slost = complex(zeros(rowsprofile))
    St = complex(zeros(rowsprofile))
    plot()
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

display(plot!(real(Slost),legend=false;title="System Losses",xlabel="Time (min)",ylabel="Power(p.u)"))
display(plot(real(St),legend=false;title="Total active power",xlabel="Time (min)",ylabel="Power (p.u)"))
display(plot(imag(St),legend=false;title="Total reactive power",xlabel="Time (min)",ylabel="power (p.u)"))
display(plot(abs.(Vbus), legend=false; title="Voltage", xlabel="Time (min)", ylabel="Voltage (p.u)"))