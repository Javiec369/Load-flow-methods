#=
Computational methods in power systems Consider the European low voltage system of 900 nodes.
    1. Load the system topology given in the "lines" table as a named_tuple.
    2. Convert the system to per unit (use only the positive sequence values).
    3. Formulate the Ybus matrix in per unit and as a sparse matrix.
    4. Show the scatter plot of the Ybus matrix.
    5. Plot the system using the positions given in the coordinates 
    table and the topology given in the lines table.
=#

using Plots
using SparseArrays
using LinearAlgebra
using DataFrames
using XLSX
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

df_lines = DataFrame(XLSX.readtable("FEEDER900.xlsx", "lines")...)
df_linecod = DataFrame(XLSX.readtable("FEEDER900.xlsx", "line_codes")...)
df_loads = DataFrame(XLSX.readtable("FEEDER900.xlsx", "loads")...);
df_coord = DataFrame(XLSX.readtable("FEEDER900.xlsx", "coordinates")...)
df_general = DataFrame(XLSX.readtable("FEEDER900.xlsx", "general")...)

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
        coord_x1 = (df_coord[Bus1, 2])
        coord_x2 = (df_coord[Bus2, 2])
        coord_y1 = (df_coord[Bus1, 3])
        coord_y2 = (df_coord[Bus2, 3])
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
Ybus = matrix_Ybus(tuple_lines)
Ybus=sparse(Ybus)
plot(spy(abs.(Ybus)))

# Función para gráficar la red
function plotting_network(vec_lines)
    # Ubicación en X1 y X2 para cada par de nodos
    X = []
    # Ubicación en Y1 e Y2 para cada par de nodos
    Y = []
    #'plot' traza una linea de ancho 'linewidth=' entre X e Y
    plot(X, Y, linewidth=1, label=false)
    for line in vec_lines
        X = [line.Cx1; line.Cx2]
        Y = [line.Cy1; line.Cy2]
        # 'plot!' traza una linea de ancho 'linewidth=' entre X e Y sobre el plot ya existente
        plot!(X, Y, linewidth=2, label=false, linecolor="blue")
    end
    x = df_coord[:,2]
    y = df_coord[:,3]
    # 'scatter!' gráfica un punto dadas sus coordenadas x e y en un gráfico ya existente 
    display(scatter!(x, y, markercolor="red", markersize=2.5))
end
grid = plotting_network(tuple_lines)