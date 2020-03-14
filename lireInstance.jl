###include("/home/jonathan/Documents/MPRO/RORT/Projet_RORT/lireInstance.jl") ###
#GROUPE 2 RORT MPRO 2020
filePath=string(@__DIR__,"/data/")
fileName="E_data.txt"
function readInstance1(filePath,fileName)
    n=0
    nArcs=0
    tauxConso=0
    tauxRech=0
    capa=0
    line_array=[]
    sommets_array=[]
    arcs_array=[]
    open(string(filePath,fileName)) do file
        for (i, line) in enumerate(eachline(file))
            if i==1
                line_array=split(line)
                n=parse(Int,line_array[1])
                nArcs=parse(Int,line_array[2])
                tauxConso=parse(Float64,line_array[3])
                tauxRech=parse(Float64,line_array[4])
                capa=parse(Float64,line_array[5])

            elseif i<=n+1
                sommets_array_line=[parse(Int,x) for x in split(line)[1:4]]
                sommets_array=push!(sommets_array,sommets_array_line)

            elseif i<=n+nArcs+1
                arcs_array_line=[parse(Int,x) for x in split(line)[1:4]]
                arcs_array=push!(arcs_array,arcs_array_line)
            end

        end
    end
    d=[[1000 for k in 1:n] for j in 1:n]
    t=[[1000 for k in 1:n] for j in 1:n]
    for i in 1:n
        d[i][i]=0
        t[i][i]=0
    end
    for i in 1:nArcs
        d[arcs_array[i][1]][arcs_array[i][2]]=arcs_array[i][3]
        t[arcs_array[i][1]][arcs_array[i][2]]=arcs_array[i][4]
    end
    s=[0 for k in 1:n]
    return n,nArcs,tauxConso,tauxRech,capa,sommets_array,d,t,s
end

# read the files from the dataset of the article
# replaces all service times (variable s) by 0 since it is our simplifying hypothesis
function readInstance2(filePath,fileName)
    tauxConso=0
    tauxRech=0
    capa=0
    sommets_array=[]
    x = []
    y = []
    # start numbers for vertexes at 2. number 1 is only for the depot.
    index_next_vertex::Int32 = 2
    arrival_depot = []
    speed = 0
    open(string(filePath,fileName)) do file
        for (i, line) in enumerate(eachline(file))
            if i==1
                continue
            end
            line_array=split(line)
            println(line_array)
            if length(line_array) == 0
                continue;
            elseif line_array[2] == "d"
                x = append!([parse(Float64, line_array[3])], x)
                y = append!([parse(Float64, line_array[4])], y)
                depot = [1, parse(Float64, line_array[6]), parse(Float64, line_array[7]), 0]
                arrival_depot = copy(depot)
                sommets_array = append!([depot], sommets_array)
            elseif line_array[2] == "f"
                push!(x, parse(Float64, line_array[3]))
                push!(y, parse(Float64, line_array[4]))
                push!(sommets_array, [index_next_vertex, parse(Float64, line_array[6]), parse(Float64, line_array[7]), 1])
                index_next_vertex += 1
            elseif line_array[2] == "c"
                push!(x, parse(Float64, line_array[3]))
                push!(y, parse(Float64, line_array[4]))
                push!(sommets_array, [index_next_vertex, parse(Float64, line_array[6]), parse(Float64, line_array[7]), 0])
                index_next_vertex += 1
            elseif line_array[1] == "r"
                m = match(r"\/([0-9]*.[0-9]*)\/", line)
                tauxConso = parse(Float64, m[1])
            elseif line_array[1] == "g"
                m = match(r"\/([0-9]*.[0-9]*)\/", line)
                tauxRech = parse(Float64, m[1])
            elseif line_array[1] == "Q"
                m = match(r"\/([0-9]*.[0-9]*)\/", line)
                capa = parse(Float64, m[1])
            elseif line_array[1] == "v"
                m = match(r"\/([0-9]*.[0-9]*)\/", line)
                speed = parse(Float64, m[1])
            end
        end
    end
    arrival_depot[1] = index_next_vertex
    # add the arrival depot as last vertex
    push!(sommets_array, arrival_depot)
    push!(x, x[1])
    push!(y, y[1])
    # nb of vertexes, counting both depot vertexes
    n = index_next_vertex
    # directed graph is complete except for the arc depot - depot
    nArcs = n * n - 1
    println("sommets_array = $sommets_array")
    println("n = $n")
    println("tauxConso = $tauxConso")
    println("tauxRech = $tauxConso")
    println("capa = $capa")
    d::Array{Array{Float64, 1}, 1}=[[1000. for k in 1:n] for j in 1:n]
    t::Array{Array{Float64, 1}, 1}=[[1000. for k in 1:n] for j in 1:n]
    for i in 1:n
        d[i][i]=0
        t[i][i]=0
    end
    for i in 1:n-1
        for j in 2:n
            if i == 1 && j == n
                continue
            end
            d[i][j] = sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)
            t[i][j] = d[i][j] / speed
        end
    end
    s=[0 for k in 1:n]
    return n,nArcs,tauxConso,tauxRech,capa,sommets_array,d,t,s
end

# read first line of file in order to know which format it is, then call appropriate function
function readInstance(filePath, fileName)
    format::Int32= 0
    open(string(filePath, fileName)) do f
        line = readline(f)
        if line[1:8] == "StringID"
            format = 2
        else
            format = 1
        end
    end
    if format == 1
        return readInstance1(filePath, fileName)
    else
        return readInstance2(filePath, fileName)
    end
end

#print(readInstance(filePath,fileName))
#Ci-dessus: script pour demander le nom de fichier de données un utilisateur
#include("lireInstance.jl")
#filePath=string(@__DIR__,"/data/")
#print("Entrer le nom du fichier de données placé dans le dossier data (ex: E_data.txt): ")
#fileName= readline(stdin)
#println("Le fichier de données est :",fileName )
#println(readInstance(filePath,fileName))
readInstance("data/evrptw_instances/", "rc206_21.txt")
#readInstance("data/", "E_data_3.txt")