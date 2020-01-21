###include("/home/jonathan/Documents/MPRO/RORT/Projet_RORT/lireInstance.jl") ###
#GROUPE 2 RORT MPRO 2020
using JuMP
filePath=string(@__DIR__,"/data/")
fileName="E_data.txt"
function readInstance(filePath,fileName)
    #readdlm("delim_file.txt", Int, '\t',Int,'\Int','\n')
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
    d=ones(Int8, n, n)*1000
    t=ones(Int8,n,n)*1000
    #for i in 1:n
    #    for j in 1:n
    #        if i~j
    #            d[i,j]=
    #        end
    #    end
    for i in 1:n
        d[i,i]=0
        t[i,i]=0
        d[arcs_array[i][1],arcs_array[i][2]]=arcs_array[i][3]
        t[arcs_array[i][1],arcs_array[i][2]]=arcs_array[i][4]
    end
    return n,nArcs,tauxConso,tauxRech,capa,sommets_array,d,t
end
print(readInstance(filePath,fileName))
