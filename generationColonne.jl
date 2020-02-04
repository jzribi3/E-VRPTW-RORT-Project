#groupe 2: J.Zribi, G.Dekeyser, P.Parmentier


# attention : l'indexation commence à 1 !! et ici notre n est égal au n+1 de l'article
# V_n+1 est entre 2 et n, V_0 désigne tous les sommets sauf le dépot d'arrivee donc de 1 à n-1
using JuMP
using GLPK

pm = Model(with_optimizer(GLPK.Optimizer))



n=4
nArc=6
r=1.0
g=1.0
Q=7
l0 =1000
sommets=[[1,0,1,0],[2,4,5,0],[3,4,5,0],[4,0,15,0]]
# indice, fenetre temps debut - fin - type sommet 0=client 1=station de recharge
s=[0,0,0,0]

d=[[0,4,2,1000],[1000,0,2,2],[1000,1,0,4],[1000,1000,1000,0]]
t=[[0,1,5,1000],[1000,0,4,5],[1000,3,0,1],[1000,1000,1000,0]]


include("heuristic.jl")

# filePath=string(@__DIR__,"/data/")
# fileName="E_data_3.txt"
# n,nArcs,tauxConso,tauxRech,capa,sommets_array,d,t = readInstance(filePath,fileName)
# V = []
# for i = 1 : length(sommets_array)
#     append!(V, [Vertex(sommets_array[i][1], sommets_array[i][2], sommets_array[i][3], sommets_array[i][4] == 1, 0)])
# end





function build_pm(dk,xk,n,nArc,r,g,Q,l0,sommets,d,t)
    K=length(dk)
    R = @variable(pm,[1:K])

    @objective(pm,Min,sum(dk[k]*R[k] for k in 1:K))

    global contraintePM
    contraintePM=Dict()

    for i in 2:(n-1)
        global contraintePM
        if sommets[i][4]==0
            contraintePM[i]=@constraint(pm, sum(sum(xk[k][i][j]*R[k] for k in 1:K) for j in 1:n)==1)
        end
    end
    println("model probleme maitre",pm)

    optimize!(pm)

    mu=Dict()

    mu[1]=0
    mu[n]=0

    for i in 2:(n-1)
        global contraintePM
        if sommets[i][4]==0
            print(JuMP.has_duals(pm))
            mu[i]=dual(contraintePM[i])
            print(mu[i])
        else
            mu[i]=0
        end
    end

    return mu, objective_value(pm)
end



function gene_col(mu,n,nArc,r,g,Q,l0,sommets,d,t)

    gc = Model(with_optimizer(GLPK.Optimizer))

    x = @variable(gc, [1:n,1:n], Bin)
    tau = @variable(gc,[1:n],Int)
    y = @variable(gc,[1:n],Int)
    selecSom = @variable(gc,[1:n],Bin)



    for i in 1:n
        @constraint(gc,sum(x[i,j] for j in 1:n)-selecSom[i]==0)
    end

    @constraint(gc,sum(x[1,j]-x[j,1] for j in 1:n)==1 )
    for i in 2:(n-1)
        @constraint(gc,sum(x[i,j]-x[j,i] for j in 1:n)==0)
    end
    @constraint(gc,sum(x[n,j]-x[j,n] for j in 1:n)==-1)

    @objective(gc,Min,sum(sum(d[i][j]*x[i,j] for i in 1:n) for j in 1:n) - sum(mu[i]*selecSom[i] for i in 1:n))



    tabClients=[]
    tabStations=[]


    for s in 2:(n-1)
        if sommets[s][4]<=0
            push!(tabClients,sommets[s][1])
        else
            push!(tabStations,sommets[s][1])
        end
    end



    nClients=length(tabClients)
    nStations=length(tabStations)

    #@constraint(gc,cons2[i in 1:nClients ],sum(x[tabClients[i],j] for j in 2:n if j!=tabClients[i])<=1)#2
    @constraint(gc,cons3[i in 1:nStations ],sum(x[tabStations[i],j] for j in 2:n if j!=tabStations[i])<=1)#3

    @constraint(gc,cons4[j in 2:(n-1)], sum(x[j,i] for i in 2:n if i!=j)-sum(x[i,j] for i in 1:(n-1) if i!=j)==0)#4

    #contrainte 5
    for i in 1:(n-1)
        if sommets[i][4]==0
            for j in 2:n
                if i!=j
                    println("i,j: ", i, ",", j)
                    println("x=1 ",tau[i]+(t[i][j]+s[i])*1-l0*(1-1)-tau[j])
                    println("x=0 ", tau[i]+(t[i][j]+s[i])*0-l0*(1-0)-tau[j])
                    @constraint(gc,tau[i]+(t[i][j]+s[i])*x[i,j]-l0*(1-x[i,j])-tau[j]<=0)
                end
            end
        end
    end

    #contrainte 6
    for i in 2:(n-1)
        if sommets[i][4]==1
            for j in 2:n
                if i!=j
                    @constraint(gc,tau[i]+(t[i][j])*x[i,j]-g*(Q-y[i])-(l0+g*Q)*(1-x[i,j])<=tau[j])
                end
            end
        end
    end

    #contrainte 7
    for j in 1:n
        @constraint(gc,sommets[j][2] <= tau[j] <= sommets[j][3])
    end

    #contrainte 8,9 : non nécessaire (concerne des livraisons)

    #contrainte 10

    for j in 2:n
        for i in 1:n
            if sommets[i][4]==0
                if i!=j
                    @constraint(gc,0<=y[j])
                    @constraint(gc,y[j]-y[i]+(r*d[i][j])*x[i,j]-Q*(1-x[i,j])<=0)
                end
            end
        end
    end

    #contrainte 11
    for j in 2:n
        for i in 1:n
            if sommets[i][4]==1 || i==1
                if i!=j
                    @constraint(gc,y[j]-Q+(r*d[i][j])*x[i,j]<=0)
                    @constraint(gc,0<=y[j])
                end
            end
        end
    end

    #println("model generation colonne",gc)

    optimize!(gc)

    d = sum(sum(d[i][j]*value(x[i,j]) for i in 1:n) for j in 1:n)



    return gc , d ,value.(x), objective_value(gc)
end

dk=[6.0,6.0]
xk=[[[0.,1.,0.,0.],[0.,0.,0.,1.],[0.,0.,0.,0.],[0.,0.,0.,0.]],[[0.,0.,1.,0.],[0.,0.,0.,0.],[0.,0.,0.,1.],[0.,0.,0.,0.]]]

mu , opt = build_pm(dk,xk,n,nArc,r,g,Q,l0,sommets,d,t)
println("mu",mu)
optimize!(pm)


gc,distroute ,xroute, coutReduit=gene_col(mu,n,nArc,r,g,Q,l0,sommets,d,t)

println("resultat",distroute, xroute )

compteur=0
coutRed = -1
while(compteur<10 && coutRed<0)
    global compteur
    global coutRed
    global distroute
    global xroute
    compteur = compteur +1
    push!(dk,distroute)
    xroute2=[ [0. for i in 1:n] for j in 1:n ]
    for i in 1:n
        for j in 1:n
            xroute2[i][j]=xroute[i,j]
        end
    end
    push!(xk,xroute2)
    mu, opt = build_pm(dk,xk,n,nArc,r,g,Q,l0,sommets,d,t)
    gc,distroute ,xroute, coutRed=gene_col(mu,n,nArc,r,g,Q,l0,sommets,d,t)
    print("cout reduit", coutRed)
    print("objectif", opt)
end


#print(termination_status(gc))
#
# println("nombre de camions :",sum(value(x[1,i]) for i in 1:n))
# println("distance :", objective_value(m))
# println("tau :", value.(tau))
# println("x ",value.(x))
# println(value(tau[3]+(t[3][2]+s[3])*value(x[3,2])-l0*(1-value(x[3,2]))-tau[2])<=0)
