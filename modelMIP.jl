#groupe 2: J.Zribi, G.Dekeyser, P.Parmentier


# attention : l'indexation commence à 1 !! et ici notre n est égal au n+1 de l'article
# V_n+1 est entre 2 et n, V_0 désigne tous les sommets sauf le dépot d'arrivee donc de 1 à n-1
using JuMP
using GLPK

m = Model(with_optimizer(GLPK.Optimizer))

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


@variable(m, x[1:n,1:n], Bin)
@variable(m,tau[1:n],Int)
@variable(m,u[1:n],Int)
@variable(m,y[1:n],Int)





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


@constraint(m,cons2[i in 1:nClients ],sum(x[tabClients[i],j] for j in 2:n if j!=tabClients[i])==1)#2
@constraint(m,cons3[i in 1:nStations ],sum(x[tabStations[i],j] for j in 2:n if j!=tabStations[i])<=1)#3

@constraint(m,cons4[j in 2:(n-1)], sum(x[j,i] for i in 2:n if i!=j)-sum(x[i,j] for i in 1:(n-1) if i!=j)==0)#4

#contrainte 5
for i in 1:(n-1)
    if sommets[i][4]==0
        for j in 2:n
            if i!=j
                println("i,j: ", i, ",", j)
                println("x=1 ",tau[i]+(t[i][j]+s[i])*1-l0*(1-1)-tau[j])
                println("x=0 ", tau[i]+(t[i][j]+s[i])*0-l0*(1-0)-tau[j])
                @constraint(m,tau[i]+(t[i][j]+s[i])*x[i,j]-l0*(1-x[i,j])-tau[j]<=0)
            end
        end
    end
end

#contrainte 6
for i in 2:(n-1)
    if sommets[i][4]==1
        for j in 2:n
            if i!=j
                @constraint(m,tau[i]+(t[i][j])*x[i,j]-g*(Q-y[i])-(l0+g*Q)*(1-x[i,j])<=tau[j])
            end
        end
    end
end

#contrainte 7
for j in 1:n
    @constraint(m,sommets[j][2] <= tau[j] <= sommets[j][3])
end

#contrainte 8,9 : non nécessaire (concerne des livraisons)

#contrainte 10

for j in 2:n
    for i in 1:n
        if sommets[i][4]==0
            if i!=j
                @constraint(m,0<=y[j])
                @constraint(m,y[j]-y[i]+(r*d[i][j])*x[i,j]-Q*(1-x[i,j])<=0)
            end
        end
    end
end

#contrainte 11
for j in 2:n
    for i in 1:n
        if sommets[i][4]==1 || i==1
            if i!=j
                @constraint(m,y[j]-Q+(r*d[i][j])*x[i,j]<=0)
                @constraint(m,0<=y[j])
            end
        end
    end
end

println("model",m)
@objective(m, Min, sum(sum(d[i][j]*x[i,j] for i in 1:(n-1) if i!=j) for j in 2:n))

optimize!(m)
print(termination_status(m))

println("nombre de camions :",sum(value(x[1,i]) for i in 1:n))
println("distance :", objective_value(m))
println("tau :", value.(tau))
println("x ",value.(x))
println(value(tau[3]+(t[3][2]+s[3])*value(x[3,2])-l0*(1-value(x[3,2]))-tau[2])<=0)
