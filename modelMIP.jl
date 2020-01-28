#groupe 2: J.Zribi, G.Dekeyser, P.Parmentier

using JuMP
using GLPK

m = Model(with_optimizer(GLPK.Optimizer))

n=4
nArc=6
r=1.0
g=1.0
Q=7
l0 =1000000

sommets=[[1,0,10,0],[2,2,10,0],[3,0,10,0],[4,0,15,0]]
#fenetre temps debut - fin - type sommet 0=client 1=station de recharge
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
print(nClients)

print(tabClients)

@constraint(m,cons2[i in 1:nClients ],sum(x[tabClients[i],j] for j in 2:n if j!=tabClients[i])==1)
@constraint(m,cons3[i in 1:nStations ],sum(x[tabStations[i],j] for j in 2:n if j!=tabStations[i])<=1)

@constraint(m,cons4[j in 2:(n-1)], sum(x[j,i] for i in 2:n if i!=j)-sum(x[i,j] for i in 1:(n-1) if i!=j)==0)

#contrainte 5
for i in 1:(n-1)
    if sommets[i][4]==0
        for j in 2:n
            if i!=j
                @constraint(m,tau[i]+(t[i][j]+s[i])*x[i,j]-l0*(1-x[i,j])<=tau[j])
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
                @constraint(m,0<=y[j]<=y[i]-(r*d[i][j])*x[i,j]+Q*(1-x[i,j]))
            end
        end
    end
end
