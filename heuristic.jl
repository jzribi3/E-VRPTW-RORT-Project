using Dates

t0=now()


struct Vertex                 # either a customer or a charging station
    i::Int32
    window_start
    window_end
    is_charging_station
    service_time              # only specified for customers
end

struct Instance
    n                         # nb of vertexes
    battery_consumption_rate  # per distance
    battery_recharging_rate   # time per unit of battery level actually
    battery_capacity
    V                         # all customers and charging stations
    D                         # distances matrix
    T                         # travel times matrix
end

struct Stop
    i::Int32
    arrival_time
    departing_time
    battery_level            # level of the battery at arrival_time
end

struct Solution
    routes                   # an array of Stop objects
end

function is_possible(stop, I)
    #if stop.departing_time > I.V[stop.i].window_end
    #    println("not possible because departing time would be ", stop.departing_time, " > ", I.V[stop.i].window_end)
    #elseif stop.battery_level < 0
    #    println("not possible because of battery level")
    #end
    return stop.departing_time <= I.V[stop.i].window_end && stop.battery_level >= 0
end

function is_route_possible(route, I)
    for position = 1 : length(route)
        if !is_possible(route[position], I)
            return false
        end
    end
    return true
end

function is_solution_possible(S, I)
    for k = 1 : length(S.routes)
        if !is_route_possible(S.routes[k], I)
            return false
        end
    end
    return true
end

function compute_next_stop(stop, u, I)
    arrival_time = stop.departing_time + I.T[stop.i][u.i]
    # if stop is at charging station then take into account the charging
    if (I.V[stop.i].is_charging_station)
        previous_battery_level = I.battery_capacity
    else
        previous_battery_level = stop.battery_level
    end
    battery_level = previous_battery_level - I.battery_consumption_rate * I.D[stop.i][u.i]
    if u.is_charging_station
        departing_time = arrival_time + I.battery_recharging_rate * (I.battery_capacity - stop.battery_level)
    elseif u.i == I.n
        departing_time = arrival_time  # does not matter for arrival depot
    else
        departing_time = max(arrival_time, u.window_start) + u.service_time
    end
    return Stop(
        convert(Int32, u.i),
        arrival_time,
        departing_time,
        battery_level
    )
end

function get_route_nb(S)
    return length(S.routes)
end

function get_route_distance(route, I)
    res = 0
    for s = 1 : length(route) - 1
        res += I.D[route[s].i][route[s + 1].i]
    end
    return res
end

function get_traveled_distance(S, I)
    res = 0
    for k = 1 : length(S.routes)
        res += get_route_distance(S.routes[k], I)
    end
    return res
end

function find_best_insertion!(route, u, I)
    # u is assumed to be a customer
    # loop over all possible insertion points
    found_position = false
    best_position = -1
    best_cost = Inf
    for position = 2:length(route)
        #println("test insert at position ", position)
        # check arrival time and battery level at customer u is fine
        next_stop = compute_next_stop(route[position - 1], u, I)
        if !is_possible(next_stop, I)
            continue
        end
        # check the next stops are all possible, up to the arrival depot
        is_possible_ = true
        for s = position:length(route)
            next_stop = compute_next_stop(next_stop, I.V[route[s].i], I)
            if !is_possible(next_stop, I)
                is_possible_ = false
                break
            end
        end
        if (!is_possible_)
            continue
        end
        #println("found insertion position")
        found_position = true
        # update best_position if the cost of this insertion is lower
        # than at previous positions
        cost = I.D[route[position - 1].i][u.i] + I.D[u.i][route[position].i] - I.D[route[position - 1].i][route[position].i]
        if (cost < best_cost)
            best_cost = cost
            best_position = position
        end
    end
    return found_position, best_position, best_cost
end

function insert!(route, u, position, I)
    # insert vertex u in route at given position
    # shift (route[position], ..., route[len(route)]) towards the right
    next_i = route[position].i
    route[position] = compute_next_stop(route[position - 1], u, I)
    for s = position + 1 : length(route)
        new_stop = compute_next_stop(route[s - 1], I.V[next_i], I)
        next_i = route[s].i
        route[s] = new_stop
    end
    append!(route, [compute_next_stop(route[length(route)], I.V[next_i], I)])
end

function get_all_permutations(list)
    if length(list) == 1
        return [[list[1]]]
    end
    res = []
    for i = 1 : length(list)
        sub_permutations = get_all_permutations(vcat(list[1:i-1], list[i+1:length(list)]))
        for s = 1 : length(sub_permutations)
            append!(res, [vcat(list[i], sub_permutations[s])])
        end
    end
    return res
end

function build_route(vertexes, I)
    # departing depot
    route = [Stop(1, 0, 0, I.battery_capacity)]
    # vertexes = customers and charging stations
    for m = 1 : length(vertexes)
        append!(route, [compute_next_stop(route[length(route)], vertexes[m], I)])
    end
    # arrival depot
    append!(route, [compute_next_stop(route[length(route)], I.V[I.n], I)])
end

function arrange_into_feasible_route(vertexes, I)
    # try every permutation of the vertexes
    # only tractable for very small number of vertexes
    permutations = get_all_permutations(vertexes)
    for i = 1 : length(permutations)
        route = build_route(permutations[i], I)
        if (is_route_possible(route, I))
            return true, route
        end
    end
    return false, []
end

function build_greedy_solution(I)
    routes = []
    is_already_assigned = falses(I.n)
    # loop over every customer in order to assign them to a route
    for i = 2 : I.n - 1
        if I.V[i].is_charging_station || is_already_assigned[i] continue end
        # loop over every route in order to find one where customer
        # can be inserted with minimum cost
        found_route = false
        best_cost = Inf
        best_position = -1
        best_route = -1
        for k = 1:length(routes)
            found_position, position, cost = find_best_insertion!(routes[k], I.V[i], I)
            if found_position && cost < best_cost
                found_route = true
                best_cost = cost
                best_position = position
                best_route = routes[k]
            end
        end
        if found_route
            insert!(best_route, I.V[i], best_position, I)
        else
            # it can happen for battery level/time window reasons or absence
            # of routes that there is no route where the customer i could be
            # inserted; in this case we make a new route for him
            new_route = build_route([I.V[i]], I)
            # check the route is feasible; sometimes the customer is so
            # far from the depot that we will need to use charging stations
            # or even other customers depots when the graph is not complete
            if !is_route_possible(new_route, I)
                # first see whether adding one charging station is enough
                # loop over charging stations
                for j = 2 : I.n - 1
                    if !I.V[j].is_charging_station continue end
                    found_position, position, cost = find_best_insertion!(new_route, I.V[j], I)
                    if (found_position)
                        insert!(new_route, I.V[j], position, I)
                        break
                    end
                end
            end
            # find out if just this worked
            if length(new_route) == 3
                # it did not so we try adding two stops at the same charging station
                for j = 2 : I.n - 1
                    if !I.V[j].is_charging_station continue end
                    new_route = build_route([I.V[j], I.V[i], I.V[j]], I)
                    if is_route_possible(new_route, I)
                        break
                    end
                end
            end
            # find out if this worked
            if !is_route_possible(new_route, I)
                # it did not so try all the possibilities with customer i, another customer j and
                # maybe one or two charging station stops
                for j = 2 : I.n - 1
                    if I.V[j].is_charging_station || i == j continue end
                    feasible, new_route = arrange_into_feasible_route([I.V[i], I.V[j]], I)
                    if feasible
                        is_already_assigned[j] = true
                        break
                    end
                    for k = 2 : I.n - 1
                        if !I.V[k].is_charging_station continue end
                        feasible, new_route = arrange_into_feasible_route([I.V[i], I.V[j], I.V[k]], I)
                        if feasible
                            is_already_assigned[j] = true
                            break
                        end
                        feasible, new_route = arrange_into_feasible_route([I.V[i], I.V[j], I.V[k], I.V[k]], I)
                        if feasible
                            is_already_assigned[j] = true
                            break
                        end
                    end
                    if feasible break end
                end
            end
            append!(routes, [new_route])
        end
        is_already_assigned[i] = true
    end
    return Solution(routes)
end

function build_instance(fileName)
    filePath=string(@__DIR__,"/data/")
    n,nArcs,tauxConso,tauxRech,capa,sommets_array,d,t = readInstance(filePath,fileName)
    println("n,nArcs,tauxConso,tauxRech,capa,t= ",n,",",nArcs,",",tauxConso,",",tauxRech,",",capa,",",t)
    V = []
    for i = 1 : length(sommets_array)
        append!(V, [Vertex(sommets_array[i][1], sommets_array[i][2], sommets_array[i][3], sommets_array[i][4] == 1, 0)])
    end
    #println("V: ", V)
    #println("D: ", d)
    #println("T: ", t)
    I = Instance(
        n,
        tauxConso,
        tauxRech,
        capa,
        V,
        d,
        t
    )
    return I
end

function build_routes_maybe_unfeasible(I, BIG_M)
    xk=[]
    dk=[]
    for i in 2:(I.n-1)
        if I.V[i].is_charging_station continue end
        xki=[[0 for j in 1:I.n] for k in 1:I.n]
        xki[1][i]=1
        xki[i][I.n]=1
        push!(xk,xki)
        new_route = build_route([I.V[i]], I)
        if is_route_possible(new_route, I)
            push!(dk,get_route_distance(new_route, I))
        else
            push!(dk,BIG_M)
        end
    end
    return dk, xk
end

#I = Instance(
#    4,
#    1,
#    1,
#    7,
#    [Vertex(1, 0, 1, false, 0), Vertex(2, 2, 10, false, 0), Vertex(3, 0, 10, false, 0), Vertex(4, 0, 15, false, 0)],
#    [[0,4,2,1000],[1000,0,2,2],[1000,1,0,4],[1000,1000,1000,0]],
#    [[0,1,5,1000],[1000,0,4,5],[1000,3,0,1],[1000,1000,1000,0]]
#)

include("lireInstance.jl")

#fileName = "E_data_3.txt"
fileName = "evrptw_instances/c103C15.txt"
I = build_instance(fileName)

# dk, xk = build_routes_maybe_unfeasible(I, 1000)
# println(dk)
# println(xk)


S = build_greedy_solution(I)
if !is_solution_possible(S, I)
    println("There is a problem with the solution, it is not feasible")
end
println("Number of routes used: ", get_route_nb(S))
println("Distance traveled: ", get_traveled_distance(S, I))
println("Routes:")

for k = 1:length(S.routes)
    println(S.routes[k])
end


tend=now()
println("deltat=",tend-t0)
