
using Dubins
using Plots

include("running_LKH.jl")
include("Dubins_interval.jl")
#include("data.jl")




function min_cost(orig_coord,orig_head,dest_coord,dest_head)    # cost to go from one point to the another point

    orig_val = orig_head # starting heading
    new_coord_orig = copy(orig_coord) # starting coordinate
    push!(new_coord_orig, orig_val) # starting coordinate plus heading
    
    dest_val = dest_head    # destination heading
    new_coord_dest = copy(dest_coord)   # destination coordinate
    push!(new_coord_dest, dest_val) # destination coordinate plus heading

    errcode, path = dubins_shortest_path(new_coord_orig, new_coord_dest, 1.)
    val = dubins_path_length(path)  # cost of the path travelled

    result = [val,orig_val,dest_val]    # cost, starting heading and final heading
        
    return result

end



function solve_dubins_tsp(coordinates,N, type,r)        #function to give upper and lower bound for the given coordinates

    n = size(coordinates)[1]        # number of targets

    M = fill(0, N*n, N*n)       # initialize the cost matrix


    # First, we obtain the cost to go from one target to another
    for i in 1:n*N # We go through every row in the cost matrix
        
        for j in 1:n*N # We now go through every column in the cost matrix

            V1 = (i-1) % N + 1 # Obtaining the vertex number corresponding to the starting target (corresponds to exiting angle)
            V2 = (j-1) % N + 1 # Obtaining the vertex number corresponding to the ending target     (corresponds to entering angle)
            C1 = div(i-1, N) + 1 # Obtaining the target number for starting target  
            C2 = div(j-1, N) + 1  # Obtaining the target number for ending target  
            V1 = (V1 + 1) % N        # adjusting the starting vertex to obtain the right cost 
            
            if C1 != C2 # We do not want to consider vertices within the same target
                
                if type == 1 # feasible solution cost
                
                    result = min_cost(coordinates[C1],(V1-1)*2*pi/N,coordinates[C2],(V2-1)*2*pi/N)
                
                elseif type==2 # lower bound cost

                    result = Dubins_interval(coordinates[C1],[(V1-1)*2*pi/N V1*2*pi/N],coordinates[C2],[(V2-1)*2*pi/N V2*2*pi/N],r)
                
                end

                M[i,j] = Int64(round(1000*result[1])) # cost assigned as per the type of solution required

            end
                
        end
    end


    Big_M=Int64(maximum(M)*n)   # calculation of Big M

    
    for i in 1:n*N
        for j in 1:n*N
            M[i,j]=M[i,j]+Big_M # augmenting the cost matrix with Big_M
        end
    end

    for i in 1:n # Running through a loop for constructing cost for going from one vertex to another within the same target
        
        M[(i-1)*N+N,(i-1)*N+1] = 0 # Connecting last vertex in the selected target to the first vertex in the same target with zero cost
        
        for j in 1:N-1 # Running a for loop to connect all intermediary vertices with the succeeding vertex with a zero cost
            M[(i-1)*N+j,(i-1)*N+j+1] = 0
        end

    end



    opt_tour = construct_tour_LKH(M, collect(1:size(M)[1]), 1.0) # run LKH to solve ATSP
    opt_tour = opt_tour[1:end-1]    # first and last point are same

    vertex = []
    coord = []
    
    for i in 1:size(opt_tour)[1]    # extracting coordinates and vertices from the optimal tour
        V1 = (opt_tour[i]-1) % N + 1        # vertex number
        C1 = div(opt_tour[i]-1,N) + 1       # target number
        push!(vertex,V1)
        push!(coord,C1)
    end


    temp = coord[1]
    global k
    k = 0 

    for i in 1 : N # find the number of vertices after which route exists the first target, that vertex represents the angle of entry of the corrosponding target
        if  temp == coord[i] 
            global k
            k = k+1
        end
    end

 

    vertex1 = []
    coord1 = []

    for i in 1:size(vertex)[1]
        if (i + (N - k)) % N == 1       # (N - k + 1) vertex represents angle of entry of the corrosponding target, hence every Nth vertex after that will also be angle of entry for the corrosponding target 
            push!(vertex1,vertex[i])
            push!(coord1,coord[i])
        end
    end

        # extract angles and cost 

    angle1 = []     # exit angle from target 1
    angle2 = []     # entrying angle at the target 2
    cost = []       # cost of traveling from target1 to target 2
    for i in 1:size(coord1)[1]-1
        C1 = coord1[i]
        C2 = coord1[i+1]
        V1 = vertex1[i]
        V2 = vertex1[i+1]

        if type == 1
                
            result = min_cost(coordinates[C1],(V1-1)*2*pi/N,coordinates[C2],(V2-1)*2*pi/N)      #  feasible solution for upper bound
        
        elseif type==2

            result = Dubins_interval(coordinates[C1],[(V1-1)*2*pi/N V1*2*pi/N],coordinates[C2],[(V2-1)*2*pi/N V2*2*pi/N],1.0)   # lower bound cost
        
        end
        
        push!(cost,result[1])
        push!(angle1,result[2])
        push!(angle2,result[3])
        
    end
        C1 = coord1[size(coord1)[1]]    # coordinate of last target
        C2 = coord1[1]      #   coordinate of first target
        V1 = vertex1[size(coord1)[1]]   # vertex of last target
        V2 = vertex1[1] # vertex of first target

        if type == 1
                
            result = min_cost(coordinates[C1],(V1-1)*2*pi/N,coordinates[C2],(V2-1)*2*pi/N) #  feasible solution for upper bound
        
        elseif type==2

            result = Dubins_interval(coordinates[C1],[(V1-1)*2*pi/N V1*2*pi/N],coordinates[C2],[(V2-1)*2*pi/N V2*2*pi/N],1.0) # lower bound cost
        
        end

        push!(cost,result[1])
        push!(angle1,result[2])
        push!(angle2,result[3])
        

    return coord1, angle1, angle2, sum(cost)    # return the sequence of targets to visit, angle of exit and entry and total cost of the path

end


function plot_dubins_tsp(coordinates,N,type,r)

    coord,angle1,angle2=solve_dubins_tsp(coordinates,N,type,r)      # solves the Dubins TSP

    x_coord=[]  # coordinate of path will be stored here
    y_coord=[]  # coordinate of path will be stored here
    x_end=[]    # coordinate of targets are stored here
    y_end=[]    # coordinate of targets are stored here
    value=[]    # cost of travelling will be stored here
    for i in 1:size(angle1)[1]-1
    
        # now we find the dubins path going from ith target to (i+1)th target

        errcode, path = dubins_shortest_path([coordinates[coord[i]][1],coordinates[coord[i]][2],angle1[i]], [coordinates[coord[i+1]][1],coordinates[coord[i+1]][2],angle2[i]], 1.)
        val = dubins_path_length(path)
        push!(value,val)
        
        global k 
        k = 0.0

        while k < val   # k is used to sample the path 
            global k
            err, temp=dubins_path_sample(path, k)
            push!(x_coord,temp[1])
            push!(y_coord,temp[2])
            k = k+0.1
        end
            push!(x_end,coordinates[coord[i]][1])
            push!(y_end,coordinates[coord[i]][2])
        
    end
    
    d=size(angle1)[1]      # last target
    
        errcode, path = dubins_shortest_path([coordinates[coord[d]][1],coordinates[coord[d]][2],angle1[d]], [coordinates[coord[1]][1],coordinates[coord[1]][2],angle2[d]], 1.)
        val = dubins_path_length(path)
        push!(value,val)
        
        global k 
        k = 0.0
        
        while k<val
            global k
            err, temp=dubins_path_sample(path, k)
            push!(x_coord,temp[1])
            push!(y_coord,temp[2])
       
            k=k+0.1
            
        end
    
            push!(x_end,coordinates[coord[d]][1])
            push!(y_end,coordinates[coord[d]][2])
    
    
        
    println(sum(value))
    plot(x_coord,y_coord, label="path")
    scatter!(x_end,y_end,label="target")
    


end