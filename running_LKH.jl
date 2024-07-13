using Formatting

# Defining the template for writing onto a .tsp file to run LKH over.
template = """NAME: name
TYPE: ATSP
COMMENT: name
DIMENSION: n_cities
EDGE_WEIGHT_TYPE: EXPLICIT
EDGE_WEIGHT_FORMAT: FULL_MATRIX
EDGE_WEIGHT_SECTION
matrix_sEOF"""

function construct_tour_LKH(matrix, vertices_covered, precision = 1.0)

    # In this function, LKH is run for the considered instance. This is the main function for running LKH.

    # Constructing the .TSP file

    filename = raw"tmp\myRoute.tsp"
    # Opening the .tsp file
    f = open(filename, "w")
    # Writing onto the .tsp file
    temp= create_tsplib_string(matrix)

    
    write(f,temp )
    close(f)
    # Running LKH
    tour = run_lkh(filename, vertices_covered, precision)

    # Saving the obtained route
    route = tour["tour"]


    return route

end

function create_tsplib_string(matrix, name = "route")
    # In this function, the .tsp file corresponding to "matrix" is constructed.

    n_cities = size(matrix)[1]

    io = IOBuffer()     # to store the matrix
    matrix_s = ""
    for i = 1: n_cities
        
        for j in 1:n_cities
            print(io, " ", matrix[i, j])
        end
            println(io) 
    end
    matrix_s = String(take!(io))  # io buffer to string
    
    # Modifying the template
    tmp = template
    tmp = replace(tmp, "name" => name)
    tmp = replace(tmp, "n_cities" => n_cities)
    tmp = replace(tmp, "matrix_s" => matrix_s)

    return tmp
end

function create_lkh_tsp_file_par_out_file(tsp_file_and_path, runs = 4)
    # In this function, the .par file and .out file for LKH will be created.
    # This is required to run LKH for the considered graph.
    
    # Removing the ".tsp" from the name of the file
    prefix, _ = split(tsp_file_and_path, ".")
    # Obtaining the name of the .par and .out files
    par_path = prefix * ".par"
    out_path = prefix * ".out"
    # Constructing the strings to be written onto the .par file
    par = """PROBLEM_FILE = str1
    RUNS = str2
    TOUR_FILE =  str3
    """
    par = replace(par, "str1" => tsp_file_and_path)
    par = replace(par, "str2" => string(runs))
    par = replace(par, "str3" => out_path)
    # Opening the .par file and writing the constructed string
    f = open(par_path, "w")
    write(f, par)
    close(f)

    return par_path, out_path

end

function run_lkh(tsp_file_and_path, vertices_covered, precision = 1.0)
    #= In this function, LKH is run over the constructed .tsp file. Here,
    tsp_file_and_path contains the relative directory in which the .tsp file is
    constructed, and the name of the .tsp file. vertices_covered denotes the list
    of vertices covered by the vehicle, including the depot =#
    
    # Obtaining the paths for the .par and .out file; also creating the .par file
    par_path, out_path = create_lkh_tsp_file_par_out_file(tsp_file_and_path)

    # Running the LKH file
    read(Cmd(["LKH-3", par_path]));

    # Opening the output file, and obtaining the tour
    meta = []
    raw = []
    f = open(out_path, "r")
    # Creating a flag variable to keep track of when the tour starts
    header = true
    lines = readlines(f) # Obtaining all lines in the output file
    solution = nothing

    for line in lines # Running through line by line
        # Checking when the tour starts
        if header
            if startswith(line, "COMMENT : Length = ")
                # Obtaining the cost of the tour, and scaling down using the
                # precision parameter. Note that the precision parameter was
                # used for integerizing the edge weights.
                solution = parse(BigInt, (last(split(line, " "))))/(10^precision)
            end
            if startswith(line, "TOUR_SECTION")
                header = false
                continue
            end
            append!(meta, line)
            # println(line)
            continue
        else
            if startswith(line, "-1") # Checking when the tour stops
                append!(raw, vertices_covered[1])
                break
            else
                append!(raw, vertices_covered[parse(Int64, strip(line, [' ']))])
                 # Appending the vertex, after removing
                # unnecessary spaces. Note that the obtained vertex is mapped
                # with "vertices_covered" here. 
                # The tour obtained from LKH corresponds to vertices being covered being 1, 2, ....
                # However, if our list of vertices to be covered is different from 1, 2, ..., we
                # need to match the corresponding vertices. For example, "1" in the tour will be
                # matched with vertices_covered[1], "2" will be matched with vertices_covered[2],
                # and so on.
            end
        end
    end

    # Collecting all initial information about the graph for which tour was constructed
    metadata = join(meta, "")

    # Closing the file
    close(f)

    # # Converting the obtained tour by matching to "vertices_covered". Note that
    # # the tour obtained from LKH corresponds to vertices being covered being 1, 2, ....
    # # However, if our list of vertices to be covered is different from 1, 2, ..., we
    # # need to match the corresponding vertices. For example, "1" in the tour will be
    # # matched with vertices_covered[1], "2" will be matched with vertices_covered[2],
    # # and so on.
    # raw_modified = []
    # for i in raw
    #     append!(raw_modified, vertices_covered[i])
    # end
    # # println(raw_modified)

    return Dict("tour" => raw, "solution" => solution, "metadata" => metadata)

end