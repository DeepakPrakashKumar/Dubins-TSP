#  Dubins-TSP
## Overview

This package provides functionality to find the lower and upper bound for the Traveling Salesman Problem(TSP) with a Dubin's vehicle. 

## Installation
To use this package, include it in your Julia environment:

```julia
include("Dubins_TSP.jl")
```

## Usage

The primary functions provided by this package are **'solve_dubins_tsp'** and **'plot_dubins_tsp'**

```julia
include("Dubins_TSP.jl")

n=10  # number of coordinates
coordinates = [rand(2)*10 for _ in 1:n]

N= 10  # number of intervals per coordinate
type=1 # 1: upper bound, 2: lower bound
r= 1.0 #   minimum turning radius

tour,angle1,angle2,cost= solve_dubins_tsp(coordinates,N,type,r)

plot_dubins_tsp(coordinates,N,type,r)


```

### Parameters:
* **coordinates**: An array of coordinates indicating the locations of the targets to be visited.
* **N**:  The number of intervals in which the heading will be divided for each target (resolution).
* **type**:  The type of solution, either lower bound or upper bound.
* **r**:  The vehicle's minimum turning radius.

### Output

#### solve_dubins_tsp
* **tour**: The sequence of targets to be visited.
* **angle1**: The heading angle of the vehicle upon entering each target.
* **angle2**: The heading angle of the vehicle upon exiting each target.
* **cost**: The total cost of visiting all the targets.

#### plot_dubins_tsp
This function plots the solution and displays the total cost of visiting all the targets.


## Contact
If you have any questions or issues, please open an issue on this repository or contact the author.
