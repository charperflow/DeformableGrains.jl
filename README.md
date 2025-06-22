# DeformableGrains.jl

A Julia package for simulating and visualizing chains of colliding deformable particles.

## Installation

You can install the package directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/charperflow/DeformableGrains.jl")
```

# Usage Example
using DeformableGrains

N = 6
v = 1
K = [50.0,1000.0,50.0,1000.0,50.0,1000.0]
M = [1.0, 0.25, 1.0, 0.25, 1.0, 0.25]
R = [1,0.5,1,0.5,1,0.5]
time = 2

Cs, T, Xs, Vs, Ws = compression_chain(N, v, K, M, R, time; g=0, Δt=0.0001, stuck=false)

# License
MIT 
