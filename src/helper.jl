#helper_vf3_4.jl
"""
    setup_dust_settled(n, v, K, M, R, h; g=9.81)

Simulates the initial setup for a system of dust particles settling under gravity and interacting through a given stiffness. 

# Arguments
- `n::Int`: The number of particles in the system.
- `v::Float64`: The initial velocity of the first particle.
- `K::Vector{Float64}`: A vector of stiffness constants for the particles.
- `M::Vector{Float64}`: A vector of masses for the particles.
- `R::Vector{Float64}`: A vector of radii for the particles.
- `h::Int`: The number of time steps to simulate.
- `g::Float64=9.81`: Gravitational acceleration (default is Earth's gravity).

# Returns
A tuple `(Xs, Vs)` where:
- `Xs::Array{Float64, 2}`: A 2D array storing the positions of particles over time. The shape is `n × h`.
- `Vs::Array{Float64, 2}`: A 2D array storing the velocities of particles over time. The shape is `n × h`.

# Example
```julia
n = 5
v = 1.0
K = [10.0, 15.0, 20.0, 25.0, 30.0]
M = [2.0, 3.0, 4.0, 5.0, 6.0]
R = [0.5, 0.6, 0.7, 0.8, 0.9]
h = 100

Xs, Vs = setup_dust_settled(n, v, K, M, R, h)
println("Initial positions: ", Xs[:, 1])
println("Initial velocities: ", Vs[:, 1])

"""
function setup_dust_settled(n,v,K,M,R,h;g=9.81)
    X = zeros(n)
    d = [2*R[i] for i in 1:n]
    

    Δx = zeros(n)
    Δx[n] = sum(M[1:n-1])*g/K[n]
    for i in 1:n-1
        Δx[i] = (sum(M[1:i])+sum(M[1:i-1]))*g/K[i]
    end

    X[n] = sum(d) - ((d[1]/2) + d[n]/2)
    d = d .- Δx   
    for i in n-1:-1:1
        X[i]=X[i+1]-(d[i+1]/2)-(d[i]/2)
    end

    V = zeros(n)
    V[1] = v
    #Arrays of positions and velocities and setting ICs
    Xs = zeros(n,h)
    Vs = zeros(n,h)
    Xs[:,1] = X
    Vs[:,1]= V

    return Xs, Vs

end #function

"Intermediate computation for the particle boundaries"
function compute_recursive_k(arr::Vector{Float64})
    # Base case: if the array has only one element, return 1
    if length(arr) == 1
        return 1.0
    else
        # Recursive call for the subarray excluding the last element
        previous_sum = compute_recursive_k(arr[1:end-1])
        # Compute the last term as a product of ratios
        last_term = prod(arr[i+1] / arr[i] for i in 1:length(arr)-1)
        # Add the last term to the recursive sum
        #println("$previous_sum + $last_term")
        return previous_sum + last_term
    end
end

"""
    solve_b0(n,chain_positions,chain_radii,chain_softness)

Solve for the leftmost boundary of a chain of particles.

Given the position, radius, and stiffness of each particle 
(`chain_positions`, `chain_radii`, `chain_softness`) in an `n`-particle chain, 
this function computes the position of the leftmost boundary.
"""
function solve_b0(n,chain_positions,chain_radii,chain_softness)
    κ = compute_recursive_k(chain_softness)
    term_1 = 0
    term_2 = 0
    #loop for term_1
    for i in 1:n-1
        c_1 = (-1)^(i+1)
        prod_1 = prod([chain_softness[j+1]/chain_softness[j] for j in 1:i])
        diff_1 = chain_radii[i+1] - chain_positions[i+1]
        term_1 = term_1 + c_1*prod_1*diff_1

        prod_2 = prod([chain_softness[j+1]/chain_softness[j] for j in 1:i])
        diff_2 = sum([(-1)^(k+1)*2*chain_positions[k] for k in 1:i])
        term_2 = term_2 + prod_2*diff_2
    end

    b_0 = (1/κ)*(term_1+term_2-chain_radii[1]+chain_positions[1])

    return b_0
end #function


"""
satisfies_newton_third_law(R_sub::Vector{Float64},
                            R_deformed::Vector{Float64},
                            K_sub::Vector{Float64};
                            tol::Float64 = 1e-8)::Bool

Ensure that the partial compression between the `i`th and `(i+1)`th particles satisfies Newton's third law.

In a chain of `n > 2` particles, the total stiffness-scaled deformation between two particles 
does not have to be equal. However, their **stiffness-scaled partial compressions** 
(the component of deformation each particle contributes to the other) must be equivalent.

# Arguments
- `R_sub::Vector{Float64}`: Relaxed radii for the current subset of contiguous particles.
- `R_deformed::Vector{Float64}`: Deformed radii for the same subset.
- `K_sub::Vector{Float64}`: Stiffness values for each particle in the subset.
- `tol::Float64 = 1e-8`: Tolerance used to ignore machine-precision-level overlaps.

# Returns
- `Bool`: `true` if Newton's third law is satisfied across the subset, `false` otherwise.
"""
function satisfies_newton_third_law(R_sub::Vector{Float64},
                                    R_deformed::Vector{Float64},
                                    K_sub::Vector{Float64};
                                    tol::Float64 = 1e-8)::Bool
    n = length(R_sub)
    δ = [R_sub[i] - R_deformed[i] for i in 1:n]

    # First particle: total compression is all due to contact with the right
    δ_right = δ[1]

    for i in 2:n
        # Left compression on i from (i-1)
        δ_left = δ_right * (K_sub[i-1] / K_sub[i])

        if δ_left > δ[i] + tol
            return false
        end

        # Remaining compression must go to the next contact (i.e., δ[i] - δ_left)
        δ_right = δ[i] - δ_left

        if δ_right < -tol
            return false
        end
    end

    return true
end

"""
find_chains(X::Vector{Float64}, R::Vector{Float64}, K::Vector{Float64}; tol::Float64 = 1e-8)

Identify which subsets of particles form chains.

Given particle positions `X`, radii `R`, and stiffnesses `K`, partition a system of `n` particles 
into subsets representing compression chains—groups of particles whose contact forces 
result in mutual deformation. The `tol` parameter accounts for machine precision when determining contact.

# Arguments
- `X::Vector{Float64}`: Center positions of each particle.
- `R::Vector{Float64}`: Relaxed (undeformed) radii of each particle.
- `K::Vector{Float64}`: Stiffness of each particle.
- `tol::Float64 = 1e-8`: Minimum overlap required to count as deformation.

# Returns
- `best_part::Vector{Vector{Int64}}`: The minimal-overlap partition of the particles into chains.
- `best_B::Dict{Int64, Vector{Float64}}`: Maps the chain index (given by first particle in each chain) to its boundary positions.
- `best_W::Vector{Float64}`: Deformed width of each particle.

# Examples
```jldoctest
julia> X = [1.0, 2.0, 5.0, 5.5];

julia> R = [1.0, 1.0, 2.0, 2.0];

julia> K = [1000.0, 1000.0, 1000.0, 1000.0];

julia> find_chains(X, R, K)
([[1, 2], [3, 4]], Dict(3 => [4.75, 5.25, 5.75], 1 => [0.5, 1.5, 2.5]), [1.0, 1.0, 0.5, 0.5])
```
"""
function find_chains(X::Vector{Float64},
                     R::Vector{Float64},
                     K::Vector{Float64};
                     tol::Float64 = 1e-8)

    N = length(X)
    best_part   = nothing
    best_B      = Dict{Int, Vector{Float64}}()
    best_W      = similar(R)
    best_overlap = Inf

    # Try every way to cut between particles
    for mask in 0:(1 << (N-1)) - 1
        # 1) Build the candidate partition
        candidate = Vector{Vector{Int}}()
        start = 1
        for j in 1:(N-1)
            if (mask >> (j-1)) & 1 == 1
                push!(candidate, collect(start:j))
                start = j+1
            end
        end
        push!(candidate, collect(start:N))

        # 2) Compute B and W for each chain, checking physics
        B_tmp = Dict{Int, Vector{Float64}}()
        W_tmp = zeros(Float64, N)
        valid = true

        for ch in candidate
            pos = X[ch]; rad = R[ch]; ksub = K[ch]
            if length(ch) == 1
                i = ch[1]
                B_tmp[i] = [X[i]-R[i], X[i]+R[i]]
                W_tmp[i] = 2R[i]
            else
                b0 = solve_b0(length(ch), pos, rad, ksub)
                n = length(ch)
                left = zeros(n); right = zeros(n)
                left[1] = b0
                for m in 1:(n-1)
                    left[m+1] = pos[m] + (pos[m] - left[m])
                    right[m]   = left[m+1]
                end
                right[end] = pos[end] + (pos[end] - left[end])
                Rdef = (right .- left) ./ 2

                # physics checks
                if any(Rdef .> rad .- tol) ||
                   !satisfies_newton_third_law(rad, Rdef, ksub; tol=tol)
                    valid = false
                    break
                end

                B_tmp[ch[1]] = vcat(left, right[end])
                W_tmp[ch]    = 2 .* Rdef
            end
        end
        if !valid
            continue
        end

        # 3) Compute relative overlap between adjacent chains
        overlap = 0.0
        count = 0
        for idx in 1:(length(candidate)-1)
            c1 = candidate[idx]
            c2 = candidate[idx+1]
            right1 = B_tmp[c1[1]][end]
            left2  = B_tmp[c2[1]][1]
            count += 1

            raw_overlap = max(0, right1 - left2)
            if raw_overlap > tol
                i1 = c1[end]  # last particle of chain 1
                i2 = c2[1]    # first particle of chain 2
                denom = R[i1] + R[i2]
                overlap += raw_overlap / denom
            end
            overlap =overlap/count
        end

        # 4) Keep the partition with minimal overlap
        if overlap < best_overlap - tol
            best_overlap = overlap
            best_part    = deepcopy(candidate)
            best_B       = deepcopy(B_tmp)
            best_W       = copy(W_tmp)
            # if perfect, stop early
            if overlap < tol
                break
            end
        end
    end

    return best_part, best_B, best_W
end

"""
    calc_forces(chain,B,X,V,K,M,R,γ)

Calculate the force on each particle.

Using a linear-dashpot model calculate the interparticular force and acceleratio of each
particle. The calculation requires the specified chain as a list of particle indexes `chain`,
the particle widths `B`, positions `X`, velocities `V`, stiffnesses `K`, masses `M`, radius
`R` and dissipation constant `γ`. 
"""
function calc_forces(chain,B,X,V,K,M,R,γ)
    mass = M
    M = size(X)[1]
    N = size(B)[1] -1
    left = chain[1]
    right= last(chain)
    squish = zeros(N)
    force_L = zeros(N) #scale squish by k/m in future cases and use these
    force_R = zeros(N)
    squish_L = zeros(N)
    squish_R = zeros(N)
    F        = zeros(M)
    A        = zeros(M)
    vs       =zeros(M)
    #rel vs left
    vs[left+1:right] = abs.(diff(V[left:right]))
    #add the rel vs to the right particle
    vs[left:right-1] = vs[left:right-1]+ abs.(diff(V[left:right])) 
    #println("Spring vs are $vs")

    XS = [X[i] for i in left:right]
    Rc = [R[i] for i in left:right]
    Kc = [K[i] for i in left:right]
    squish = [2*Rc[i]-(B[i+1]-B[i]) for i in 1:N]
    squish_L[1] = 0
    force_L[1] = 0
 
    for i in 1:N-1
        squish_R[i]=min(squish[i]-squish_L[i],2*Rc[i]) #particle cant be compressed more than its total diameter
        force_R[i] =(Kc[i])*(squish_R[i])
        force_L[i+1] = force_R[i]
        squish_L[i+1]=(1/Kc[i+1])*force_R[i]
    end
    squish_R[N]=0
    force_R[N] = 0 
    F[left:right] = force_L-force_R # these would be forces in k,n ≠ 1 case
    A = [F[i]/mass[i] for i in 1:M]
    #println("A is $A")
    return A

end#function


"""
    F(I,K,M,R,W,g;γ=0)

The master function for the RK4.
"""
function F(I,K,M,R,W,g;γ=0)
    #println("we're in F")
    N = Int(size(I)[1]/2)
    X = I[1:N]
    V = I[N+1:2*N]

    dX = V
    dV = zeros(N)
    #println("we're diving into find_chains")
    chains, B_dict, W = find_chains(X,R,K)
    J = size(chains)[1]
    #println("we're about to calculate forces, J is $J")
    for j in 1:J
        chain = chains[j]
        left = chain[1]
        right = last(chain)
        B = B_dict[left]
        
        A = calc_forces(chain,B,X,V,K,M,R,γ)
        dV = dV + A
    end
    dV = dV .+ g #constant of acceleration for gravity
    
    return [dX; dV]

end#function


