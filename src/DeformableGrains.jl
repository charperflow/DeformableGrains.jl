module DeformableGrains

greet() = print("Hello World!") #This is cute
include("helper.jl")
include("collision_model.jl")
include("visualize_results.jl")

export compression_chain, animate_particles

end # module DeformableGrains
