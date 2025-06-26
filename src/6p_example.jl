#6p_example.jl
#===
We start by loading packages that we will use in postprocessing
=#

using Plots, DataFrames, CSV 

#========
Setting up and running the collison model. 
The inputs are:
- N::Int  number of particles
- v::Float64  initial velocity of the leftmost particle (m/s)
- K::Vector{Float64} stiffness of the particles (N/m) 
- M::Vector{Float64} particles masses (kg)
- R::Vector{Float64} particle raddii

- Cs::Vector{Any} contains the contact structures of every particle
- T::Vector{Float64} is a vector containing every timestep
- Xs::Vector{Float64} contain the positions at each timestep
- Vs::Vector{Float64} particle velocities at every timestep
- Ws::Vector{Float64} the widths at every timestep
========#
N = 6 
v = 1
K = [100,100,1000.0,1000.0,1000.0,100]
M = [1.0,1.0,0.25,0.25,0.25,1.0]
R = [1.0,1.0,.25,.25,.25,1.0]
time = 2
Cs,T,Xs,Vs,Ws = compression_chain(N,v,K,M,R,time;g=0,Δt=0.0001,stuck=false);

#=====================
Helpful plots!
In general this also shows a good way to access the relevant data if you 
want to perform any additional postprocessing.
=====================#

plot_velocity = plot(T,[Vs[1,:] Vs[2,:] Vs[3,:] Vs[4,:] Vs[5,:] Vs[6,:]],
	xlabel = "time (s)",
	ylabel = "vlocity (m/s)",
	title = "Velocities no g - fixed end",
	xlim=(0,2),
	linewidth=3,
	label=["particle 1" "particle 2" "particle 3" "particle 4" "particle 5" "particle 6"],
	legend=:topright)

plot_widths = plot(T,[Ws[1,1] .- Ws[1,:] Ws[2,1] .- Ws[2,:] Ws[3,1] .-  Ws[3,:] Ws[4,1] .- Ws[4,:] Ws[5,1] .- Ws[5,:] Ws[6,1] .-  Ws[6,:] ],
	xlabel = "time (s)",
	ylabel = "deformation",
	title = "Deformations no g - fixed end",
	xlim=(0,2),
	linewidth=5,
	linestyle=:dash,
	alpha = 0.6,
	label=["particle 1" "particle 2" "particle 3" "particle 4" "particle 5" "particle 6"],
	legend=:topright)


C_x1 = [c[1] for c in filter(c->!isempty(c),Cs)]
contacts_x1 = [length(c)-1 for c in C_x1]
plot(T[2:end],contacts_x1,
	xlim = (0,1.5),
	linewidth=3,
	title = "Length of First Compression Chain Through Time",
	xlabel = "time (s)",
	ylabel = "Length of X_c1",
	legend = :false
	)

#=======
Writing out the data in a a format to be fed into the animation function
======#

H = diff(T)[1]

df = DataFrame(x1=Xs[1,:],x2=Xs[2,:],x3=Xs[3,:],x4=Xs[4,:],x5=Xs[5,:],
               x6=Xs[6,:],
               v1=Vs[1,:],v2=Vs[2,:],v3=Vs[3,:],v4=Vs[4,:],v5=Vs[5,:],
               v6=Vs[6,:],
               e1=Ws[1,:],e2=Ws[2,:],e3=Ws[3,:],e4=Ws[4,:],e5=Ws[5,:],
               e6=Ws[6,:],
               time=T,h=H)

CSV.write("6_particles.csv",df)

animate_particles(df,"6p_test.gif",frame_step=100)
