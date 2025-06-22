using DeformableGrains, DataFrames, CSV

# Example test run
N = 6
v = 1
K = [50.0,1000.0,50.0,1000.0,50.0,1000.0]
M = [1.0, 0.25, 1.0, 0.25, 1.0, 0.25]
R = [1,0.5,1,0.5,1,0.5]
time = 2

Cs, T, Xs, Vs, Ws = compression_chain(N, v, K, M, R, time; g=0, Δt=0.0001, stuck=false)

H = diff(T)[1]

df = DataFrame(x1=Xs[1,:],x2=Xs[2,:],x3=Xs[3,:],x4=Xs[4,:],x5=Xs[5,:],
               x6=Xs[6,:],
               v1=Vs[1,:],v2=Vs[2,:],v3=Vs[3,:],v4=Vs[4,:],v5=Vs[5,:],
               v6=Vs[6,:],
               e1=Ws[1,:],e2=Ws[2,:],e3=Ws[3,:],e4=Ws[4,:],e5=Ws[5,:],
               e6=Ws[6,:],
               time=T,h=H)

CSV.write("6_particles.csv",df)

#animate_particles(df,"test.gif",frame_step=100)
