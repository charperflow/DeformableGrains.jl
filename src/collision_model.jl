#collision_model_2_f.jl

function compression_chain(N,v,K,M,R,time;Δt = 0.0001,g=9.81,γ=0,stuck=true)   
    h = convert(Int64,ceil(time/Δt))+1 #plus one to allow for start at zero
    T       = [(i-1)*Δt for i in 1:h]
    Xs, Vs  = setup_dust_settled(N,v,K,M,R,h;g=g)
    Ws      = Matrix{Float64}(undef,N,h)
    Cs      = []
    push!(Cs,[])
    for i in 1:h
        Ws[:,i] = 2*R
    end
    
    
    #println("===============================")
    #println("===============================")
    for i in 2:h
    
        prcnt_done = round((i/h)*100,digits=2)
        time_done = T[i]
        #println("================")
        
        if mod(i,1000)==0
            print("We are %$prcnt_done done \r")
        end
        
    
       # println("We are $time_done seconds")
       # println("================")
        #sleep(0.5)
        
    
        X = Xs[:,i-1]
        V = Vs[:,i-1]
        W = Ws[:,i-1]
        #println("into EYYEE")
        I   = [X;V]
        k1  = F(I,K,M,R,W,g)
        b   = I+(Δt/2)*k1
        k2  = F(b,K,M,R,W,g)
        c   = I+(Δt/2)*k2
        k3  = F(c,K,M,R,W,g)
        d   = I+(Δt*k3)
        k4  = F(d,K,M,R,W,g)
    
        I_new = I+(Δt/6)*(k1+2*k2+2*k3+k4)
        X_new = I_new[1:N]
        #X_new[N] = (N-1)*(2*r_0)
        V_new = I_new[N+1:2*N]
        if stuck
            V_new[N] = 0 
        end
        W_new = [2*R[i] for i in 1:N]
        #println("about to dive")
        # inside your time‐loop, after you compute Xs[:, t] and before you call find_chains
        #=
        if i == 2863-1 || i == 2863
            @info "DEBUG SNAPSHOT t=$i" X = X_new
            chains, B_dict, W = find_chains(X_new, R, K)
            @info "  chains=" chains
            @info "  widths=" W
        end
        =#
        chains,Bs, W_new = find_chains(X_new,R,K)
    
        push!(Cs,chains)
       #=
        for i in 1:N
            W_new[i] = (W_new[i]>1 ? 1 : W_new[i])
        end
        =#
        Xs[:,i]=X_new
        Vs[:,i]=V_new
        Ws[:,i]=W_new

    end
    
    #println("===============================")
    #println("===============================")    
    
    H = [Δt for i in 1:h]

    return Cs, T, Xs, Vs, Ws
    
end #function
