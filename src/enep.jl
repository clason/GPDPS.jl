function gpdps_enep!(u,v,data,params)
    σ,τ,ω,maxit = params
    S,Sᵀ,z₁,z₂,α,a,b,ue₁,ue₂ = data 

    # preallocate
    uold = similar(u)
    ξ = similar(u)
    u₁,u₂ = view(u,:,1), view(u,:,2)
    v₁,v₂ = view(v,:,1), view(v,:,2)
    ξ₁,ξ₂ = view(ξ,:,1), view(ξ,:,2)
    ρ₁,ρ₂ = similar(u₁), similar(v₁)
    Kᵤ,Kᵥ = similar(u), similar(v)
    for it = 1:maxit
        # Kx
        yᵥᵤ = S(v₁,u₂)
        yᵤᵥ = S(u₁,v₂)
        yᵤᵤ = S(u₁,u₂)
        @. ρ₁ = 2yᵤᵤ - yᵤᵥ - z₁
        @. ρ₂ = 2yᵤᵤ - yᵥᵤ - z₂
        Kᵤ .= Sᵀ(ρ₁,ρ₂) .+ α.*u

        # print output
        Psi = (norm(yᵥᵤ-z₁)^2/2 + α/2*norm(v₁) - norm(yᵤᵤ-z₁)^2/2 - α/2*norm(u₁) +
               norm(yᵤᵥ-z₂)^2/2 + α/2*norm(v₂) - norm(yᵤᵤ-z₂)^2/2 - α/2*norm(u₂)) /
              length(z₁)
        err = (norm(u₁-ue₁)^2 + norm(u₂-ue₂)^2)/length(z₁)
        @printf("it %3d: Psi = %1.3e\t Err = %1.3e\n", it, Psi, err)
        
        # x = prox_G(x-τKx)
        @. uold = u
        @. u = clamp(u-τ*Kᵤ,a,b)
        @. ξ = u + ω*(u-uold)
    
        # Ky
        yᵥᵤ = S(v₁,ξ₂)
        yᵤᵥ = S(ξ₁,v₂)
        @. ρ₁ = z₁ - yᵤᵥ
        @. ρ₂ = z₂ - yᵥᵤ
        Kᵥ .= Sᵀ(ρ₁,ρ₂) .- α.*v

        # y = prox_F*(y+σKy)
        @. v = clamp(v+σ*Kᵥ,a,b)
    end
end

function setup_enep(N,α,a,b)
    x₁,x₂,A = setup_fdm(N)

    B₁ = x₂.<0.5
    B₂ = x₂.>0.5

    y  = @. sin(π*x₁)*sin(π*x₂)
    p₁ = @. sin(2π*x₁)*sin(2π*x₂)
    p₂ = @. sin(3π*x₁)*sin(3π*x₂)

    z₁ = y .+ A*p₁
    z₂ = y .+ A*p₂

    u₁ = clamp.(B₁.*p₁./α,a,b)
    u₂ = clamp.(B₂.*p₂./α,a,b)

    f = A*y .- B₁.*u₁ .- B₂.*u₂

    C  = cholesky(A)
    S  = (u₁,u₂) -> C\(B₁.*u₁ .+ B₂.*u₂ .+ f)
    Sᵀ = (q₁,q₂) -> hcat(B₁.*(C\q₁),B₂.*(C\q₂))

    return S,Sᵀ,z₁,z₂,u₁,u₂
end

function setup_fdm(N)
    dx = 1/(N+1)
    xm = range(0,stop=1,length=N)
    xx,yy = reshape(xm, 1, N), reshape(xm, N, 1)
    x₁,x₂ = repeat(xx,outer=(N,1))[:], repeat(yy,outer=(1,N))[:]

    e  = fill(1/dx^2,N)
    D₂ = spdiagm(-1=>-e[2:end],0=>2e,1=>-e[2:end])
    Id = sparse(I,N,N)
    A  = kron(D₂,Id)+kron(Id,D₂)
    return x₁,x₂,A
end

