function gpdps_potts!(u,v,data,params;
                      uref=0*u,
                      vref=0*v)
    σ,τ,ω,maxit     = params
    α,γ,f,isotropic = data

    fᵧ = t -> 2t^2/(2t^2 + γ)
    if isotropic
        Fᵧ = (u₁,u₂) -> sum(@. fᵧ(sqrt(u₁^2+u₂^2)))
        κ₁ = (u₁,u₂,v₁,v₂) -> 2(1-u₁*v₁-u₂*v₂)
        κ₂ = (u₁,u₂,v₁,v₂) -> 2(1-u₁*v₁-u₂*v₂)
    else
        Fᵧ = (u₁,u₂) -> sum(@. fᵧ(u₁)+fᵧ(u₂))
        κ₁ = (u₁,u₂,v₁,v₂) -> 2(1-u₁*v₁)
        κ₂ = (u₁,u₂,v₁,v₂) -> 2(1-u₂*v₂)
    end

    # preallocate
    uold = similar(u)
    ξ    = similar(u)
    v₁,v₂ = view(v,:,:,1), view(v,:,:,2)
    u₁,u₂ = similar(u), similar(u)
    z₁,z₂ = similar(u), similar(u)
    Kᵤ,Kᵥ = similar(u), similar(v)
    for it = 1:maxit
        # Kx = ∇ᵀ(κ(∇u⋅v)⋅v)
        u₁,u₂ = ∇!(u₁,u₂,u)
        @. z₁ = κ₁(u₁,u₂,v₁,v₂)*v₁
        @. z₂ = κ₂(u₁,u₂,v₁,v₂)*v₂
        Kᵤ = ∇ᵀ!(Kᵤ,z₁,z₂)

        # print output
        if mod(it,1000) == 0
            J   = norm(u.-f)^2/(2α) + Fᵧ(u₁,u₂)
            err = norm(u.-uref)^2+norm(v.-vref)^2
            @printf("it %6d: J = %1.6e, err = %1.6e\n", it, J, err)
            imshow(u', colormap=2)
        end

        # x = prox_τG(x-τKx)
        @. uold = u
        @. u = (u-τ*Kᵤ+τ/α*f)/(1+τ/α)
        @. ξ = u + ω*(u-uold)

        # Ky = κ(∇ξ⋅v)⋅∇ξ
        u₁,u₂ = ∇!(u₁,u₂,ξ)
        @. Kᵥ[:,:,1] = κ₁(u₁,u₂,v₁,v₂)*u₁
        @. Kᵥ[:,:,2] = κ₂(u₁,u₂,v₁,v₂)*u₂

        # y = prox_σF*(y+σKy)
        @. v = (v+σ*Kᵥ)/(1+γ*σ)
    end
end

function ∇!(u₁,u₂,u)
    @. @views begin
        u₁[1:(end-1), :] = u[2:end, :] - u[1:(end-1), :]
        u₁[end, :] = 0

        u₂[:, 1:(end-1)] = u[:, 2:end] - u[:, 1:(end-1)]
        u₂[:, end] = 0
    end
    return u₁, u₂
end

function ∇ᵀ!(v,v₁,v₂)
    @. @views begin
        v[2:(end-1), :] = v₁[1:(end-2), :] - v₁[2:(end-1), :]  
        v[1, :] = -v₁[1, :]
        v[end, :] = v₁[end-1, :]

        v[:, 2:(end-1)] += v₂[:, 1:(end-2)] - v₂[:, 2:(end-1)]
        v[:, 1] += -v₂[:, 1]
        v[:, end] += v₂[:, end-1]
    end
    return v
end

function estimate_steplengths(f,α,γ,isotropic;expected_max_jump,γ_for_my)
    dyn_range = maximum(abs.(f))-minimum(abs.(f))
    expected_max_jump = 1.0
    γ_for_my = 10

    L = √8
    R_K = 2*L

    if isotropic
        mx = √2*expected_max_jump*dyn_range
    else
        mx = expected_max_jump*dyn_range
    end

    fmy = t -> 2t/(2t^2+γ_for_my)
    # fmy has its maximum in the interval [0, mx] either at mx or at √(γ_for_my/2)
    my = max(fmy(mx), fmy(√(γ_for_my/2)))

    # This estimate ignores the difficult-to-estimate ρ_y term, essentially
    # assuming ρ_y=0. To offset this, we set δ≪1.
    δ = 0.1
    μ = (1+δ)/2
    γ̃ = γ/2
    α̃ = α*2
    ξx = 1/α-1/α̃
    @assert(ξx>2*L*my^2)
    λx = 2*L^2*my^4/(ξx-2*L*my^2)
    λy = mx^2
    τ₀ = δ/λx
    τ₁ = 2*γ̃*α̃/(λy+√(λy^2+4*γ̃*α̃*(R_K^2/(1-μ)+2*λy/α̃)))
    τ  = min(τ₀, τ₁)
    σ = τ/(γ̃*α̃)
    ω = 1/(1+τ/α̃)

    return τ,σ,ω
end
