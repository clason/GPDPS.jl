module GPDPS

using Printf
using FileIO, JLD2
using ImageMagick, GR
using LinearAlgebra, SparseArrays

export test_enep
export test_potts

include("enep.jl")
include("potts.jl")

function test_enep(N=128)
    α = 1
    a = -0.5
    b = 0.5

    S,Sᵀ,z₁,z₂,ue₁,ue₂ = setup_enep(N,α,a,b)
    data = S,Sᵀ,z₁,z₂,α,a,b,ue₁,ue₂

    σ = 1.0
    τ = 0.99
    ω = 1.0
    maxit = 10
    params = σ,τ,ω,maxit

    u = zeros(N*N,2)
    v = zeros(N*N,2)
    gpdps_enep!(u,v,data,params)
    return u,v
end

function test_potts(α=1, γ=1e-3;
                    image     = "blobs",
                    isotropic = true,
                    maxit     = 100000)

    # load test image
    imgfile = joinpath(dirname(pathof(GPDPS)), "..", "img", image*".tif")
    f = Float64.(load(imgfile))
    data = α,γ,f,isotropic

    # compute step size estimates
    τ,σ,ω = estimate_steplengths(f,α,γ,isotropic,expected_max_jump=1.0,γ_for_my=10)
    @printf("Step length parameters: τ = %g, σ = %g, ω = %f\n", τ, σ, ω)
    params = σ,τ,ω,maxit

    # starting points
    u = copy(f)
    v = zeros((size(f)...,2))

    # load reference file to compute errors
    h = string(hash((isotropic,α,γ)));
    ref_file = "ref_" * image * "_" * h * ".jld"
    if isfile(ref_file)
        println("Loading reference file "*ref_file)
        @load ref_file uref vref
    else
        println("No reference file, using uref=(0,0)")
        uref = zero(u)
        vref = zero(v)
    end

    gpdps_potts!(u,v,data,params, uref=uref, vref=vref)

    # save reference file to compute errors
    println("Saving reference file " * ref_file)
    uref = u; vref = v;
    @save ref_file uref vref

    return u,v
end

end
