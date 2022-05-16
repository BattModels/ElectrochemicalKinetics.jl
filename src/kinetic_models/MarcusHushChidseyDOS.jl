"""
    MarcusHushChidseyDOS(A=1.0, λ, dos)
    MarcusHushChidseyDOS(A=1.0, λ, dos_file)

Computes Marcus-Hush-Chidsey + DOS kinetics as described in Kurchin and Viswanathan: 10.1063/5.0023611 
"""
struct MarcusHushChidseyDOS <: IntegralModel
    A::Float64
    λ::Float64
    dos::DOSData
end

function Base.show(io::IO, mhcd::MarcusHushChidseyDOS)
    s = repr(typeof(mhcd)) * "("
    s *= "A=$(round(mhcd.A, sigdigits=3)), λ=$(round(mhcd.λ, sigdigits=3)))"
    print(io, s)
end

# default prefactor is 1
MarcusHushChidseyDOS(λ, dd::DOSData) = MarcusHushChidseyDOS(1.0, λ, dd)

MarcusHushChidseyDOS(A, λ, dos_file::Union{Matrix,String}; kwargs...) =
    MarcusHushChidseyDOS(A, λ, DOSData(dos_file; kwargs...))

MarcusHushChidseyDOS(λ, dos_file::Union{Matrix,String}; kwargs...) = MarcusHushChidseyDOS(1.0, λ, DOSData(dos_file; kwargs...))

function integrand(
    mhcd::MarcusHushChidseyDOS,
    V_dl,
    ox::Bool;
    kT::Real = 0.026,
    V_q = 0.0,
)
    function marcus_term(E)
        local exp_arg
        if ox
            exp_arg = -(( (mhcd.λ .- V_dl) .+ E) .^ 2) ./ (4 * mhcd.λ * kT)
        else
            exp_arg = -(( (mhcd.λ .+ V_dl) .- E) .^ 2) ./ (4 * mhcd.λ * kT)
        end
        exp.(exp_arg)
    end
    fd(E) = ox ? 1 .- fermi_dirac(E; kT = kT) : fermi_dirac(E; kT = kT)
    # TODO: add two dispatches as with MHC above
    E -> mhcd.A .* mhcd.dos.interp_func.(E .+ V_q) .* fd(E) .* marcus_term(E)
end