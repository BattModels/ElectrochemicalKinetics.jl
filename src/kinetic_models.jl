"""
    fermi_dirac(E, kT=0.026)

Compute the value of the Fermi-Dirac distribution at energy `E` (relative to the Fermi energy) and thermal energy `kT`.
"""
fermi_dirac(E, kT=.026) = 1 / (1 + exp(E/kT))

abstract type KineticModel end 

# checks to catch missed dispatches for new types
integrand(km::KineticModel, Reduction, V_app, direction, kT=.026) = error("A kinetic model must dispatch the `integrand` function in the reductive direction!")
integrand(km::KineticModel, Oxidation, V_app, direction, kT=.026) = error("A kinetic model must dispatch the `integrand` function in the oxidative direction!")

# dispatch for net rates
integrand(km::KineticModel, V_dl, kT=.026) = abs(integrand(km, V_dl, true; kT=kT) - integrand(km, V_dl, false, kT))

struct Marcus <: KineticModel
    λ::Float64
end

function integrand(m::Marcus, V_dl<:Real, ox::Bool, kT=.026) 
    if ox # oxidative direction
        arg = E -> -(m.λ-V_dl+E)^2 / (4*m.λ*kT)
    else # reductive direction
        arg = E -> -(m.λ+V_dl-E)^2 / (4*m.λ*kT)
    end
    E -> exp(arg(E))
end

struct MarcusHushChidsey <: KineticModel
    λ::Float64
    average_dos::Float64
end

function integrand(mhc::MarcusHushChidsey, V_dl<:Real, ox::Bool, kT=.026)
    # TODO: check dispatch of Net and add assertion if necessary
    marcus = E -> integrand(Marcus(mhc.λ), V_dl, ox, kT)
    fd(E) = ox ? 1-fermi_dirac(E, kT) : fermi_dirac(E, kT)
    return E -> mhc.average_dos * marcus(E) * fd(E)
end

struct MarcusHushChidseyDOS <: KineticModel
    λ::Float64
    dos::DOSData
end



"""
    marcus_integrand(E, λ, V_dl, ox=true; kT=.026)

Compute the value of the Marcus integrand, given by:
``-(λ - V_{dl} + E)^2 / (4 λ kT)``
in the oxidative direction (when `ox==true`), and
``-(λ + V_{dl} - E)^2 / (4 λ kT)``
in the reductive direction (when `ox==false`)

# Arguments (all energies are in eV)
* `E`: energy coordinate at which to evaluate integrand
* `λ`: reorganization energy parameter
* `V_dl`: voltage at the interface (equal to applied voltage in absence of quantum capacitance, see [`calculate_Vdl_interp`](@ref))
* `ox`: boolean flag for oxidative vs. reductive direction of reaction
* `kT`: thermal energy, defaults to room temperature value
"""
function marcus_integrand(E, λ, V_dl, ox; kT=.026)
    if ox # oxidative direction
        arg = -(λ-V_dl+E)^2 / (4*λ*kT)
    else # reductive direction
        arg = -(λ+V_dl-E)^2 / (4*λ*kT)
    end
    exp(arg)
end

"""
    MHC_integrand(E, λ, V_dl, average_dos, ox; kT=.026)
    MHC_integrand(E, λ, V_dl, dd::DOSData, ox; kT=.026)

Compute the value of the Marcus-Hush-Chidsey integrand, given by the value of the Marcus integrand (see [`marcus_integrand`](@ref)) multiplied by the appropriate Fermi-Dirac distribution (for electrons or holes).

# Arguments (all energies are in eV)
* `E`: energy coordinate at which to evaluate integrand
* `λ`: reorganization energy parameter
* `V_dl`: voltage at the interface (equal to applied voltage in absence of quantum capacitance, see [`calculate_Vdl_interp`](@ref))
* `average_dos`: average value of the density of states (this is merely a scaling factor, useful for direct comparisons to the results of the full MHC+DOS method)
* `ox`: boolean flag for oxidative vs. reductive direction of reaction
* `kT`: thermal energy, defaults to room temperature value
"""
function MHC_integrand(E, λ, V_dl, average_dos, ox; kT=.026)
    marcus = marcus_integrand(E, λ, V_dl, ox; kT=kT)
    fd = ox ? 1-fermi_dirac(E; kT=kT) : fermi_dirac(E; kT=kT)
    average_dos * marcus * fd
end

MHC_integrand(E, λ, V_dl, dd::DOSData, ox; kT=.026) = MHC_integrand(E, λ, V_dl, dd.average_value, ox; kT=kT)

"""
    MHC_integrand_net(E, λ, V_dl, average_dos; kT=.026)
    MHC_integrand_net(E, λ, V_dl, dd::DOSData; kT=.026)

Convenience function to compute the value of the net Marcus-Hush-Chidsey integrand (magnitude of difference between oxidative and reductive directions).

See also: [`MHC_integrand`](@ref)
"""
MHC_integrand_net(E, λ, V_dl, average_dos; kT=.026) = sign(V_dl) * (MHC_integrand(E, λ, V_dl, average_dos, true) - MHC_integrand(E, λ, V_dl, average_dos, false))
MHC_integrand_net(E, λ, V_dl, dd::DOSData; kT=.026) = MHC_integrand_net(E, λ, V_dl, dd.average_value; kT=kT)

"""
    MHC_DOS_integrand(E, λ, V_dl, dos_func, ox=true; kT=.026)

Compute the value of the Marcus-Hush-Chidsey integrand, given by the value of the Marcus integrand (see [`marcus_integrand`](@ref)) multiplied by the appropriate Fermi-Dirac distribution (for electrons or holes).

# Arguments (all energies are in eV)
* `E`: energy coordinate at which to evaluate integrand
* `λ`: reorganization energy parameter
* `V_dl`: voltage at the interface (equal to applied voltage in absence of quantum capacitance, see [`calculate_Vdl_interp`](@ref))
* `dos_func`: Callable object capable of evaluating the density of states at any value of energy within relevant bounds (can be produced from a data file using [`get_dos`](@ref))
* `ox`: boolean flag for oxidative vs. reductive direction of reaction
* `kT`: thermal energy, defaults to room temperature value
* `V_q`: voltage offset due to quantum capacitance [HOLDEN TO CHECK THIS DESCRIPTION]
"""
MHC_DOS_integrand(E, λ, V_dl, dos_func, ox; kT=.026, V_q=0) = dos_func(E+V_q) * MHC_integrand(E, λ, V_dl, 1, ox; kT=kT)

"""
    MHC_integrand_net(E, λ, V_dl, dos_func; kT=.026, V_q=0)

Convenience function to compute the value of the net Marcus-Hush-Chidsey + DOS integrand (magnitude of difference between oxidative and reductive directions).

See also: [`MHC_DOS_integrand`](@ref)
"""
MHC_DOS_integrand_net(E, λ, V_dl, dos_func; kT=.026, V_q=0) = dos_func(E+V_q) * MHC_integrand_net(E, λ, V_dl, 1; kT=kT)
