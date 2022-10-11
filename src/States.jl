"""
    time_propagate

Returns a new vector with the given state evolved in time. This is free evolution, i.e., the Hamiltonian associated with the
system is assumed to be time-independent. 

# Arguments

| Argument  |       Type      |                                       Description                                                                  |
|:---------:|:----------------|:-------------------------------------------------------------------------------------------------------------------| 
| state     | Vector          | Vector representing the state to be evolved, it must be expressed in the eigenbasis of the system.                 |
| system    | QuantumWire     | QuantumWire object - constaints all static information of the system.                                              |
| t         | Number/Quantity | Time for which the system must be evolved into. If a number, unit is taken as `ps`.                                |
"""
function time_propagate(state::Vector, system::QuantumWire, t::Number)

    ħ = ustrip(u"eV*ps", CODATA2018.PlanckConstant) / 2π

    expv = [exp(-im*E*t/ħ) for E in system.evals]

    # newstate = state * exp(-iEt/h)
    return expv .* state
end

function time_propagate(state::Vector, system::SQuantumWire, t::Number)

    ħ = ustrip(u"eV*ps", CODATA2018.PlanckConstant) / 2π

    # Vector to accumulate terms of the Taylor series
    accum = similar(state)

    # Vector to hold a specific term of the Taylor series
    term = similar(state)

    n = 1
    # Common factor 
    fac = -im*t/ħ
    U = -(im*t/ħ) .* system.H

    # Zeroth term is just the initial vector
    accum .= state

    # First term of the series -iδt/ħ H|s⟩
    #term = fac * (system.H * state)
    term = U*state
    accum += term

    println("Stating time propagation...")
    while √(sum(abs2.(term))) > 1e-10
        n += 1
        if n > 100
            break
        end

        term = fac/n * (system.H * term)
        #term = term .- dot(term, accum)
        #term = U^n * state / factorial(big(n))

        accum += term
        #normalize!(accum)
        println("Ite $(n) term residue $(√(sum(abs2.(term))))")
    end
    return accum
end

# Support for units
function time_propagate(state::Vector, system::QuantumSystem, t::Quantity) 
    time_propagate(state, system, ustrip(u"ps", t))
end

"""
    create_exciton_wavepacket

Creates a localized excited state as Gaussian wavepacket centered at μ with standard deviation of σ.
It represents the probability of finding the excitation energy in a specific molecule, which in turn has a well defined position.

# Arguments

| Argument |       Type      |                                       Description                                                                  |
|:--------:|:----------------|:-------------------------------------------------------------------------------------------------------------------| 
| μ        | Number/Quantity | Average position of the wave packet. If a number, unit is taken as `nm`.                                           |
| σ        | Number/Quantity | Standard deviation for the Gaussian wave packet. If a number, unit is taken as `nm`.                               |
| system   | QuantumWire     | QuantumWire object - constaints all static information of the system.                                              |
"""
function create_exciton_wavepacket(μ::Number, σ::Number, system::QuantumSystem)

    N = length(system.mol_energies) + length(system.phot_energies)

    locstate = zeros(ComplexF64, N)

    # Compute expansion coefficients for the molecular states as √P
    # where P is a normal Gaussian distribution
    for (i, xn) in zip(system.mol_range, system.mol_positions)
        locstate[i] = sqrt.(1/√(2π*σ^2) * exp(- (xn - μ)^2 / 2σ^2))
    end

    # Normalize
    normalize!(locstate)

    if typeof(system) <: QuantumWire
        # Convert localized (uncoupled) basis to eigenbasis
        return system.Uix' * locstate
    end
    
    # If SQuantumState, return vector in the local basis
    return locstate
end

# Support for units
function create_exciton_wavepacket(μ::Quantity, σ::Quantity, system::QuantumSystem) 
    create_exciton_wavepacket(ustrip(u"nm", μ), ustrip(u"nm", σ), system)
end

# Support for units
function create_exciton_wavepacket(μ::Quantity, σ::Number, system::QuantumSystem) 
    create_exciton_wavepacket(ustrip(u"nm", μ), ustrip(u"nm", σ * unit(μ)), system)
end