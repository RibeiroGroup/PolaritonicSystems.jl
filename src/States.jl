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
    out = similar(state)
    time_propagate!(out, state, system, t)
    return out
end

function time_propagate!(out::Vector, state::Vector, system::QuantumWire, t::Number)
    Threads.@threads for i in eachindex(out)
        out[i] = exp(imħ*system.evals[i]*t) * state[i]
    end
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
| q        | Number/Quantity | KEYWORD ARGUMENT - Average wavenumber of the exciton - adds complex component.                                     |
"""
function create_exciton_wavepacket(μ::Number, σ::Number, system::QuantumSystem; q::Number=0.0)

    N = length(system.mol_energies) + length(system.phot_energies)

    locstate = zeros(ComplexF64, N)

    # Compute expansion coefficients for the molecular states as √P
    # where P is a normal Gaussian distribution
    for (i, xn) in zip(system.mol_range, system.mol_positions)
        locstate[i] = sqrt.(1/√(2π*σ^2) * exp(- (xn - μ)^2 / 2σ^2)) * exp(im*q*xn)
    end

    # Normalize
    normalize!(locstate)

    # Convert localized (uncoupled) basis to eigenbasis
    return system.Uix' * locstate
end

# Support for units
function create_exciton_wavepacket(μ::Quantity, σ::Quantity, system::QuantumSystem; q::Quantity=0.0nm^-1) 
    create_exciton_wavepacket(ustrip(u"nm", μ), ustrip(u"nm", σ), system, q=ustrip(u"nm^-1", q))
end

# Support for units
function create_exciton_wavepacket(μ::Quantity, σ::Number, system::QuantumSystem; q::Quantity=0.0nm^-1)  
    create_exciton_wavepacket(ustrip(u"nm", μ), ustrip(u"nm", σ * unit(μ)), system, q=ustrip(u"nm^-1", q))
end

function energy_filter!(state::Vector, sys::QuantumSystem, ω0, δω)
    for i in eachindex(state)
        ωi = sys.evals[i]
        state[i] *= exp(-(ωi-ω0)^2 / (2*δω^2))
    end
    normalize!(state)
end

"""
    get_q_transformation

From a system object (`sys`) returns a transformation matrix that takes a state vector in the *local* basis
and changes into the k (or q) representation, that is, reciprocal space. 

!! Not thouroughly tested for when the number of photon modes is not equal the number of sites. 
"""
function get_q_transformation(sys::QuantumSystem)
    Nm = length(sys.mol_range)
    Nc = length(sys.phot_range)
    Um = zeros(ComplexF64, Nc, Nm)

    for i = 1:Nc
        q = sys.phot_wavevectors[i]
        for j = 1:Nm
            r = sys.mol_positions[j]
            Um[i,j] = exp(-im*q*r)
        end
    end

    U = zeros(ComplexF64, Nc+Nm, Nc+Nm)

    U[1:Nc, 1:Nc] .= I(Nc)

    U[(Nc+1):end, (Nc+1):end] .= Um ./ √(Nm)

    return U
end
