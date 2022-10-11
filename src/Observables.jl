"""
    mean_square_disp

Computes the average squared displacement (with respect to an initial position μ): ⟨(x - x0)²⟩

# Arguments

| Argument |       Type      |                                       Description                                                                  |
|:--------:|:----------------|:-------------------------------------------------------------------------------------------------------------------| 
| state    | Vector          | Vector representing the state to be analyzed.                                                                      |
| x0       | Number/Quantity | Position for which the displacement is computed against. If a number, unit is taken as `nm`.                       |
| system   | QuantumWire     | QuantumWire object - constaints all static information of the system.                                              |
"""
function mean_square_disp(state::Vector, x0::Number, system::QuantumWire)

    # Transform to localized (uncoupled) state
    unc_state = system.Uix * state

    # Loop through molecules
    x2 = 0.0
    n = 1
    for i = system.mol_range
        # Contribution of each molecule is P * (x-x0)^2, where P is the probability of the state
        # being localized at that molecule and x is its position
        x2 += abs2(unc_state[i]) * (system.mol_positions[n] - x0)^2
        n += 1
    end

    # Renormalize (conditional probability) with the probability of finding the molecular exciton
    Pmol = prob_any_mol(state, system)
    return x2 / Pmol
end

# Support for units
function mean_square_disp(state::Vector, x0::Quantity, system::QuantumWire)
    mean_square_disp(state, ustrip(u"nm", x0), system)
end

"""
    prob_any_phot

For a given state vector, computes the probability of finding a photon in any mode.

# Arguments

| Argument |    Type     |                                       Description                                                                  |
|:--------:|:------------|:-------------------------------------------------------------------------------------------------------------------| 
| state    | Vector      | Vector representing the state to be analyzed.                                                                      |
| system   | QuantumWire | QuantumWire object - constaints all static information of the system.                                              |
"""
function prob_any_phot(state::Vector, system::QuantumWire)
    locstate = system.Uix * state 
    return abs2.(locstate[system.phot_range]) |> sum
end

"""
    prob_any_mol

For a given state vector, computes the probability of finding an exciton on any molecule.

# Arguments

| Argument |    Type     |                                       Description                                                                  |
|:--------:|:------------|:-------------------------------------------------------------------------------------------------------------------| 
| state    | Vector      | Vector representing the state to be analyzed.                                                                      |
| system   | QuantumWire | QuantumWire object - constaints all static information of the system.                                              |
"""
function prob_any_mol(state::Vector, system::QuantumWire)
    locstate = system.Uix * state 
    return abs2.(locstate[system.mol_range]) |> sum
end

"""
    get_exciton_prob

For a given state vector, returns a renormalized vector with probabilities of finding the excitation energy on each molecule.

# Arguments

| Argument |    Type     |                                       Description                                                                  |
|:--------:|:------------|:-------------------------------------------------------------------------------------------------------------------| 
| state    | Vector      | Vector representing the state to be analyzed.                                                                      |
| system   | QuantumWire | QuantumWire object - constaints all static information of the system.                                              |
"""
function get_exciton_prob(state::Vector, system::QuantumWire)
    loc = system.Uix * state
    loc = abs2.(loc[system.mol_range]) 

    # Renormalize 
    return loc ./ sum(loc)
end