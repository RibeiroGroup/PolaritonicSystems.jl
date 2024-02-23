export effective_group_velocity, group_velocity, extract_LP, extract_UP


const hbar_eVps = 0.0006582119569509067 # eV⋅ps

"""
    PolaritonBranch

Data structure that organizes information from QuantumSystem by branches. 

# Fields

| Field/Attribute  |   Type   |                                            Description                                               |
|:----------------:|:---------|:-----------------------------------------------------------------------------------------------------| 
| eval             | Vector   | Vector containing the energies of the branch.                                                        |
| qvals            | Vector   | Vector containing the wavevectors of the branch.                                                     |
| mol_cont         | Vector   | Vector containing the molecular (material, dipolar, etc.) content of each eigenstate in the branch.  |
"""
struct PolaritonBranch{T}
    evals::Vector{T}
    qvals::Vector{T}
    mol_cont::Vector{T}
end

"""
    phot_cont(sys, i)

Returns the photonic content of the i-th eigenstate in the system `sys`.
"""
function phot_cont(sys::QuantumSystem, i)
    return sum(abs2.(sys.Uix[sys.phot_range, i]))
end

"""
    mol_cont(sys, i)

Returns the molecular (dipolar, material, etc) content of the i-th eigenstate in the system `sys`.
"""
function mol_cont(sys::QuantumSystem, i)
    return 1 - phot_cont(sys, i)
end

"""
    _numerical_derivative(x, y)

From an array of x values and another with y values, computes dy/dx numerically using a three point formula. 
This function is used internally to compute group velocities. 
"""
function _numerical_derivative(x, y)
    n = length(x)
    
    # Initialize the derivative array
    dy_dx = zeros(n)  

    # Use three-point formula for data in the middle
    for i in 2:n-1
        dy_dx[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
    end

    # Calculate the endpoints using a forward and a backward difference
    dy_dx[1] = (y[2] - y[1]) / (x[2] - x[1])
    dy_dx[n] = (y[n] - y[n-1]) / (x[n] - x[n-1])

    return dy_dx
end

"""
    extract_LP

Simple function that extracts the LP branch of a `sys` object returning the data as a `PolaritonBranch` object.
This object contains three arrays with energies, wavevectors, and molecular content (see the struct descriptions). 
All the arrays are sorted by the wavevector values, which allows for derivitatives to be taken more easily.

This code is used a simple tool to ultimately compute group velocities. Currently, it can only handle cases
where the number of molecules is equal the number of photonic modes.
"""
function extract_LP(sys)

    # This simple code can only handle the case where the number of photonic modes
    # is equal the number of molecules. 
    @assert length(sys.phot_wavevectors) == length(sys.mol_energies)

    N = length(sys.evals) ÷ 2

    # Get energies
    LP = sys.evals[1:N]

    # Get molecular content of each eigenstate of the LP branch
    Π = zeros(N)
    for i in eachindex(Π)
        Π[i] = mol_cont(sys, i)
    end

    # Sort the output by the wavevector (q) value
    m = sortperm(sys.phot_wavevectors)

    # Return results wrapped in the `PolaritonBranch` struct
    return PolaritonBranch(LP[m], sys.phot_wavevectors[m], Π[m])
end

"""
    extract_UP

Simple function that extracts the UP branch of a `sys` object returning the data as a `PolaritonBranch` object.
This object contains three arrays with energies, wavevectors, and molecular content (see the struct descriptions). 
All the arrays are sorted by the wavevector values, which allows for derivitatives to be taken more easily.

This code is used a simple tool to ultimately compute group velocities. Currently, it can only handle cases
where the number of molecules is equal the number of photonic modes.
"""
function extract_UP(sys)

    # This simple code can only handle the case where the number of photonic modes
    # is equal the number of molecules. 
    @assert length(sys.phot_wavevectors) == length(sys.mol_energies)

    N = length(sys.evals) ÷ 2

    # Get energies
    UP = sys.evals[(N+1):end]

    # Get molecular content of each eigenstate of the UP branch
    Π = zeros(N)
    for i in eachindex(Π)
        Π[i] = mol_cont(sys, i+N)
    end

    # Sort the output by the wavevector (q) value
    m = sortperm(sys.phot_wavevectors)

    # Return results wrapped in the `PolaritonBranch` struct
    return PolaritonBranch(UP[m], sys.phot_wavevectors[m], Π[m])
end

function group_velocity(branch::PolaritonBranch)
    # Compute group velocities as ∂E/∂q * ħ⁻¹
    return _numerical_derivative(branch.qvals, branch.evals) ./ hbar_eVps 
end

function effective_group_velocity(branch::PolaritonBranch)
    vg = group_velocity(branch)
    return branch.mol_cont .* vg
end
