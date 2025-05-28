module HallPIC

using DataInterpolations: DataInterpolations, LinearInterpolation
using DocStringExtensions
using Random: Random
using DelimitedFiles 


#===================================================
Nondimensionalization notes

# Base quantities
m_0 = 1 amu 		[mass]
q_0 = q_e 			[charge]
n_0 = 1e14			[volume -> length]
phi_0 = q T_0 = q E_0 = 1 V or 1 eV  [energy -> time]

x_0 = 1 / (n_0^1/3)
u_0 = sqrt(q_e * phi_0 / m_0)
t_0 = x_0 / u_0
E_0 = phi_0 / x_0
===================================================#

# Physical constants and base units
const q_e = 1.608e-19
const N_A = 6.022e26
const m_0 = 1 / N_A
const n_0 = 1e12
const phi_0 = 1.0
const m_e = 0.00054858 #normalized electron mass 


# Derived quantities
const x_0 = 1 / cbrt(n_0)
const u_0 = sqrt(q_e * phi_0 / m_0)
const t_0 = x_0 / u_0
const E_0 = phi_0 / x_0

#======================================================
Gas and Species definitions
======================================================#

@kwdef struct Gas
	mass::Float64
	name::Symbol
end

"""
$(TYPEDEF)
A single chemical species (charge + mass)
$(TYPEDFIELDS)
"""
@kwdef struct Species
	gas::Gas
    reactions::Vector{UInt8}
	charge::Int8
end

@inline (g::Gas)(Z::Integer) = Species(gas=g, charge=Int8(Z), reactions=UInt8[])

# Helper methods
@inline mass(s::Species) = s.gas.mass
@inline charge(s::Species) = s.charge
@inline charge_to_mass(s::Species) = s.charge / s.gas.mass

#======================================================
ParticleContainer definitions
======================================================#

"""
$(TYPEDEF)
Per-particle information for a single species.
`T` is expected to be either Float32 or Float64
$(TYPEDFIELDS)
"""
struct ParticleContainer{T<:AbstractFloat}
    weight::Vector{T}
    pos::Vector{T}
    vel::Vector{T}
    acc::Vector{T}
    inds::Vector{Int64}
	species::Species
    n_d::Int32
end

function ParticleContainer{T}(N, species::Species) where T<:AbstractFloat
	weight = zeros(T, N)
	pos = zeros(T, N)
	vel = zeros(T, N)
	acc = zeros(T, N)
    inds = zeros(Int64, N)
    n_d = 500 # hard code for now, should be a simulation hyperparamter
	return ParticleContainer{T}(
		weight, pos, vel, acc, inds, species, n_d 
	)
end

Base.length(pc::ParticleContainer) = length(pc.pos)
Base.firstindex(pc::ParticleContainer) = firstindex(pc.pos)
Base.lastindex(pc::ParticleContainer) = lastindex(pc.pos)
Base.eachindex(pc::ParticleContainer) = firstindex(pc):lastindex(pc)

"""
$(TYPEDSIGNATURES)

Resize `pc` to contain `n` elements.
If `n` is smaller than the current length, the first `n` elements will be retained.
If `n` is larger, the new particles are not guaranteed to be initialized.
"""
function Base.resize!(pc::ParticleContainer, n::Integer)
    resize!(pc.weight, n)
    resize!(pc.pos, n)
    resize!(pc.vel, n)
    resize!(pc.acc, n)
    resize!(pc.inds, n)
end

@inline function swap!(v::AbstractArray, ind_1::Integer, ind_2::Integer)
    v[ind_1], v[ind_2] = v[ind_2], v[ind_1]
    return v
end

"""
$(TYPEDSIGNATURES)

Swap the particle in position `ind_1` with that in position `ind_2`.
"""
function swap!(pc::ParticleContainer, ind_1::Integer, ind_2::Integer)
    swap!(pc.weight, ind_1, ind_2)
    swap!(pc.pos, ind_1, ind_2)
    swap!(pc.vel, ind_1, ind_2)
    swap!(pc.acc, ind_1, ind_2)
    swap!(pc.inds, ind_1, ind_2)
end

"""
$(TYPEDSIGNATURES)

Returns `true` if the particle at position `i` should be removed
"""
@inline should_remove(pc::ParticleContainer, i::Integer) = pc.weight[i] <= 0

"""
$(TYPEDSIGNATURES)

Partition a collection based on a predicate.
At completion, elements for which `pred(v, index)` returns true will be at the back of the array.
Requires a function `swap(v, index_1, index_2)` that swaps the elements of `v` at the given indices.
Returns the index of the last `false` element.

Requires that `eachindex`, `firstindex`, and `lastindex` are implemented for the given collection.
"""
function partition!(v, pred=getindex, swap=swap!)
    i = firstindex(v)
    j = lastindex(v)
    while i <= lastindex(v)
        # Increment `i` until we find a `true` element
        if pred(v, i)
            # Decrement `j` until we find a `false` element
            while j > i && pred(v, j)
                j -= 1
            end

            # If `i == j`, we're done and can return
            if j <= i
                break
            end

            swap(v, i, j)
        end
        i += 1
    end

    return i - 1
end

"""
$(TYPEDSIGNATURES)

Remove particles with zero weight.
"""
function remove_flagged_particles!(pc::ParticleContainer{T}) where T
    new_size = partition!(pc, should_remove, swap!)
    resize!(pc, new_size)
    return pc
end


"""
$(TYPEDSIGNATURES)
Add particles to a `ParticleContainer`.
"""
function add_particles!(pc::ParticleContainer{T}, x::Vector{T}, v::Vector{T}, w::Vector{T}, index::Vector{Int64}) where T
	# check that new arrays have same length
	M = length(x)
	N = length(pc)
	@assert M == length(v) && M == length(w)
	# append position, velocity, weight to pc arrays
	append!(pc.pos, x)
	append!(pc.vel, v)
	append!(pc.weight, w)
    append!(pc.inds, index)
	# extend acceleration array to correct size and fill with zeros
	resize!(pc.acc, N+M)
	for i in 1:M
		pc.acc[N+i] = zero(T)
	end
	return pc
end

function push_pos!(pc::ParticleContainer, dt::AbstractFloat)
    @inbounds for i in eachindex(pc.pos)
        pc.pos[i] = muladd(pc.vel[i], dt, pc.pos[i])
    end
	return pc
end

function push_vel!(pc::ParticleContainer, dt::AbstractFloat)
    @inbounds for i in eachindex(pc.vel)
        pc.vel[i] = muladd(pc.acc[i], dt, pc.vel[i])
    end
	return pc
end


"""
$(TYPEDEF)

Particles leave the domain through this boundary.
"""

struct OpenBoundary 

end





"""
$(TYPEDEF)

Particles hit a wall and accommodate to the given temperature.
Charged particles neutralize when they hit the wall.

$(TYPEDFIELDS)
"""
mutable struct WallBoundary
    "The temperature of the wall, non-dimensionalized"
    temperature::Float32
    "The number of particles of a single species that have hit the wall this timestep"
    weight::Float32
    function WallBoundary(temperature::AbstractFloat)
        return new(Float32(temperature), Float32(0.0))
    end
    
end


const HeavySpeciesBoundary = Union{OpenBoundary, WallBoundary}

@enum BoundaryLoc begin
    Left
    Right
end


"""
$(TYPEDSIGNATURES)

Flag particles in the given cell for deletion by setting their weight to zero.
The `remove_flagged_particles!` function can then be called to remove particles with zero weight.
"""
function flag_particles_in_cell!(pc::ParticleContainer{T}, cell_index) where T
    for (i, ic) in enumerate(pc.inds)
        pred = abs(ic) != cell_index
        pc.weight[i] *= pred
    end
end


#======================================================
SpeciesProperties and CellProperties definitions
======================================================#

"""
$(TYPEDEF)
T: Float32 or Float64
Holds fluid properties (density, velocity, temperature) for a single species.
$(TYPEDFIELDS)
"""
struct SpeciesProperties{T<:AbstractFloat}
    dens::Vector{T}
    vel::Vector{T}
    temp::Vector{T}
    avg_weight::Vector{T}
    N_particles::Vector{UInt64}
	species::Species
end

"""
$(TYPEDEF)
Contains fluid information for a single cell.
$(TYPEDFIELDS)
"""
struct CellProperties{T<:AbstractFloat}
	dens::T
	vel::T
	temp::T
	avg_weight::T
    N_particles::UInt64
end

function SpeciesProperties{T}(num_cells, species::Species) where T<:AbstractFloat
	dens = zeros(T, num_cells)
	vel = zeros(T, num_cells)
	temp = zeros(T, num_cells)
	avg_weight = zeros(T, num_cells)
    N_particles = zeros(UInt64, num_cells)
	return SpeciesProperties{T}(dens, vel, temp, avg_weight, N_particles, species)
end


#======================================================
  Iterator interface definition for SpeciesProperties
======================================================#

Base.length(sp::SpeciesProperties) = length(sp.dens)

function Base.getindex(sp::SpeciesProperties{T}, i) where T
	return CellProperties{T}(
		sp.dens[i], sp.vel[i], sp.temp[i], sp.avg_weight[i], sp.N_particles[i]
	)
end

Base.isdone(sp::SpeciesProperties, i) = i > length(sp)
Base.eltype(::SpeciesProperties{T}) where T = CellProperties{T}

Base.iterate(sp::SpeciesProperties) = length(sp) > 0 ? (sp[1], 2) : nothing 
function Base.iterate(sp::SpeciesProperties, i)
	i > length(sp) && return nothing
	return sp[i], i+1
end

#======================================================
Geometry
======================================================#

"""
$(TYPEDSIGNATURES)
Create a particle container and fill it with particles matching the macroscopic properties
from a provided SpeciesProperties object.

Accounting note: the edges of cell i are edges[i] and edges[i+1]
"""
function initialize_particles(sp::SpeciesProperties{T}, grid, particles_per_cell) where T
	pos_buf = zeros(T, particles_per_cell)
	vel_buf = zeros(T, particles_per_cell)
	weight_buf = zeros(T, particles_per_cell)
    inds_buf = zeros(Int64, particles_per_cell)
	@assert length(grid.cell_centers) == length(sp)

	pc = ParticleContainer{T}(0, sp.species)

    n_itp = LinearInterpolation(sp.dens, grid.cell_centers)
    u_itp = LinearInterpolation(sp.vel, grid.cell_centers)
    T_itp = LinearInterpolation(sp.temp, grid.cell_centers)

    for i in 2:length(grid.cell_centers) - 1
		z_L = grid.face_centers[i]
		z_R = grid.face_centers[i+1]
		dz = z_R - z_L
        dz_part = dz / particles_per_cell

        # Load particles in uniformly and evenly-spaced
        # In addition to being smoother than random particles, it also causes the particles to be sorted by position.
		pos_buf .= range(z_L + dz_part/2, z_R - dz_part/2, length=particles_per_cell)

        inds_buf = copysign.(i, grid.cell_centers[i] .- pos_buf)

		Random.randn!(vel_buf)
		@. vel_buf *= sqrt(T_itp(pos_buf) / pc.species.gas.mass)
		@. vel_buf = vel_buf + u_itp(pos_buf)

        V = grid.cell_volumes[i]
		@. weight_buf = n_itp(pos_buf) / particles_per_cell * V

		add_particles!(pc, pos_buf, vel_buf, weight_buf, inds_buf)
	end

	return pc
end


struct Grid
    dz::Float64
    cell_centers::Vector{Float64}
    cell_volumes::Vector{Float64}
    face_centers::Vector{Float64}
    face_areas::Vector{Float64}
    left_boundary::HeavySpeciesBoundary
    right_boundary::HeavySpeciesBoundary
end

function Grid(num_cells::Integer, left, right, area, left_boundary = WallBoundary(0.5), right_boundary = WallBoundary(0.5))
    dz = (right - left) / num_cells

    face_centers = collect(range(left - dz, right + dz, step = dz))
    face_areas = fill(area, length(face_centers))

    cell_centers = zeros(num_cells+2)
    cell_volumes = zeros(num_cells+2)

    for i in eachindex(cell_centers)
        z_L, z_R = face_centers[i], face_centers[i+1]
        cell_centers[i] = 0.5 * (z_L + z_R)
        cell_volumes[i] = area * (z_R - z_L)

    end

    return Grid(dz, cell_centers, cell_volumes, face_centers, face_areas, left_boundary, right_boundary)
end

"""
$(TYPEDSIGNATURES)
"""
function locate_particles!(pc::ParticleContainer{T}, grid::Grid) where T
    @inbounds for (i, x) in enumerate(pc.pos)
        
        cell_index = floor(Int64, x / grid.dz) + 2
        center_pos = grid.cell_centers[cell_index]
        pc.inds[i] = copysign(cell_index, x-center_pos)
    end

    return pc.inds
end


function find_cell_indices(cell_index, grid::Grid)

    # don't deposit into ghost cells
    if cell_index == -2 || cell_index == length(grid.cell_centers)-1
        s, ic = 0, abs(cell_index)
    else
        s, ic = sign(cell_index), abs(cell_index)
    end

    return s, ic
end


function apply_boundary!(pc::ParticleContainer, grid::Grid, ::OpenBoundary, flag::Int8)
    # flagging particles, then call remove flagged particles at the end of every step  
    if flag == -1
        for (ip, pos) in enumerate(pc.pos)
            if (pos - grid.face_centers[2])/grid.dz <= -0.5  
                pc.weight[ip] = 0.0
            end
        end
    else
        for (ip, pos) in enumerate(pc.pos)
            if (pos - grid.face_centers[end-1])/grid.dz >= 0.5  
                pc.weight[ip] = 0.0
            end
        end
    end
    return pc
end


function apply_boundary!(pc::ParticleContainer, grid::Grid, ::WallBoundary, flag::Int8)
    """
    if flag == -1
        for (ip, pos) in enumerate(pc.pos)
            if (pos - grid.face_centers[2])/grid.dz <= -0.5  
                pc.vel[ip] *= -1.0
            end
        end
    else
        for (ip, pos) in enumerate(pc.pos)
            if (pos - grid.face_centers[end-1])/grid.dz >= 0.5  
                pc.vel[ip] *= -1.0
            end
        end
    end

    """
    return pc
end


"""
$(TYPEDSIGNATURES)
Push particle container to next timestep using Leapfrog scheme
"""
function push!(pc::ParticleContainer, dt::AbstractFloat, grid::Grid)
    push_vel!(pc, dt)
    push_pos!(pc, dt)

    # enforce boundary conditions 
    apply_boundary!(pc, grid, grid.left_boundary, Int8(-1))
    apply_boundary!(pc, grid, grid.right_boundary, Int8(1))

    # cleanup 
    remove_flagged_particles!(pc)

    return pc
end

"""
$(TYPEDSIGNATURES)

Deposit particle properties into a gridded SpeciesProperties object
"""
function deposit!(fluid_properties::SpeciesProperties{T}, particles::ParticleContainer, grid::Grid, avg_interval=50) where T

    # zero density, velocity, temperature, particle count
    for i in eachindex(fluid_properties.dens)
        fluid_properties.dens[i] = 0
        fluid_properties.vel[i] = 0
        fluid_properties.temp[i] = 0
        fluid_properties.N_particles[i] = 0
    end
    # loop over particles, deposit density and momentum density
    for (ip, z_part) in enumerate(particles.pos)

        # The cell index is positive if the particle is right of the cell center,
        # or negative if it is left of the cell center
        s, ic = find_cell_indices(particles.inds[ip], grid)

        # Grid information and cell weighting
        inv_vol_c = 1 / grid.cell_volumes[ic]
        inv_vol_s = 1 / grid.cell_volumes[ic+s]
        z_cell = grid.cell_centers[ic]
        t = abs((z_part - z_cell) / grid.dz)

        # Particle velocity and weight
        up = particles.vel[ip]
        wp = particles.weight[ip]

        # Density contribution to each cell
        dens_c = (1 - t) * wp * inv_vol_c
        dens_s = t * wp * inv_vol_s
        fluid_properties.dens[ic]   += dens_c 
        fluid_properties.dens[ic+s] += dens_s

        # Momentum contribution to each cell
        fluid_properties.vel[ic]    += dens_c * up
        fluid_properties.vel[ic+s]  += dens_s * up

        # add to the counter of the number of particles in the cells this particle touches 
        fluid_properties.N_particles[ic]   += 1
        fluid_properties.N_particles[ic+s] += 1
    end

    # Calculate cell-centered fluid properties
    for ic in 2:length(grid.cell_centers)-1 
        density = fluid_properties.dens[ic]
    
        # Get avg. velocity from momentum density
        fluid_properties.vel[ic] /= density

        # Compute smoothed average weight in cell
        new_avg_wt = density / fluid_properties.N_particles[ic] * grid.cell_volumes[ic]
        old_avg_wt = fluid_properties.avg_weight[ic]
        fluid_properties.avg_weight[ic] = ((avg_interval - 1) * old_avg_wt + new_avg_wt) / avg_interval
    end

    # Zero temperature again
    fluid_properties.temp .= 0.0

    # Calculate energy density
    for (ip, z_part) in enumerate(particles.pos)

        # The cell index is positive if the particle is right of the cell center,
        # or negative if it is left of the cell center
        s, ic = find_cell_indices(particles.inds[ip], grid)

        # Grid information and cell weighting
        inv_vol_c = 1.0 / grid.cell_volumes[ic]
        inv_vol_s = 1.0 / grid.cell_volumes[ic+s]
        z_cell = grid.cell_centers[ic]
        t = abs((z_part - z_cell) / grid.dz)

        # Particle velocity and weight
        up = particles.vel[ip]
        wp = particles.weight[ip]

        # Cell average velocity 
        vel_c = fluid_properties.vel[ic]
        vel_s = fluid_properties.vel[ic+s]

        # Accumulate energy density
        dens_c = (1 - t) * wp * inv_vol_c
        dens_s = t * wp * inv_vol_s
        fluid_properties.temp[ic] += dens_c * (up - vel_c)^2
        fluid_properties.temp[ic+s] += dens_s * (up - vel_s)^2
    end

    m = fluid_properties.species.gas.mass

    # Convert energy density to nondimensional temperature
    for ic in 2:length(grid.cell_centers)-1 
        fluid_properties.temp[ic] *= m / (fluid_properties.dens[ic])
    end

    # Linearly extrapolate bulk properties to ghost cells
    # enforce maximum of order of magnitude drop to prevent negative density 
    n2, n3 = fluid_properties.dens[2], fluid_properties.dens[3]
    n1 = maximum([2 * n2 - n3, n2/10])
    dens_ratio = n2 / n1
    fluid_properties.dens[1] = n1
    fluid_properties.vel[1] = fluid_properties.vel[2] * dens_ratio
    fluid_properties.temp[1] = fluid_properties.temp[2] * dens_ratio

    n2, n3 = fluid_properties.dens[end-1], fluid_properties.dens[end-2]
    n1 = maximum([2 * n2 - n3, n2/10])
    dens_ratio = n2 / n1
    fluid_properties.dens[end] = n1
    fluid_properties.vel[end] = fluid_properties.vel[end-1] * dens_ratio
    fluid_properties.temp[end] = fluid_properties.temp[end-1] * dens_ratio
    
    return fluid_properties
end

function gather!(pc::ParticleContainer{T}, grid::Grid , E::Vector{T}) where T
	q_m = charge_to_mass(pc.species)
	@inbounds for i in eachindex(pc.pos)
        idx = abs(pc.inds[i])
        dist = (grid.face_centers[idx+1] - pc.pos[i])/ grid.dz
		pc.acc[i] = T(q_m * (E[idx+1] * (1 - dist) + E[idx] * dist))
	end
	return pc
end

function calc_electron_density_and_avg_charge!(n_e::Vector{T}, avg_charge::Vector{T}, species::Vector{SpeciesProperties{T}}) where T
    n_e .= zero(T)
    avg_charge .= zero(T)
    for sp in species
        Z = sp.species.charge
        if Z == 0
            continue
        end
        for (i, n_s) in enumerate(sp.dens)
            n_e[i] += n_s * Z
            avg_charge[i] += n_s
        end
    end
    @. avg_charge = n_e / avg_charge
    return n_e, avg_charge
end

function boltzmann_electric_field_and_potential!(E, phi, n_e, T_e, grid)

    # note: E is stored on edges
    # V - V_0 = Te ln(n) / ln(n_0)
    # -dPhi/dz = E = Te d(ln(n))/dz
    phi_0 = 0.0
    n_0 = n_e[2]

    for (i, n) in enumerate(n_e)
        phi[i] = phi_0 + T_e[i] * log(n / n_0)
    end

    for i in 2:length(E)-1
        # i -> left edge of cell i
        iL, iR = i-1, i
        z_L = grid.cell_centers[iL]
        z_R = grid.cell_centers[iR]
        E[i] = -(phi[iR] - phi[iL]) / (z_R - z_L)
    end

    return E, phi
end

#======================================================
Reactions definitions
======================================================#

struct Reaction{T<:AbstractFloat, L<:LinearInterpolation}
    rate_table::L
    products::Vector{Species}
    delta_n::Vector{T}
    delta_n_remainder::Vector{T}
    products_idx::Vector{UInt8}
    product_coefficients::Vector{UInt8}
    reactant::Species
    threshold_energy::T
    reactant_idx::UInt8
end



function deplete_reactant!(reaction::Reaction, reactant::ParticleContainer{T}, reactant_properties::SpeciesProperties{T}, products::Vector{ParticleContainer{T}}, product_properties::Vector{SpeciesProperties{T}}, grid::Grid, electron_properties::SpeciesProperties{T}, dt) where T

    # for each (non-ghost) cell, calculate the change in number density due to the reaction (delta n)
    for i in 2:length(grid.cell_centers)-1

        # find the rate from a lookup table 
        rate = reaction.rate_table(electron_properties.temp[i])
        # calculate the number of particles produced
        reaction.delta_n[i] = (reactant_properties.dens[i] * electron_properties.dens[i] * rate * dt) + reaction.delta_n_remainder[i]
    end

    # loop over the products to compute the weight consumed 
    # for each cell, add particles 
    fill!(reaction.delta_n_remainder, Inf)
    for (ip, product) in enumerate(product_properties)

        n_desired = products[ip].n_d # number of desired particles touching cell, hyperparameter from the simulation

        for i in 2:length(grid.cell_centers)-1
            # determine the number of partices to generate
            w_gen = product.avg_weight[i] * (product.N_particles[i]/n_desired) # check for particles in the cell 
            real_particles_produced = reaction.product_coefficients[ip] * reaction.delta_n[i] 
            n_gen = floor(Int32, real_particles_produced / (w_gen))

            reaction.delta_n_remainder[i] = minimum([reaction.delta_n_remainder[i], n_gen * w_gen])            
        end
    end

    for i in 2:length(grid.cell_centers)-1
        n_consumed = reaction.delta_n_remainder[i]
        # set the remainder and the delta n 
        reaction.delta_n_remainder[i] = reaction.delta_n[i] - n_consumed
        reaction.delta_n[i] = n_consumed
    end 


    # now that we have the delta_n, adjust particle weights 
    for i in 2:length(reactant.pos)-1

        # pull the index 
        s, ic = find_cell_indices(reactant.inds[i], grid)

        # for each cell the particle touches  
        for ii in 0:1

            # find the index 
            idx = ic + ii*s
            # calculate the shape function for portion of particle that touches the cell 
            shape = abs(reactant.pos[i] - grid.cell_centers[idx]) / grid.dz
            # relative proportion of weight in the cell 
            weight_proportion = reaction.delta_n[idx]/reactant_properties.dens[idx]
            # reduce the weight 
            reactant.weight[i] -= reactant.weight[i] * weight_proportion * shape
        end
        
    end

    return reaction, reactant 
end

"""
Note that the products vector should be sorted the same way as in the reaction structure 
"""
function generate_products!(products::Vector{ParticleContainer{T}}, product_properties::Vector{SpeciesProperties{T}},  reaction::Reaction{T}, reactant_properties::SpeciesProperties{T}, grid::Grid) where T

    pos_buf = T[]
    vel_buf = T[]
    weight_buf = T[]
    inds_buf = Int64[]

    dz = grid.face_centers[2:end] - grid.face_centers[1:end-1]
    # for each product 
    for (ip, product) in enumerate(products) 
        n_desired = product.n_d # number of desired particles touching cell, hyperparameter from the simulation
        # for each cell, add particles 
        for i in 2:length(grid.cell_centers)-1

            # determine the number of partices to generate
            w_gen = product_properties[ip].avg_weight[i] * (product_properties[ip].N_particles[i]/n_desired) # check for particles in the cell 
            real_particles_produced = reaction.product_coefficients[ip] * reaction.delta_n[i] 
            n_gen = Int32(floor(real_particles_produced / (w_gen)))

            # update the generation weight slightly to consume the full delta_n 
            w_gen = real_particles_produced / n_gen 

            # resize the buffers 
            resize!(pos_buf, n_gen)
            resize!(vel_buf, n_gen)
            resize!(weight_buf, n_gen)

            # generate the particles 
            # position 
            z_L = grid.face_centers[i]
            Random.rand!(pos_buf)
            pos_buf .= dz[i] * pos_buf .+ z_L
            inds_buf = copysign.(i, grid.cell_centers[i] .- pos_buf)

            # velocity, fix to use reactant mass
            v_th = sqrt(reactant_properties.temp[i] / product.species.gas.mass)
            Random.rand!(vel_buf)
            vel_buf .= vel_buf .* v_th .+ reactant_properties.vel[i]

            # weight
            weight_buf .= w_gen 
            # actual generation  
            add_particles!(product, pos_buf, vel_buf, weight_buf, inds_buf)
        end
    end 

    return products
end 

function read_reaction_rates(filepath)

    energy, values = open(filepath) do file 
        energy = parse(Float64, strip(split(readline(file), ':')[2]))
        values = readdlm(file, skipstart = 1)
        energy, values 
    end

    # include the n_0 normalization for how the rates are used 
    # to maintain the normalization, the two densities need to be multiplied by n_0 
    # and the entire addition divided by n_0, results in a net multiplication of n_0
    # timestep is assumed to be normalized by t_0, so multiply by t_0 as well 
    interp_object = LinearInterpolation(values[:,2] * n_0 * t_0 * 0.0, values[:,1])
    return energy, interp_object
end


function update_particles!(particles::Vector{ParticleContainer{T}}, reactions::Vector{R}, bulk_properties::Vector{SpeciesProperties{T}}, electrons::SpeciesProperties{T},  grid::Grid, dt::T) where {T, R<:Reaction}

    # reaction step 
    for reaction in reactions

        # find the relevant containers and properties 
        reactant = particles[reaction.reactant_idx]
        reactant_properties = bulk_properties[reaction.reactant_idx]
        products = particles[reaction.products_idx]
        product_properties = bulk_properties[reaction.products_idx]

        # first compute the density change and deplete the reactant  
        deplete_reactant!(reaction, reactant, reactant_properties, products, product_properties, grid, electrons ,dt)

        # now add additional particles 
        generate_products!(products, product_properties, reaction, reactant_properties, grid)

        # ensure that the overall containers are updated
        particles[reaction.reactant_idx] = reactant 
        particles[reaction.products_idx] = products

    end

    # locate and push 
    for container in particles
        push!(container, dt, grid)
        locate_particles!(container, grid)
    end 

    return particles, reactions
end

function update_electron_density!(electrons::SpeciesProperties{T}, heavy_species_properties::Vector{SpeciesProperties{T}}) where T
    
    # reset the electron density 
    electrons.dens .=0

    for (is, species) in enumerate(heavy_species_properties)
        # pull the density
        density = species.dens
        Z = species.species.charge 

        for i in 1:length(density)
            electrons.dens[i] += Z * density[i]
        end
    end

    return electrons 
end


function iterate!(particles::Vector{ParticleContainer{T}}, reactions::Vector{R}, bulk_properties::Vector{SpeciesProperties{T}}, electrons::SpeciesProperties{T},E_array::Vector{T}, Phi::Vector{T}, grid::Grid, dt::T) where {T, R<:Reaction}

    # first update the particles
    update_particles!(particles, reactions, bulk_properties, electrons, grid, dt)

    # deposit 
    for (ip, particle) in enumerate(particles)
        deposit!(bulk_properties[ip], particle, grid)
    end

    # update the electrons 
    # for now this is just the boltzmann relation 
    update_electron_density!(electrons, bulk_properties)

    boltzmann_electric_field_and_potential!(E_array, Phi, electrons.dens, electrons.temp, grid)

    # gather step 
    for (ip, particle) in enumerate(particles)
        gather!(particle, grid, E_array)
    end

    return particles, reactions, bulk_properties, electrons, E_array, Phi 
end


#======================================================
Output definitions
======================================================#
struct Output{T<:AbstractFloat}
    densities::Array{T}
    velocities::Array{T}
    temperatures::Array{T}
    phi::Matrix{T}
    E::Matrix{T}
    call_count::Vector{Int64}

end

function Output{T}(grid::Grid, n_species::Int64 , save_iterations::Int64) where T

    n_cells = length(grid.cell_centers)
    n_edges = length(grid.face_centers)

    cell_buf = Matrix{T}(undef,save_iterations, n_cells)
    edge_buf = Matrix{T}(undef,save_iterations, n_edges)

    species_buf = zeros(T, save_iterations, n_cells, n_species )
    

    return Output(species_buf, copy(species_buf), copy(species_buf), cell_buf, edge_buf, [1])

end

function save_output!(output::Output{T}, fluid_properties::Vector{SpeciesProperties{T}}, phi::Vector{T}, E::Vector{T}) where T
    
    row = output.call_count[1]

    # set the field and potential
    output.phi[row, :] .= phi
    output.E[row, :] .= E

    # save properties for all fluids 
    for (f, properties) in enumerate(fluid_properties)
        output.densities[row,:, f] .= properties.dens
        output.velocities[row,:, f] .= properties.vel
        output.temperatures[row,:, f] .= properties.temp
    end

    # iterate the call counter 
    output.call_count[1] = row + 1 

    return output 
end


end
