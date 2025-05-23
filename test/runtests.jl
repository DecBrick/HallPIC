using HallPIC: HallPIC as hp
using Documenter
using Test
using DataInterpolations: LinearInterpolation

Documenter.doctest(hp)

const Hydrogen = hp.Gas(name=:H, mass=1)
const Xenon = hp.Gas(name=:Xe, mass=131.293)

@testset "Grid construction" begin
	N = 100
	left_boundary = hp.OpenBoundary()
	right_boundary = hp.WallBoundary(0.5)
	area = 1.0
	x0 = 0
	x1 = 1
	grid = hp.Grid(N, x0, x1, area, left_boundary, right_boundary)

	@test length(grid.cell_centers) == N+2
	@test length(grid.face_centers) == N+3
	@test grid.left_boundary == left_boundary
	@test grid.right_boundary == right_boundary
	dz = (x1 - x0) / N
	@test grid.dz ≈ dz  
	@test all(volume ≈ dz * area for volume in grid.cell_volumes)
	@test grid.face_centers[end-1] - grid.face_centers[2] == (x1 - x0)
end

function test_add_particles(::Type{T}) where T 
	# try initalizing the particles object 
	particles = hp.ParticleContainer{T}(0, Xenon(1))

	# check the initialization 
	@test hp.charge(particles.species) == 1
	@test hp.mass(particles.species) == Xenon.mass
	@test isempty(particles.pos)
	@test isempty(particles.vel)
	@test isempty(particles.acc)
	@test isempty(particles.weight)

	# try adding two particles 
	pos = [T(1.57), T(0.51)]
	vel = [T(19823.1041), T(981.471)]
	weight = [T(1.589e16), T(1.589e16)]
	inds = Int64.([1,2])

	particles = hp.add_particles!(particles, pos, vel, weight, inds)

	# check them 
	@test length(particles.pos) == 2
	@test particles.pos[1] == pos[1]
	@test particles.pos[2] == pos[2]

	@test length(particles.vel) == 2
	@test particles.vel[1] == vel[1]
	@test particles.vel[2] == vel[2]

	@test length(particles.weight) == 2
	@test particles.weight[1] == weight[1]
	@test particles.weight[2] == weight[2]

	@test length(particles.acc) == 2
	@test particles.acc[1] == 0
	@test particles.acc[2] == 0

	@test length(particles.inds) == 2
	@test particles.inds[1] == inds[1]
	@test particles.inds[2] == inds[2]
end

@testset "Add particles" begin
	for T in [Float32, Float64]
		test_add_particles(T)
	end
end

@testset "Boundary conditions" begin 

	# make a dummy particles object 
	particles = hp.ParticleContainer{Float64}(0, Xenon(1))
	pos = [1.57, 0.51, -1.00]
	vel = [19823.1041, 981.471, 4319.42]
	weight = [1.589e16, 1.589e16, 1.589e16]
	inds = Int64.([0, 0, 0])
	particles = hp.add_particles!(particles, pos, vel, weight, inds)
	og_particles = hp.add_particles!(particles, pos, vel, weight, inds)
	# Grid object
	N=10
	left_boundary = hp.OpenBoundary()
	right_boundary = hp.WallBoundary(0.5)
	area = 1.0
	x0 = 0
	x1 = 1
	grid = hp.Grid(N, x0, x1, area, left_boundary, right_boundary)


	# first test, initialize and call a wall boundary 
	left_boundary = hp.OpenBoundary()
	hp.apply_boundary!(particles, grid, left_boundary, Int8(-1))

	# check that the first two haven't been marked and the third has 
	@test particles.weight[1] ≈ 1.589e16
	@test particles.weight[2] ≈ 1.589e16
	@test particles.weight[3] ≈ 0.0

	# now check the right boundary 
	particles.weight[3] = 1.589e16
	right_boundary = hp.OpenBoundary()
	hp.apply_boundary!(particles, grid, right_boundary, Int8(0))

	# check that the last two haven't been marked and the first has 
	@test particles.weight[1] ≈ 0.0
	@test particles.weight[2] ≈ 1.589e16
	@test particles.weight[3] ≈ 1.589e16


	# now try a wall boundary 
	particles.weight[1] = 1.589e16
	boundary =  hp.WallBoundary(0.5)
	hp.apply_boundary!(particles, grid, boundary, Int8(0))

	# check that it's unaffected 
	@test particles.weight ≈ og_particles.weight

end

function test_leapfrog(::Type{T}) where T
	function leapfrog_gather!(pc::hp.ParticleContainer)
		@inbounds for i in eachindex(pc.pos)
			pc.acc[i] = -pc.pos[i]
		end
	end

	dt = 0.2
	x0 = 1.0
	v0 = 0
	tmax = 5*pi
	num_steps = ceil(Int, tmax / dt)
	t = range(0.0, tmax, length=num_steps)

	pos = zeros(num_steps)
	vel = zeros(num_steps)

	pc = hp.ParticleContainer{T}(1, Hydrogen(1))

	# Test 1: Initially-stationary particle remains at rest
	pc.pos[1] = 0.0
	pc.vel[1] = 0.0
	pc.weight[1] = 1e16

	# No acceleration at x = 0
	leapfrog_gather!(pc)
	@test pc.acc[1] == 0

	hp.push_vel!(pc, -dt/2)

	# No movement
	grid = hp.Grid(1, 0, 1, 1, hp.WallBoundary(0.5), hp.WallBoundary(0.5))
	hp.push!(pc, dt, grid)
	@test pc.pos[1] == 0
	@test pc.vel[1] == 0

	# Test 2: simple harmonic oscillator starting at (x,v) = (1,0)
	# Check that leapfrog is centered as expected
	pc.pos[1] = x0
	pc.vel[1] = v0
	leapfrog_gather!(pc)

	# half step backward in velocity
	hp.push_vel!(pc, -dt/2)
	pos[1] = pc.pos[1]
	vel[1] = pc.vel[1]

	for i in 2:num_steps
		leapfrog_gather!(pc)
		hp.push!(pc, dt, grid)

		pos[i] = pc.pos[1]
		vel[i] = pc.vel[1]
	end

	# Check that leapfrog is centered
	min, max = extrema(pos)
	@test min >= -1
	@test max <= 1

	min, max = extrema(vel)
	@test min >= -1
	@test max <= 1
end

@testset "Harmonic oscillator" begin
	for T in [Float32, Float64]
		test_leapfrog(T)
	end
end

function test_linear_drop(::Type{T}, V) where T
	pc = hp.ParticleContainer{T}(1, Hydrogen(1))
	M = 100
	L = 0.05
	x = LinRange(0, L, M)
	E = V / L * ones(T, M)

	x_norm = x ./ hp.x_0
	E_norm = T.(E ./ hp.E_0)

	#extrap = hp.DataInterpolations.ExtrapolationType.Constant
	#E_itp = hp.LinearInterpolation(E_norm, x_norm, extrapolation = extrap)

	q = hp.charge(pc.species) * hp.q_e
	m = hp.mass(pc.species) * hp.m_0

	u_expected = sqrt(2 * q * V / m)
	dt = 0.5 * L / u_expected / M
	dt_norm = dt / hp.t_0

	# pushback step
	grid = hp.Grid(1, 0, 1, 1, hp.WallBoundary(0.5), hp.WallBoundary(0.5))
	pc.weight[1] = 1e16
	pc.inds[1] = 1
	hp.gather!(pc, grid, E_norm)
	hp.push_vel!(pc, -dt_norm/2)

	t = 0.0
	while pc.pos[1] < x_norm[end] 
		t += dt_norm
		hp.gather!(pc, grid, E_norm)
		hp.push!(pc, dt_norm, grid)
	end

	hp.push_vel!(pc, dt_norm/2)

	a = hp.charge(pc.species) * E_norm[1] / hp.mass(pc.species)
	u_exact = T(a * t)
	x_exact = T(0.5 * a * t^2)

	@test x_exact ≈ pc.pos[1]
	@test u_exact ≈ pc.vel[1]
end

@testset "Linear potential drop" begin
	for T in [Float32, Float64]
		test_linear_drop(T, 600.0)
	end
end


@testset "SpeciesProperties" begin
	for T in [Float32, Float64]
		N = 5
		sp = hp.SpeciesProperties{T}(N, Hydrogen(1))
		@test length(sp) == N

		base = collect(1.0:N)
		@. sp.dens = 1.0 * base
		@. sp.vel = 2.0 * base
		@. sp.temp = 3.0 * base
		@. sp.avg_weight = 4.0*base

		@test eltype(sp) == hp.CellProperties{T}

		for (i, cell) in enumerate(sp)
			@test cell.dens == sp.dens[i] && cell.dens == sp[i].dens
			@test cell.vel == sp.vel[i] && cell.vel == sp[i].vel
			@test cell.temp == sp.temp[i] && cell.temp == sp[i].temp
			@test cell.avg_weight == sp.avg_weight[i] == sp[i].avg_weight
			@test cell.N_particles == 0
			@test cell == sp[i]
		end
		c = collect(sp)
		@test eltype(c) == hp.CellProperties{T}
		@test length(c) == N

		sp_empty = hp.SpeciesProperties{T}(0, Hydrogen(1))
		for _ in sp_empty
		end
		c = collect(sp_empty)
		@test eltype(c) == hp.CellProperties{T}
		@test length(c) == 0
		@test length(sp_empty) == 0
	end
end


@testset "Particle location" begin
	N = 100
	x0 = 0.
	x1 = 1.0*N
	grid = hp.Grid(N, x0, x1, 1.0)
	dx = (x1 - x0) / N
	pc = hp.ParticleContainer{Float32}(7, Hydrogen(1))

	pc.pos[1] = 50.5
	pc.pos[2] = 50.1
	pc.pos[3] = 50.7
	pc.pos[4] = x0
	pc.pos[5] = -0.6
	pc.pos[6] = dx
	pc.pos[7] = x1

	hp.locate_particles!(pc, grid)

	@test pc.inds[1] == 52
	@test pc.inds[2] == -52
	@test pc.inds[3] == 52
	@test pc.inds[4] == -2
	@test pc.inds[5] == -1
	@test pc.inds[6] == -3
	@test pc.inds[7] == -1*(N+2)
end

function test_deposition(::Type{T}, n_cell, profile = :uniform, rtol = 0.01) where T
	# initialize species 
	xenon = hp.Gas(name=:Xe, mass=131.293)
	species = xenon(0)

	# define the cell properties
	L = 1.0
	grid = hp.Grid(n_cell, 0, L, 1.0)
	vol = grid.cell_volumes[1]
	particles_per_cell = 500_000 # per cell 

	# Initialize fluid property arrays
	n_0 = 1e18 / hp.n_0
	u_0 = 3.0 
	T_0 = 50.0
	w_0 = n_0 / particles_per_cell * vol 
	avg_interval = 10

	uniform(x) = 1
	linear(x) = 0.5 * (1+x/L)

	fluid_properties = hp.SpeciesProperties{Float64}(n_cell+2, species)
	prof_z = if profile == :linear
		linear.(grid.cell_centers)
	else
		uniform.(grid.cell_centers)
	end
	@. fluid_properties.dens = prof_z * n_0
	@. fluid_properties.vel = prof_z * u_0 
	@. fluid_properties.temp = prof_z * T_0 
	fluid_properties.avg_weight .= 0.0
	fluid_properties.N_particles .= 0

	n_exact = copy(fluid_properties.dens)
	u_exact = copy(fluid_properties.vel)
	T_exact = copy(fluid_properties.temp)
	w_exact = n_exact / particles_per_cell * vol

	# Create particles from fluid properties
	particles = hp.initialize_particles(fluid_properties, grid, particles_per_cell)

	# check that particles are in bounds and limits are correct
	@test size(particles.pos)[1] == n_cell * particles_per_cell
	min, max = extrema(particles.pos)
	@test max <= L
	@test min >= 0

	# verify that particle positions are sorted
	@test issorted(particles.pos)

	# verify that particles are evenly-spaced
	dz = diff(particles.pos)
	@test all(_dz ≈ dz[1] for _dz in dz)

	min, max = extrema(particles.weight)
	w_min = 0.5 * (w_exact[1] + w_exact[2])
	@test isapprox(max, w_0; rtol)
	if profile == :uniform
		@test isapprox(min, w_0; rtol)
	else
		@test isapprox(min, w_min; rtol)
	end

	min, max = extrema(particles.vel)
	thermal_speed = sqrt(T_0 / Xenon.mass)
	@test max <= u_0 + 8 * thermal_speed
	@test min >= u_0 - 8 * thermal_speed

	# deposit to the cells 
	hp.locate_particles!(particles, grid)
	hp.deposit!(fluid_properties, particles, grid, avg_interval)

	# check that interior cells have correct properties and ghost cells are still zero
	inds = (firstindex(n_exact)+1) : (lastindex(n_exact)-1)
	@test all(
		isapprox.(fluid_properties.dens[inds], n_exact[inds]; rtol)
	)

	@test all(
		isapprox.(fluid_properties.avg_weight[inds], w_exact[inds]/(2*avg_interval); rtol)
	)

	@test all(
		isapprox.(fluid_properties.vel[inds], u_exact[inds]; rtol)
	)

	@test all(
		isapprox.(fluid_properties.temp[inds], T_exact[inds]; rtol)
	)

	# check ghost cell extrapolation
	if fluid_properties.dens[2] > fluid_properties.dens[1]
		@test fluid_properties.vel[2] < fluid_properties.vel[1]
	else
		@test fluid_properties.vel[2] >= fluid_properties.vel[1]
	end

	if fluid_properties.dens[end-1] > fluid_properties.dens[end]
		@test fluid_properties.vel[end-1] < fluid_properties.vel[end]
	else
		@test fluid_properties.vel[end-1] >= fluid_properties.vel[end]
	end
end

@testset "Deposition" begin
	num_cells = 10
	profiles = [:uniform, :linear]
	tolerances = [2e-2, 5e-2]
	for T in [Float32, Float64]
		for (prof, tol) in zip(profiles, tolerances)
			test_deposition(T, num_cells, prof, tol)
		end
	end
end

@testset "Electron density and average charge" begin
	species = [
		Xenon(0), Xenon(1), Xenon(2), Xenon(3), Xenon(-1)
	]

	n_0 = 1f18 / hp.n_0

	densities = (length(species):-1:1) .* n_0
	n_e_exact = sum(sp.charge * dens for (sp, dens) in zip(species, densities))
	n_i_exact = sum(dens for (sp, dens) in zip(species, densities) if sp.charge != 0)
	avg_Z_exact = n_e_exact / n_i_exact

	num_cells=10

	for T in [Float32, Float64]
		fluid_containers = [
			hp.SpeciesProperties{T}(num_cells, sp) for sp in species
		]

		for (fluid, dens) in zip(fluid_containers, densities)
			fluid.dens .= dens
		end

		n_e = zeros(T, num_cells)
		avg_Z = zeros(T, num_cells)
		hp.calc_electron_density_and_avg_charge!(n_e, avg_Z, fluid_containers)

		@test all(n_e .== T(n_e_exact))
		@test all(avg_Z .== T(avg_Z_exact))
	end
end

@testset "Boltzmann relation 1" begin
	num_cells = 10
	L = 1.0
	grid = hp.Grid(num_cells, 0, L, 1.0)

	T = Float32

	n_0 = 1e18 / hp.n_0
	n_e = ones(T, num_cells+2) * n_0
	T_e = ones(T, num_cells+2) * 5.0

	# check that uniform density causes phi == 0 and E == 0 everywhere
	phi = zeros(T, num_cells+2)
	E = zeros(T, num_cells+3)

	hp.boltzmann_electric_field_and_potential!(E, phi, n_e, T_e, grid)

	@test all(phi .== zero(T))
	@test all(E .== zero(T))

	# check that quadratic density causes large phi in middle and E pointing outward
	quadratic(x) = 1 - 2 * (x/L - 0.5)^2
	@. n_e = quadratic(grid.cell_centers) * n_0

	hp.boltzmann_electric_field_and_potential!(E, phi, n_e, T_e, grid)
	atol = 1e-12
	@test isapprox(phi[2], zero(T); atol)
	@test isapprox(phi[num_cells+1], zero(T); atol)
	n_ratio = maximum(n_e) / n_e[2]
	@test maximum(phi) ≈ T_e[1] * log(n_ratio)
	# note: this is true if we have an even number of cells,
	# since an edge will lie exactly in the middle
	midpt = length(E) ÷ 2 + 1
	@test isapprox(E[midpt], zero(T); atol)
	@test all(E[2:midpt-1] .< 0)
	@test all(E[midpt+1:end-1] .> 0)
end

@testset "Partition" begin
	N = 100
	num_trials = 10
	for _ in 1:num_trials
		vec = rand(Bool, N)
		num_zero = count(==(0), vec)
		num_one = length(vec) - num_zero
		new_size = hp.partition!(vec)
		@test new_size == num_zero
		@test all(!vec[i] for i in 1:new_size)
		@test all(vec[i] for i in new_size+1:length(vec))

		# test partitioning an already partitioned list
		vec2 = copy(vec)
		size_2 = hp.partition!(vec2)
		@test new_size == size_2
		@test all(v1 == v2 for (v1, v2) in zip(vec, vec2))
	end

	# test on all trues
	vec = ones(Bool, N)
	new_size = hp.partition!(vec)
	@test new_size == 0
	@test all(vec)

	# test on all falses
	vec = zeros(Bool, N)
	new_size = hp.partition!(vec)
	@test new_size == N
	@test all(!v for v in vec)
end

@testset "Removing particles" begin
	N = 100
	particles_per_cell = 50
	T = Float32
	grid = hp.Grid(N, 0, 1, 1)
	dens = fill(T(1e6), N+2)
	vel = ones(T, N+2)
	temp = ones(T, N+2)
	weights = ones(T, N+2)
	N_particles = 2 * particles_per_cell * ones(UInt64,N+2)
	species = Xenon(1)
	fc = hp.SpeciesProperties{T}(dens, vel, temp, weights, N_particles, species)
	
	pc = hp.initialize_particles(fc, grid, particles_per_cell)
	hp.locate_particles!(pc, grid)

	@test length(pc) == particles_per_cell * N
	@test firstindex(pc) == 1
	@test lastindex(pc) == length(pc)

	# try removing all particles in a specific cell
	cell_ind = 23
	pre_count = count(x->abs(x) == cell_ind, pc.inds)
	@test pre_count == particles_per_cell

	pre_flag_count = count(==(0), pc.weight)
	@test pre_flag_count == 0

	hp.flag_particles_in_cell!(pc, cell_ind)

	flag_count = count(==(0), pc.weight)
	@test flag_count == particles_per_cell

	hp.remove_flagged_particles!(pc)
	@test length(pc) == N*particles_per_cell - particles_per_cell
	@test length(pc) == length(pc.inds)
	@test length(pc) == length(pc.weight)
	@test length(pc) == length(pc.vel)
	@test length(pc) == length(pc.pos)
	@test length(pc) == length(pc.acc)
	post_count = count(x->abs(x) == cell_ind, pc.inds)
	@test post_count == 0
end

function test_read_reaction_table()
	# load the reaction 
	filepath = "../reactions/ionization_Xe_Xe+.dat"
	threshold_energy, table = hp.read_reaction_rates(filepath)

	# ensure that the values are correct 
	@test threshold_energy ≈ 12.1298437

	@test table.t[1] ≈ 0.3878E-01
	@test table.t[100] ≈ 499.5

	@test table.u[1] ≈ 0.000
	@test table.u[100] ≈ 0.3526E-12*hp.n_0*hp.t_0

end


function test_initialize_reaction(::Type{T}) where T

	# load the rate table 
	filepath = "../reactions/ionization_Xe_Xe+.dat"
	threshold_energy, table = hp.read_reaction_rates(filepath)

	# define the species  
	reactant = Xenon(0)
	product = [Xenon(1)]

	# initialize the struct 
	Xe_ionization = hp.Reaction{T, LinearInterpolation}(table, product, [0.0, 0.0, 0.0, 0.0],[0.0, 0.0, 0.0, 0.0],  UInt8.([0]),UInt8.([1]), reactant, threshold_energy, UInt8(0))

	# actually test 
	@test Xe_ionization.reactant == reactant
	@test Xe_ionization.products == product 
	@test Xe_ionization.product_coefficients == [1]
	@test Xe_ionization.threshold_energy ≈ threshold_energy
	@test Xe_ionization.rate_table == table
	@test Xe_ionization.delta_n ≈ zeros(4) 


end


function test_reaction_step(reactant_gas, product_gases, product_coefficients, reaction_path, TeV, ::Type{T}, rtol=0.01) where T
	
	# first set up the plasma 
	# seed properties/initialize cell arrays 
	n_cell = 10 
	n_n = 500
	grid = hp.Grid(n_cell, 0, 1.0, 1.0, hp.OpenBoundary(), hp.OpenBoundary())
	n_particles = n_n * n_cell 
	
	# initialize the reactant 
	reactant_properties = hp.SpeciesProperties{T}(n_cell+2, reactant_gas)
	reactant_properties.dens .= 1e18 / hp.n_0# 1/m^3
	reactant_properties.vel .= 300 ./ hp.u_0
	reactant_properties.temp .= (500 / 11604) # eV
	reactant_properties.avg_weight .= 0
	reactant_properties.N_particles .= 0
	reactant = hp.initialize_particles(reactant_properties, grid, n_n)
	# deposit to grid 
	hp.locate_particles!(reactant, grid)
	hp.deposit!(reactant_properties, reactant, grid)

	# copy properties for tests 
	old_reactant_density = copy(reactant_properties.dens)
	old_weights = copy(reactant.weight)

	product_properties = Vector{hp.SpeciesProperties{T}}(undef, length(product_gases))
	products = Vector{hp.ParticleContainer{T}}(undef, length(product_gases))

	for (ip, product) in enumerate(product_gases)

		# initialize 
		properties = hp.SpeciesProperties{T}(n_cell+2, product)
		properties.dens .=  1e16 / hp.n_0# 1/m^3
		properties.vel .= 5000 / hp.u_0
		properties.temp .=  0.1 # eV
		properties.avg_weight .= 0
		properties.N_particles .= 0
		

		product = hp.initialize_particles(properties, grid, n_n)
		

		# deposit to grid 
		hp.locate_particles!(product, grid)
		hp.deposit!(properties, product, grid) 


		# save to array 
		product_properties[ip] = properties
		products[ip] = product 
	end
	
	# now can initialization for reaction properties 
	# load the rate table 
	threshold_energy, table = hp.read_reaction_rates(reaction_path)
	reaction = hp.Reaction{T, LinearInterpolation}(table, product_gases,  zeros(n_cell+2), zeros(n_cell+2), zeros(UInt8, length(product_gases)), product_coefficients, reactant_gas, threshold_energy, 0)

	# initialize some electron properties
	electron = hp.Gas(name=:e, mass=hp.m_e)
	electron_properties = hp.SpeciesProperties{T}(n_cell+2, electron(-1))
	electron_properties.temp .= TeV # choose 10eV for now 
	electron_properties.dens .= product_properties[1].dens # quasineutrality 

	# reduce weights 
	dt = 15e-9 / hp.t_0
	reaction, reactant = hp.deplete_reactant!(reaction, reactant, reactant_properties, products, product_properties, grid, electron_properties, dt) 

	# check that number is conserved 
	hp.deposit!(reactant_properties, reactant, grid)
	rate = reaction.rate_table(TeV)
	delta_ns = zeros(T, n_cell+2)
	delta_n_vals = zeros(T, n_cell+2)
	for i in 2:n_cell+1

		delta_n = dt * electron_properties.dens[i] * old_reactant_density[i] * rate
		# find the minimum 
		n_consumed = Inf
        for (ip, product) in enumerate(product_properties)

            n_desired = products[ip].n_d # number of desired particles touching cell, hyperparameter from the simulation
        
            # determine the number of partices to generate
            w_gen = product.avg_weight[i] * (product.N_particles[i]/n_desired) # check for particles in the cell 
            real_particles_produced = reaction.product_coefficients[ip] * delta_n
            n_gen = Int32(floor(real_particles_produced / (w_gen)))

            n_consumed = minimum([n_consumed, n_gen * w_gen])            

        end
		delta_ns[i] = n_consumed
		delta_n_vals[i] = delta_n
	end
		
	@test all(reaction.delta_n .≈ delta_ns)
	@test all(isapprox(reactant_properties.dens, old_reactant_density .- delta_n_vals; rtol))

	# check that weights are reduced as expected 
	@test all(old_weights .>= reactant.weight)
	

	# add particles 
	products = hp.generate_products!(products, product_properties,reaction, reactant_properties, grid)

	
	# do the tests for each product 
	for (ip, product) in enumerate(products)

		# check that number is conserved 
		old_density = copy(product_properties[ip].dens)
		hp.locate_particles!(product, grid)
		hp.deposit!(product_properties[ip], product, grid) 
		
		@test all(isapprox(product_properties[ip].dens, old_density + delta_ns; rtol))

		# final check that the number of ions has expanded 
		@test length(product.pos) > n_particles 
	end
end

@testset "Reactions" begin
	test_read_reaction_table()

	"""
	First test: Xe ionization 
	Test to ensure that the xenon neutral to singly charged reaction works correctly which consists of: 
	1. Ensure that the correct reaction rate is determined
	2. Ensure that the neutral weight is removed according to this rate
	3. Ensure that the ion weight is added according to this weight (mass is conserved)
	4. Ensure that ion particles are added 

	This test case is to ensure the overall reaction functions are behaving as expected 
	"""
	tol = 1e-2

	for T in [Float32, Float64]
		test_initialize_reaction(T)
		test_reaction_step(Xenon(0), [Xenon(1)], [1], "../reactions/ionization_Xe_Xe+.dat", 10, T, tol)
	end

	"""
	Second test N2 dissociation
	The purpose of this test is to ensure that the reaction functions
	correctly handle reactions where multiple instances of the same
	product are produced 

	Tests follow the structure of the first test 

	Rates & threshold_energy from J. Chem. Phys. 98, 9544–9553 (1993)
	https://doi.org/10.1063/1.464385

	"""

	N2 = hp.Gas(name=:N2, mass=28.014)
	N = hp.Gas(name=:N, mass=14.007)

	for T in [Float32, Float64]
		test_reaction_step(N2(0), [N(0)], [2], "../reactions/Dissociation_N2.dat", 10, T, tol)
	end

	"""
	Third test NO dissociation
	The purpose of this test is to ensure that the reaction functions
	correctly handle reactions where multiple products are produced 

	Tests follow the structure of the first test 

	Rates from M Castillo et al 2004 Plasma Sources Sci. Technol. 13 343
	Threshold energy from  M Castillo et al 2004 Plasma Sources Sci. Technol. 13 39

	"""

	NO = hp.Gas(name=:NO, mass=30.006)
	O = hp.Gas(name=:O, mass=15.999)

	for T in [Float32, Float64]
		test_reaction_step(NO(0), [N(0), O(0)], [1,1], "../reactions/Dissociation_NO.dat", 60, T, tol)
	end
end

