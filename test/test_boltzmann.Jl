using HallPIC: HallPIC as hp 
using Plots
using Statistics
using JSON3
using BenchmarkTools

const Xenon = hp.Gas(name=:Xe, mass=131.293)

function Boltzmann_basic_test()
	T = Float64
	# general initialization 
	n_iterations = 1000
	n_cell = 100
	n_n = 100
	dt = T(75e-9 / hp.t_0)
	grid = hp.Grid(n_cell, 0, 0.1 / hp.x_0, 0.1 / hp.x_0^2, hp.OpenBoundary(), hp.OpenBoundary())
	E = zeros(T, length(grid.face_centers))
	phi = zeros(T, length(grid.cell_centers))

	# initialize neutrals 
	neutral_properties = hp.SpeciesProperties{T}(n_cell+2, Xenon(0))
	neutral_properties.dens .= 1e19 / hp.n_0# 1/m^3
	neutral_properties.vel .= 0.0
	neutral_properties.temp .= (500 / 11604) # eV
	neutral_properties.avg_weight .= 0
	neutral_properties.N_particles .= 0
	neutrals = hp.initialize_particles(neutral_properties, grid, n_n)
	# deposit to grid 
	hp.locate_particles!(neutrals, grid)
	hp.deposit!(neutral_properties, neutrals, grid)


	# initialize ions 
	ion_properties = hp.SpeciesProperties{T}(n_cell+2, Xenon(1))
	ion_properties.dens .=  1e16 / hp.n_0# 1/m^3
	ion_properties.vel .= 0.0
	ion_properties.temp .=  0.1 # eV
	ion_properties.avg_weight .= 0
	ion_properties.N_particles .= 0

	ions = hp.initialize_particles(ion_properties, grid, n_n)

	# deposit to grid 
	hp.locate_particles!(ions, grid)
	hp.deposit!(ion_properties, ions, grid) 


	# initialize some electron properties
	electron = hp.Gas(name=:e, mass=hp.m_e)
	electron_properties = hp.SpeciesProperties{T}(n_cell+2, electron(-1))
	electron_properties.temp .= 30 # choose 10eV for now 
	electron_properties.dens .= ion_properties[1].dens # quasineutrality 

	# initalize the reaction 
	# load the rate table 

	#threshold_energy, table = hp.read_reaction_rates("../reactions/ionization_Xe_Xe+.dat")
	threshold_energy, table = hp.read_reaction_rates("C:/Users/brickd/.julia/dev/HallPIC/reactions/ionization_Xe_Xe+.dat")
	reaction = hp.Reaction{T}(table, [Xenon(1)], zeros(n_cell+2), zeros(n_cell+2), [2], [1], Xenon(0), threshold_energy, 1)


	# put everything into the format expected by the iterator 
	particles = [neutrals, ions]
	bulk_properties = [neutral_properties, ion_properties]
	reactions = [reaction]

	# initialize output object
	outputs = hp.Output{T}(grid, 2, n_iterations)

	# run loop 
	for i in 1:n_iterations 
		hp.iterate!(particles, reactions, bulk_properties, electron_properties, E, phi, grid, dt)

		hp.save_output!(outputs, bulk_properties, phi, E)
	end


	#display(plot(grid.cell_centers, Statistics.mean(outputs.phi, dims=1)[:]))

	#open("boltzmann_output.json", "w") do f 
#		JSON3.print(f, outputs)
#	end
end 

#Boltzmann_basic_test()
#@benchmark Boltzmann_basic_test() seconds = 20
@profview Boltzmann_basic_test()