using HallPIC: HallPIC as hp 
using Plots
using Statistics
using JSON3
using BenchmarkTools

#constants 
ϵ_0 = 8.85418782e-12 
μ_0 = 1.25663706e-6 

# Species 
ion = hp.Gas(name=:ion, mass=100*hp.m_e)
electron = hp.Gas(name=:e, mass=hp.m_e)

# general initialization 
n_cell = 256 

# the following is to compute the cell spacing, based on the Alven speed 
B0 = 0.1
v_A = 3e5
n_p = (B0 / v_A)^2 / (μ_0 * hp.m_0 * (ion.mass + hp.m_e)) 
w_pi = sqrt( hp.q_e^2 * n_p / (hp.m_0 * ion.mass * ϵ_0))
li = 3e8 / w_pi 
dz = (1.0 / 6.0) * l_i
L_z = dz * n_cell

grid = hp.Grid(n_cell, 0, L_z / hp.x_0, L_z / hp.x_0^2, hp.PeriodicBoundary(), hp.PeriodicBoundary())

E = zeros(T, length(grid.face_centers))
phi = zeros(T, length(grid.cell_centers))

#time initialization
w_ci = hp.q_0 * B0 / (ion.mass*hp.m_0)
LT = 40.0 
DT = 1e-3 

n_iterations = Int(LT / DT)
dt = DT / w_ci / hp.t_0

#damping properties 
epsilon = 0.03
m = 4 
L_z = 

    


# initialize ions 
n_ions = 8192 

ion_properties = hp.SpeciesProperties{T}(n_cell+2, Xenon(1))
ion_properties.dens .=  1e16 / hp.n_0# 1/m^3
ion_properties.vel .= 0.0
ion_properties.temp .=  1 # eV
ion_properties.avg_weight .= 0
ion_properties.N_particles .= 0

ions = hp.initialize_particles(ion_properties, grid, n_ions)

# deposit to grid 
hp.locate_particles!(ions, grid)
hp.deposit!(ion_properties, ions, grid) 


# initialize some electron properties

electron_properties = hp.SpeciesProperties{T}(n_cell+2, electron(-1))
electron_properties.temp .= 30 # choose 10eV for now 
electron_properties.dens .= ion_properties[1].dens # quasineutrality 

# initalize the reaction 
# load the rate table 

#threshold_energy, table = hp.read_reaction_rates("../reactions/ionization_Xe_Xe+.dat")
threshold_energy, table = hp.read_reaction_rates("C:/Users/brickd/.julia/dev/HallPIC/reactions/ionization_Xe_Xe+.dat")
reaction = hp.Reaction{T, typeof(table)}(table, [Xenon(1)], zeros(n_cell+2), zeros(n_cell+2), [2], [1], Xenon(0), threshold_energy, 1)


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