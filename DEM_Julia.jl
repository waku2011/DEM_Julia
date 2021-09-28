#
# DEM simulatior for Julia lang. practice
#
using Plots

const global num_particles = 10
const global paricle_type = 2 # particle and wall material

const global ms = 1.0  # mass of particle
const global rs = 0.01 # radius of particle
const global Is = fill(2.0/5.0*ms*rs*rs, 3)  # moment of Inertia of particle
const global g = [0, 0, -9.81] # gravity

# particle state
xs = rand(Float64, num_particles, 3) .- 0.5 # position vector 
vs = zeros(Float64, num_particles, 3)       # velocity vector
ws = zeros(Float64, num_particles, 3)       # rotational vector

# forces and torques


# Plots
plt3d = plot(xs[:,1],xs[:,2], xs[:,3], seriestype=:scatter, markersize = 5, tilte="Particle", label="")
display(plt3d)
readline()
