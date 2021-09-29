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
xs = rand(Float64, num_particles, 3) .- 0.5 # position vector xyz:[-0.5,0.5]
vs = zeros(Float64, num_particles, 3)       # velocity vector
ws = zeros(Float64, num_particles, 3)       # rotational velocity vector
ts = zeros(Float64, num_particles, 3)       # rotation vector

xsnew = xs
vsnew = vs
wsnew = ws
tsnew = ts

# forces and torques

force  = zeros(Float64, num_particles, 3)
force[1:num_particles,1] .= ms*g[1]
force[1:num_particles,2] .= ms*g[2]
force[1:num_particles,3] .= ms*g[3]
torque = zeros(Float64, num_particles, 3)

# time stepping (symplectic Euler scheme)
deltaT = 0.0001 # s
endTime = 1.0   # s
timeSteps = Integer(endTime/deltaT)

a = Animation()

for i in 1:timeSteps
  
  # velocity update
  global vsnew = vs + force/ms * deltaT
  global wsnew[:,1] .= ws[:,1] .+ torque[:,1]./Is[:,1]*deltaT
  global wsnew[:,2] .= ws[:,2] .+ torque[:,2]./Is[:,2]*deltaT
  global wsnew[:,3] .= ws[:,3] .+ torque[:,3]./Is[:,3]*deltaT
  
  # position and rotation angle update
  global xsnew = xs .+ vsnew * deltaT
  global tsnew = ts .+ wsnew * deltaT

  # time increment
  global xs = xsnew
  global vs = vsnew
  global ws = wsnew
  global ts = tsnew 
   
  # Plots
  plt3d = plot(xs[:,1],xs[:,2], xs[:,3], seriestype=:scatter, markersize=5, title="Particle", label="")
  frame(a, plt3d)
  
end

gif(a)
