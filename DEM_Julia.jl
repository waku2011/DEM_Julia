#
# DEM simulatior for Julia lang. practice
#

using Printf
using Plots
ENV["GKSwstype"] = "nul"

const global num_particles = 10
const global paricle_type = 2 # particle and wall material

const global ms = 1.0  # mass of particle
const global rs = 0.01 # radius of particle
const global Is = fill(2.0/5.0*ms*rs*rs, num_particles, 3)  # moment of Inertia of particle
const global g = [0, 0, -9.81] # gravity

vb = [-1,+1];

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
deltaT = 0.01 # s
endTime = 1.0  # s
timeSteps = Integer(endTime/deltaT)

anim = Animation()

for i in 1:timeSteps
  
  # velocity update
  global vsnew = vs + force/ms * deltaT
  
  
  for n=1:num_particles
    global wsnew[n,1] = ws[n,1] + torque[n,1]./Is[n,1]*deltaT
    global wsnew[n,2] = ws[n,2] + torque[n,2]./Is[n,2]*deltaT
    global wsnew[n,3] = ws[n,3] + torque[n,3]./Is[n,3]*deltaT
  end
     
  # position and rotation angle update
  global xsnew = xs + vsnew * deltaT
  global tsnew = ts + wsnew * deltaT

  #@printf("%d %f  \n", i, minimum(xsnew))
  
  # time increment
  global xs = xsnew
  global vs = vsnew
  global ws = wsnew
  global ts = tsnew 
   
  # Plots
  plot(xs[:,1],xs[:,2], xs[:,3], xlim=vb, ylim=vb, zlim=vb, seriestype=:scatter, markersize=6, title="Particle", legend=:none)
  #display(plt3d)
  frame(anim)
  
end
gif(anim)
