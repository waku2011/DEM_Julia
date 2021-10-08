#
# DEM simulatior for Julia lang. practice
#

using Printf
using Plots
ENV["GKSwstype"] = "nul"

const global num_active_particles = 20
const global num_boundary_particles = 100*100
const global num_particles = num_active_particles + num_boundary_particles
const global particle_type = 2 # active particle and boundary material

const global ms = 1.0  # mass of particle
const global rs = 0.01 # radius of particle
const global Is = fill(2.0/5.0*ms*rs*rs, num_particles, 3)  # moment of Inertia of particle
const global g = [0, 0, -9.81] # gravity
const global pID = fill(1, num_particles)

vb = [-1,+1];

# initialize particle state 

xs = rand(Float64, num_particles, 3) .- 0.5 # position vector xyz:[-0.5,0.5]
n = num_active_particles
for i in 1:100, j in 1:100
   global n +=1
   xs[n,1] = -0.5 + (1.0/100) * (i-0.5)  
   xs[n,2] = -0.5 + (1.0/100) * (j-0.5)  
   xs[n,3] = -0.5
end
vs = zeros(Float64, num_particles, 3)       # velocity vector
ws = zeros(Float64, num_particles, 3)       # rotational velocity vector
ts = zeros(Float64, num_particles, 3)       # rotation vector

force  = zeros(Float64, num_particles, 3)
torque = zeros(Float64, num_particles, 3)

xsnew = xs
vsnew = vs
wsnew = ws
tsnew = ts

# time stepping (symplectic Euler scheme)
deltaT = 0.01 # s
endTime = 1.0  # s
timeSteps = Integer(endTime/deltaT)

# ------------------------

function forces_update(force, state="init")
   if state == "init"
     # @printf("init \n")
     force  = zeros(Float64, num_particles, 3)
   else
     # @printf("update \n")
     force[1:num_active_particles,1] .= ms*g[1]
     force[1:num_active_particles,2] .= ms*g[2]
     force[1:num_active_particles,3] .= ms*g[3]
   end
end

function torque_update(torque, state="init")
   if state == "init"
     #@printf("init \n")
     torque  = zeros(Float64, num_particles, 3)
   else
     #@printf("update \n")
     torque[1:num_particles,1] .= 0.001
     torque[1:num_particles,2] .= 0.001
     torque[1:num_particles,3] .= 0.001
   end
end

# ------------------------

anim = Animation()

forces_update(force, "init")
torque_update(torque, "init")

for i in 1:timeSteps
  
  # forces and torque  update
  forces_update(force,"runtime")
  torque_update(torque,"runtime")
  
  # velocities update
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
