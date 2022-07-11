using CUDA, Plots

# Define the parameters

# Synaptic parameters

kappa = 0.1
lambda = 0.1 # syn mm^2
S0 = 1 # syn mm^-2
s1 = 0.3 # mm
s2 = 0.4 # mm

W0 = 1 # syn mm ^ -2
W1 = 0.1 # syn mm^-2
r1 = 0.1 # mm

# Plasticity 

tp = 0.56 # s
H0 = 1 # syn

# Activation functions
retinal_activation = Dict(:rate => 100, :threshold=>0.5, :slope=>0.01)
colliculus_activation = Dict(:rate => 90, :threshold=>0.5, :slope=>0.01)

# Length and Time Scales 
tau_long = 1 # s
tau_short = 0.01 # s chosen to average bins of 10 in the PRANAS output of spikes to get a NFT smoothed profile

Tact = round(Int, 3/tau_short) # 3 seconds to cover a few length scales
Tsyn = round(Int, 3 * 85000 / Tact) # 3 days of activity 
L = 1000; # length scale is 1mm => a cell is 0.001mm

