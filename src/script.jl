using CUDA, Plots

# Define the parameters

# Synaptic parameters

kappa = 0.1
lambda = 0.1 # syn mm^2
S0 = 1 # syn mm^-2
s1 = 0.3 # mm
s2 = 0.4 # mm

W0 = 1 syn mm ^ -2
W1 = 0.1 # syn mm^-2
r1 = 0.1 # mm

# Length and Time Scales 
tau_long = 1 # s
tau_short = 0.01 # s chosen to average bins of 10 in the PRANAS output of spikes to get a NFT smoothed profile

Tact = round(Int, 3/tau_short) # 3 seconds to cover a few length scales
Tsyn = round(Int, 3 * 85000 / Tact) # 3 days of activity 
L = 1000; # length scale is 1mm => a cell is 0.001mm


# Plasticity 

tp = 0.56 # s
H0 = 1 # syn

# Initialise
len_gen = 0:1/L:(L-1)
tact_gen = -1:1/Tact:1
retinal_mask = map((x,y) -> 1.0 * (x^2 + y^2 < 1), len_gen, len_gen);
collicular_mask = retinal_mask # simplicity
#THESE ARE WRONG
W = map((x,y) -> W0 * exp(-(x^2+y^2)/(2 * r0 ^ 2) + W1 * exp(-(x^2+y^2)/(2 * r1^2)), len_gen, len_gen)
S = map((x,y) -> S0 * exp(-(x^2+y^2)/(2 * s0 ^ 2) + S1 * exp(-(x^2+y^2)/(2 * s1^2)), len_gen, len_gen)
H = map(t -> sign(t) * exp(-abs(t)/tp), tact_gen)

# Logic

		
