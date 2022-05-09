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

function nft_step(a, ui, W, S, tau_short)
	"""Solve the (i+1)th time-step using a Range-Kutta method. du/dt = -u + W * u + S * a"""
	input = S * a
	
	k1 = tau_short .* (-ui .+ W * ui)
	k2 = tau_short .* (-(ui .+ 1/2 .* k1) .+ W * (ui .+ 1/2 .* k1) .+ input)
	k3 = tau_short .* (-(ui .+ 1/2 .* k2) .+ W * (ui .+ 1/2 .* k2) .+ input)
	k4 = tau_short .* (-(ui .+ 1/2 .* k3) .+ W * (ui .+ 1/2 .* k3) .+ input)
	
	return ui .+ 1/6 .* k1 .+ 1/3 .* k2 .+ 1/3 .* k3 .+ 1/6 .* k4
end

function dSact(at, ut, Tact, H)
	"""Compute the convolutions for each retinal-colliculus location pair."""
	
	ds = 0
	for t = 1:Tact
		ind = Tact+2-t
		ds += a[t] * @view H[ind:ind+Tact]'*ut
	end		
	return 1/Tact * dS
end

function dS_step(at, ut, Si, dSua, tau, kappa, lambda)
	"""Solve the (i+1)th synaptic time-step using Range-Kutta4 method. dS/dT = -lamda S + dSua + noise."""
	noise = kappa * randn(size(Si))
	
	k1 = tau .* (-lambda .* S .+ noise .+ dSua) 
	k2 = tau .* (-lambda .* (S .+ 1/2 .* k1) + noise .+ dSua)	
	k3 = tau .* (-lambda .* (S .+ 1/2 .* k2) + noise .+ dSua)	
	k4 = tau .* (-lambda .* (S .+ 1/2 .* k3) + noise .+ dSua)	
	return Si + 1/6 .* k1 .+ 1/3 .* k2 .+ 1/3 .* k3 .+ 1/6 .* k4
end

function solveS(retinal_inputs, Tact, Tsyn, W, S, tau_short, tau_long, lambda, kappa, H)
	"""Solve and record synaptic weight changes to S under retinal inputs with assumed recurrent connectivity W."""
	Srecord = zeros(size(S)[1], size(S)[2], Tsyn)
	Srecord[:,:,1] .= Array(S)
	u = CuArray(zeros(size(S)[1], Tact))
	ri = CuArray(retinal_inputs);
	for s = 1:Tsyn
		# solve the activity in this step
		a = @view ri[:,s]
		ufinal = u[:,Tact]
		u = CuArray(zeros(size(S)[1], Tact))
		u[:,1] .= ufinal
		for t = 2:Tact
			u[:,t] .= nft_step(a[:,t-1], u[:,t-1], W, S, tau_short)
		end
		
		# calculate synaptic changes
		dSua = dSact(a, u, Tact, H)
		S .= dS_step(a, u, dSau, tau, kappa, lambda)
		Srecord[:,:,s] .= S
	end
end
		
