# Define the methods used to solve the activity partial-differential integral equations and the synaptic evolution equation. 

function nft_step(a, ui, W, S, tau_short, frs, frw)
	"""Solve the (i+1)th time-step using a Range-Kutta method. du/dt = -u + f_W(W * u) + f_S(S * a)"""
	input = act.(S * a; frs...)
	
	k1 = tau_short .* (-ui .+ act.(W * ui; frw...))
	k2 = tau_short .* (-(ui .+ 1/2 .* k1) .+ act.(W * (ui .+ 1/2 .* k1); frw...) .+ input)
	k3 = tau_short .* (-(ui .+ 1/2 .* k2) .+ act.(W * (ui .+ 1/2 .* k2); frw...) .+ input)
	k4 = tau_short .* (-(ui .+ 1/2 .* k3) .+ act.(W * (ui .+ 1/2 .* k3); frw...) .+ input)
	
	return ui .+ 1/6 .* k1 .+ 1/3 .* k2 .+ 1/3 .* k3 .+ 1/6 .* k4
end

function dSact(a, ut, Tact, H)
	"""Compute the convolutions for each retinal-colliculus location pair."""
	
	ds = 0
	for t = 1:Tact
		ind = Tact+2-t
		ds += a[t] * H[ind:ind+Tact]'*ut
	end		
	return 1/Tact * dS
end

function dS_step(Si, dSua, tau, kappa, lambda)
	"""Solve the (i+1)th synaptic time-step using Range-Kutta4 method. dS/dT = -lamda S + dSua + noise."""
	noise = kappa * randn(size(Si))
	
	k1 = tau .* (-lambda .* S .+ noise .+ dSua) 
	k2 = tau .* (-lambda .* (S .+ 1/2 .* k1) + noise .+ dSua)	
	k3 = tau .* (-lambda .* (S .+ 1/2 .* k2) + noise .+ dSua)	
	k4 = tau .* (-lambda .* (S .+ 1/2 .* k3) + noise .+ dSua)	
	return Si + 1/6 .* k1 .+ 1/3 .* k2 .+ 1/3 .* k3 .+ 1/6 .* k4
end

function solveS(retinal_inputs, Tact, Tsyn, W, S, Smask, tau_short, tau_long, lambda, kappa, H)
	"""Solve and record synaptic weight changes to S under retinal inputs with assumed recurrent connectivity W."""
	Srecord = zeros(size(S)[1], size(S)[2], Tsyn)
	Srecord[:,:,1] .= Array(S)
	u = CuArray(zeros(size(S)[1], Tact))
	ri = CuArray(retinal_inputs);
	for s = 1:Tsyn
		# solve the activity in this step
		ufinal = u[:,Tact]
		u = CuArray(zeros(size(S)[1], Tact))
		u[:,1] .= ufinal
		for t = 2:Tact
			u[:,t] .= nft_step(ri[:,t-1,s], u[:,t-1], W, S, tau_short)
		end
		# calculate synaptic changes
		dSua = dSact(a, u, Tact, H)
		S .= dS_step(S, dSua, tau_long, kappa, lambda) .* Smask
		Srecord[:,:,s] .= S
	end
end
