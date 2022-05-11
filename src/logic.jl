module 2DNFT

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

end # module