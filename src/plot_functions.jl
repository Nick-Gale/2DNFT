function receptive_fields(retinal_hash, colliculus_hash, params, img_dir)
	"""Plot the retinal location which a given collicular location is most receptive under stimulation."""
	Tact = params["Tact"]
	base_a = params["retinal_activation"]["threshold"] .+ 0.1*slope
	frs = params["retinal_activation"]
	frw = params["colliculus_activation"]

	records = zeros(length(collicular_hash), length(retinal_hash))
	for (i, loc) in retinal_hash
		a = zeros(length(retinal_hash), Tact)
		u = zeros(length(colliculus_hash), Tact)
		a[i,:] .= base_a
		for t = 1:(Tact-1)
			u[:,t+1] = nft_step(a[:,t], u[:,t], W, S, tau_short, frw, frs)
		end
		for (j, loc) in colliculus_hash
			records[j, i] = u[j, end]
		end
	end

	maximal_ret_inds = map(x -> argmax(records[x, :]), 1:size(records[1]))

	color_map = map(x -> RGB(retinal_hash[x][1], 0.45, retinal_hash[x][2]), maximal_ret_inds)

	plt = plot(st=:scatter, DPI=400)
	for (j, loc) in collicular_hash
		plot!(plt, loc[1], loc[2], c=color_map[j])
	end
	
	savefig(img_dir, plt)
	return nothing
end

function plot_direction_selectivity()
	"""Plot the direction which a colliculus cell is most responsive to. #NEEDS TO BE FIXED _ CURRENTLY YANKED"""
	Tact = params["Tact"]
	base_a = params["retinal_activation"]["threshold"] .+ 0.1*slope
	frs = params["retinal_activation"]
	frw = params["colliculus_activation"]

	records = zeros(length(collicular_hash), length(retinal_hash))
	for (i, loc) in retinal_hash
		a = zeros(length(retinal_hash), Tact)
		u = zeros(length(colliculus_hash), Tact)
		a[i,:] .= base_a
		for t = 1:(Tact-1)
			u[:,t+1] = nft_step(a[:,t], u[:,t], W, S, tau_short, frw, frs)
		end
		for (j, loc) in colliculus_hash
			records[j, i] = u[j, end]
		end
	end

	maximal_ret_inds = map(x -> argmax(records[x, :]), 1:size(records[1]))

	color_map = map(x -> RGB(retinal_hash[x][1], 0.45, retinal_hash[x][2]), maximal_ret_inds)

	plt = plot(st=:scatter, DPI=400)
	for (j, loc) in collicular_hash
		plot!(plt, loc[1], loc[2], c=color_map[j])
	end
	
	savefig(img_dir, plt)
end
