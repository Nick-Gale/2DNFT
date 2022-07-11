# Some important biological functions

function act(val; rate=100, slope=0.01, threshold=0.5)
    """A sigmoidal activation function returning spikes/s specified with a base rate, a firing threshold in mV, and a slope in mV."""
    return rate/(1 + exp(-(val - threshold)/slope))
end

function Wkernel(collicular_hash; R1=1, R2=0.5, r1=0.1, r2=0.35, rebufscale=0.5)
    """Specifies the recurrent connections as a difference of two Gaussians with units of syn mm^-2 and length scales in mm."""
    n = length(collicular_hash)
    mat = zeros(n, n)
    for i = 1:n
        for j = 1:n
            d2 = sum(collicular_hash[i] .^ 2 .+ collicular_hash[j] .^ 2) * rebufscale
            mat[i, j] = R1 * exp(-d2/2/r1^2) - R2 * exp(-d2/2/r2^2)
        end
    end
    return mat
end

function S0(retina_hash, colliculus_hash; S1=0.5, S2=0.25, s1=0.4, s2=0.6, buf=0.5)
    """Creates the adjacency matrix for each retinal site specified as a difference of two Gaussians. Assumes the retina is in register with colliculus and rescales according to major and minor axes."""
    maxy_ret = reduce((x,y) -> retina_hash[x][2] <= retina_hash[y][2] ? y : x, keys(retina_hash))
    maxx_ret = reduce((x,y) -> retina_hash[x][1] <= retina_hash[y][1] ? y : x, keys(retina_hash))

    maxy_col = reduce((x,y) -> colliculus_hash[x][2] <= colliculus_hash[y][2] ? y : x, keys(colliculus_hash))
    maxx_col = reduce((x,y) -> colliculus_hash[x][1] <= colliculus_hash[y][1] ? y : x, keys(colliculus_hash))
    adj = zeros(length(retina_hash), length(colliculus_hash))

    for (i, valret) in retina_hash
        for (j, valcol) in colliculus_hash
            d2 = sum((valret ./ (maxx_ret, maxy_ret) .* (maxx_col * buf, maxy_col * buf) .- valcol).^2)
            adj[i,j] =  S1 * exp(-d2/2/s1^2) - S2 * exp(-d2/2/s2^2)
        end
    end

    return adj
end

function plasticity_window(A0, tp, T; case="STDP")
    """The plasticity function with a window of tp milliseconds averaging T milliseconds of activity. STDP is specfied as case=1 and CDP as case=2."""
    if case == "STDP"
        i = 1
    elseif case == "CDP"
        i = 2
    else
        throw("""You need to specify a valid case: \"STDP\" or \"CDP\". These are strings.""")
    end

    iter = -T/2:T/2
    return sign.(iter) .^ i .* A0 .* exp.(-abs.(iter) ./ tp)
end
