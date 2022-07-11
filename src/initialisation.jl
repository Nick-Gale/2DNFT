function hash_cartesian(L, interior_boundary; params)
    """Create a linear index hash for location of active cells in the enviroment."""
    count = 0
    hash = Dict()
    for i = 1:L
        for j = 1:L
            x = (i - 1)/L - 0.5
            y = (j - 1)/L - 0.5
            if interior_boundary(x, y; params...)
                count+=1
                hash[count] = (x, y)
            end
        end
    end
    return hash
end

function hash_retina(x, y; radius=0.5)
    """Interior boundary function for hashing cartesian coordinates of retina. Assumes circular cross-section."""
    return sqrt(x^2 + y^2 < radius)
end

function hash_colliculus(x, y; xradius=0.5, yradius=0.25)
    """Interior boundary function for hashing cartesian coordinates of colliculus. Assumes colliclus has been rotated and rescaled."""
    return (x/xrad)^2 + (y/yrad)^2 < 1
end

function link_retina_colliculus(ret_mask, colliculus_mask, buf=0.5)
    """Link colliculus sites which can be activated by retinal stimulation. The linkage is all-to-all inside a buffered region of colliculus and 0 outside the buffer. The buffer with a rapidly decaying colliculus recurrent kernel removes the edge effects but is a computational convience - the buffer region interior boundary should be considered as the colliculus exterior boundary."""
    nret = length(ret_mask)
    ncol = length(colliculus_mask)
   
    maxy = reduce((x,y) -> colliculus_mask[x][2] <= colliculus_mask[y][2] ? y : x, keys(colliculus_mask))
    maxx = reduce((x,y) -> colliculus_mask[x][1] <= colliculus_mask[y][1] ? y : x, keys(colliculus_mask))

    Smask = zeros(ncol, nret)
    for j = 1:nret
        for i = 1:ncol
            Smask[i, j] = (colliculus_mask[i][1] < buf * maxx) && (colliculus_mask[i][2] < buf * maxy) 
        end
    end
    return Smask
end

function intialise(parameters; buf=0.5)
    """Initialise all the things required for the simulation"""
    max_ellipse_radius = maximum([initSdict[:collicular_x_radius], initSdict[:collicular_y_radius]])
    rebuffed_scale = max_ellipse_radius / buf # the multiplier to convert from distances in the collicular hash to real mm distances.

    retina_init = Dict(:radius => parameters[:retinal_radius])
    collicular_init = Dict(:xradius => parameters[:collicular_x_radius] / rebuffed_scale, :yradius => initSdict[:collicular_y_radius] / rebuffed_scale)
    S0init = Dict(:S1 => parameters[:S1], :S2 => parameters[:S2], :s1 => parameters[:s1], :s2 => parameters[:s2], :buf=>buf)
    Winit = Dict(:R1 => parameters[:R1], :R2 => parameters[:R2], :r1 => parameters[:r1], :r2 => parameters[:r2], :rebufscale=>rebuffed_scale)
        
    
    retinal_hash = hash_cartesian(L, hash_retina; retina_init...) 
    collicular_hash = hash_cartesian(L, hash_colliculus; collicular_init...) 

    Smask = link_retina_colliculus(retinal_hash, collicular_hash, buf)
    maskedS0 = S0(retinal_hash, collicular_hash; S0init...) .* Smask

    W = Wkernel(collicular_hash; Winit...)

    return rebuffed_scale, retinal_hash, collicular_hash, Smask, maskedS0, W, H
end