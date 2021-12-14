using Distributions
using PaddedViews
using Logging

verboseLevel = 1
logger = Logging.SimpleLogger(stderr, verboseLevel)
Logging.global_logger(logger)

"""
    Stimulus(k::Int64)

# Examples

```jldoctest
julia> s = Stimulus(10);

julia> s.k == 10
true
```
"""
mutable struct Stimulus
    k::Int64
end

"""
    Area(name...)
    Area(name, n, k, beta=0.05)

Area class. Arguments are described below. 
`Area(name, n, k, beta=0.05)` initializes an area with all omitted argument empty. 

# Arguments

- `name::String`: Area ID
- `n::Int64`: Number of neurons
- `k::Int64`: Number of winners
- `beta::Float64`: Default plasticity
- `stimulus_beta::Dict{String,Float64}`: Plasticity for connections from stimuli to this area
- `area_beta::Dict{String,Float64}`: Plasiticity for connections from areas to this area
- `w::Int64`: Number of winners
- `winners::Vector{Int64}`: List of winners
- `new_w::Int64`: Number of new winners
- `new_winners::Vector{Int64}`: List of new winners
- `saved_winners::Vector{Vector{Int64}}`: List of lists of winners at each round
- `saved_w::Vector{Int64}`: List of numbers of winners at each round
- `num_first_winners::Int64`: Initial number of winners
- `fixed_assembly::Bool`: Whether to fix (freeze) the assembly (winners) in this area
- `fire`: Whether to fully simulate this area: if the area is explicit all neurons
- `explicit::Bool`: Whether to run the simulation explicitly. 
# TODO: the concept of explicit simulation should be explained somewhere...

# Examples

```jldoctest
julia> s = Stimulus(10);

julia> s.k == 10
true
```
"""
mutable struct Area
    name::String
    n::Int64
    k::Int64
    beta::Float64
    stimulus_beta::Dict{String,Float64}
    area_beta::Dict{String,Float64}
    w::Int64
    winners::Vector{Int64}
    new_w::Int64
    new_winners::Vector{Int64}
    saved_winners::Vector{Vector{Int64}}
    saved_w::Vector{Int64}
    num_first_winners::Int64
    fixed_assembly::Bool
    explicit::Bool
    function Area(name::String, n::Int64, k::Int64, beta::Float64=0.05)
        return new(name, n, k, beta, Dict{String,Float64}(), Dict{String,Float64}(), 0, [], 0, [], [], [], 0, false, false)
    end
end

function update_winners(a::Area)
    a.winners = a.new_winners
    if (!a.explicit)
        a.w = a.new_w
    end
end

function update_stimulus_beta(a::Area, stim_name::String, new_beta::Real)
    a.stimulus_beta[stim_name] = new_beta
end

function update_area_beta(a::Area, from_area_name::String, new_beta::Real)
    a.area_beta[from_area_name] = new_beta
end

function fix_assembly(a::Area)
    (length(a.winners) > 0) || throw(ArgumentError("Area $(a.name) does not have assembly (cannot fix)."))
    a.fixed_assembly = true 
end

function unfix_assembly(a::Area)
    a.fixed_assembly = false
end

"""
    Brain(areas...)

The brain. 

# Parameters

* `areas`: Dictionary of areas in the brain
* `stimuli`: Dictionary of stimuli in the brain
* `stimuli_connectomes`: For each stimulus in the brain, dictionary of connections from it to each area of the brain
* `p`: probability parameter for the edges of the underlying Erdos-Renyi graph 

# Example 
```jldoctest
julia> b = Brain(.1); 
julia> b.areas
Dict{String, Area}()
```
"""
mutable struct Brain
    areas::Dict{String,Area}
    stimuli::Dict{String,Stimulus}
    stimuli_connectomes::Dict{String,Dict{String,Vector{Float64}}}
    connectomes::Dict{String,Dict{String,Matrix{Float64}}}
    p::Float64
    save_size::Bool
    save_winners::Bool
    no_plasticity::Bool
    # Constructor
    function Brain(p, save_size=true, save_winners=false)
        return new(Dict{String,Area}(), Dict{String,Stimulus}(), Dict{String,Dict{String,Vector{Float64}}}(), Dict{String,Dict{String,Matrix{Float64}}}(), p, save_size, save_winners, false)
    end
end


# Add stimulus to brain. The weights of the connections to an explicit area of the brain are set randomly 
# according to a binomial distribution with parameters k and the value p of the brain. The plasticity of
# these connections are set equal to the default plasticity of the area.
function add_stimulus(b::Brain, stim_name::String, k::Int)
    b.stimuli[stim_name] = Stimulus(k)
    new_connectomes = Dict{String,Vector{Float64}}()
    for (area_name, area) in b.areas
        if (area.explicit)
            new_connectomes[area_name] = rand(Binomial(k, b.p),  area.n) .* 1.0
        else
            new_connectomes[area_name] = Vector{Float64}()
        end
        area.stimulus_beta[stim_name] = area.beta
    end
    b.stimuli_connectomes[stim_name] = new_connectomes
end

# Add a non-explicit area to the brain. Since the area is not explicit, all connections from the
# stimuli of the brain and from/to all the areas of the brain are initially set to empty (they
# will be set during the project phase in order to improve performance).
function add_area(b::Brain, area_name::String, n::Int, k::Int, beta::Real)
    b.areas[area_name] = Area(area_name, n, k, beta)
    for (stim_name, stim_connectomes) in b.stimuli_connectomes
        stim_connectomes[area_name] = Vector{Float64}()
        b.areas[area_name].stimulus_beta[stim_name] = beta
    end
    new_connectomes = Dict{String,Matrix{Float64}}()
    for (other_area_name, other_area) in b.areas
        other_area_size = 0
        if (other_area.explicit)
            other_area_size = other_area.n
        end
        new_connectomes[other_area_name] = Matrix{Float64}(undef, 0, other_area_size)
        if (other_area_name != area_name)
            b.connectomes[other_area_name][area_name] = Matrix{Float64}(undef, other_area_size, 0)
        end
        other_area.area_beta[area_name] = other_area.beta
        b.areas[area_name].area_beta[other_area_name] = beta
    end
    b.connectomes[area_name] = new_connectomes
end

# Add an explicit area to the brain. Since the area is explicit, the weights of all connections 
# from a stimuli of the brain the new area are initially set randomly according to a 
# binomial distribution with parameters the value k of the stimulus and the value p of the brain
# (with the default plasticity). The weights of all connections from/to an explicit area of the brain
# to the new area are initially and fully set randomly according to a binomial distribution 
# with parameters 1 and the value p of the brain. The weights of all connections from/to a 
# non-explicit area of the brain to the new area are initially set to empty. 
# In all cases, the plasticity of the connections is set to the default plasticity.
# The number of winners of the new area is set equal to the number of its neurons.
function add_explicit_area(b::Brain, area_name::String, n::Int, k::Int, beta::Real)
    b.areas[area_name] = Area(area_name, n, k, beta)
    b.areas[area_name].explicit = true
    for (stim_name, stim_connectomes) in b.stimuli_connectomes
        stim_connectomes[area_name] = rand(Binomial(stimuli[stim_name].k, b.p),  n) .* 1.0
        b.areas[area_name].stimulus_beta[stim_name] = beta
    end
    new_connectomes = Dict{String,Matrix{Float64}}()
    for (other_area_name, other_area) in b.areas
        if (other_area_name == area_name)
            new_connectomes[other_area_name] = rand(Binomial(1, b.p),  n, n) .* 1.0
        else
            if (other_area.explicit)
                other_n = other_area.n
                new_connectomes[other_area_name] = rand(Binomial(1, b.p),  n, other_n) .* 1.0
                b.connectomes[other_area_name][area_name] = rand(Binomial(1, b.p), other_n,  n) .* 1.0
            else
                new_connectomes[other_area_name] = Matrix{Float64}(undef, n, 0)
                b.connectomes[other_area_name][area_name] = Matrix{Float64}(undef, 0, n)
            end
        end
        other_area.area_beta[area_name] = other_area.beta
        b.areas[area_name].area_beta[other_area_name] = beta
    end
    b.connectomes[area_name] = new_connectomes
    b.areas[area_name].w = n
end

# add_explicit_area(b::Brain, area_name::String, n::Int, k::Int, beta::Real) = add_explicit_area(b, area_name, Int64(n), Int64(k), beta)

# Update the plasticities of the connections between stimuli and areas. Each area update consists of
# of the destination area and a list of update rules: each update rule specifies the source area
# and the new plasticity. Each stimulus update consists of the destination area and a list of update rules:
# each update rule specifies the source stimulus and the new plasticity. 
function update_plasticities(b::Brain, area_update_map::Dict{String,Dict{String,Float64}}=Dict{String,Dict{String,Float64}}(), stim_update_map::Dict{String,Dict{String,Float64}}=Dict{String,Dict{String,Float64}}())
    for (to_area, update_rules) in area_update_map
        for (from_area, new_beta) in update_rules
            b.areas[to_area].area_beta[from_area] = new_beta
        end
    end
    for (to_area, update_rules) in stim_update_map
        for (stim, new_beta) in update_rules
            b.areas[to_area].stimulus_beta[stim] = new_beta
        end
    end
    b
end

# Compute the potential new k winners of the area to which we are going to project.
# To this aim compute the threshold alpha for inputs that are above (n-k)/n percentile, 
# use the normal approximation, between alpha and total_k, round to integer, and 
# create the k potential_new_winners.
function compute_potential_new_winners(b::Brain, to_area::Area, total_k::Int)
    effective_n = to_area.n - to_area.w
    alpha = quantile(Binomial(total_k, b.p), (effective_n - to_area.k) / effective_n)
    @info alpha
    std = sqrt(total_k * b.p * (1.0 - b.p))
    mu = total_k * b.p
    l = (alpha - mu) / std
    r = (total_k - mu) / std
    @info l r mu std
    potential_new_winners = rand(TruncatedNormal(0, 1, l, r), to_area.k) .* std
    for i in 1:to_area.k
        potential_new_winners[i] = potential_new_winners[i] + mu
        potential_new_winners[i] = round(potential_new_winners[i], digits=0)
    end
    @info length(potential_new_winners)
    @info maximum(potential_new_winners)
    @info potential_new_winners
    return potential_new_winners
end

# TODO: the following should be refactored
function project_into(brain::Brain, to_area::Area, from_stimuli::Set{String}, from_areas::Set{String})
    @info "Projecting" from_stimuli from_areas to_area.name
    for from_area in from_areas
        @assert ((length(brain.areas[from_area].winners) != 0) && (brain.areas[from_area].w != 0)) "Projecting from area ($from_area)  with no winners."
    end
    # Compute previous inputs from the winners of the areas from which we project to the
    # current w winners of the area to which we project (indexed from 1 to w). In particular, 
    # for each winner i of the area to which we project, its total input is computed by summing 
    # the weights of the connectomes connecting either a stimulus or a winner of an area from 
    # which we project to it.
    name = to_area.name
    prev_winner_inputs = zeros(to_area.w)
    for stim in from_stimuli
        stim_inputs = brain.stimuli_connectomes[stim][name]
        for i in 1:to_area.w
            prev_winner_inputs[i] = prev_winner_inputs[i] + stim_inputs[i]
        end
    end
    for from_area in from_areas
        connectome = brain.connectomes[from_area][name]
        for w in brain.areas[from_area].winners
            for i in 1:to_area.w
                prev_winner_inputs[i] = prev_winner_inputs[i] + connectome[w,i]
            end    
        end
    end
    @info length(prev_winner_inputs) maximum(prev_winner_inputs; init=-Inf) prev_winner_inputs
    # Simulate to_area.k potential new winners if the area is not explicit.
    if (!to_area.explicit)
        # Compute the number of input stimuli and areas, the total number of input connectomes,
        # and the number of input connectomes for each input stimulus and area.
        total_k = 0
        input_sizes = []
        num_inputs = 0
        for stim in from_stimuli
            total_k = total_k + brain.stimuli[stim].k
            push!(input_sizes, brain.stimuli[stim].k)
            num_inputs = num_inputs + 1
        end
        for from_area in from_areas
            # @assert (length(brain.areas[from_area].winners) == brain.areas[from_area].k) "The number of winners is different from k"
            effective_k = length(brain.areas[from_area].winners)
            total_k = total_k + effective_k
            push!(input_sizes, effective_k)
            num_inputs = num_inputs + 1
        end
        @assert (num_inputs == length(from_stimuli) + length(from_areas)) "The number of inputs should be the sum of winners in source stimuli and areas"
        @info total_k input_sizes
        # Compute the potential new k winners of the area to which we are going to project.
        # To this aim compute the threshold alpha for inputs that are above (n-k)/n percentile, 
        # use the normal approximation, between alpha and total_k, round to integer, and 
        # create the k potential_new_winners.
        effective_n = to_area.n - to_area.w
        alpha = quantile(Binomial(total_k, brain.p), (effective_n - to_area.k) / effective_n)
        @info alpha
        std = sqrt(total_k * brain.p * (1.0 - brain.p))
        mu = total_k * brain.p
        a = (alpha - mu) / std
        b = (total_k - mu) / std
        @info a b mu std
        potential_new_winners = rand(TruncatedNormal(0, 1, a, b), to_area.k) .* std
        for i in 1:to_area.k
            potential_new_winners[i] = potential_new_winners[i] + mu
            potential_new_winners[i] = round(potential_new_winners[i], digits=0)
        end
        @info length(potential_new_winners) maximum(potential_new_winners; init=-Inf) potential_new_winners
        all_potential_winners = vcat(prev_winner_inputs, potential_new_winners)
    else
        all_potential_winners = prev_winner_inputs
    end
    @info all_potential_winners
    # At the end I have decided to adopt the simplest solution using sortperm and not partialsortperm!.
    # However, the two solutions should be equivalent (maybe the one using partialsortperm! is more efficient).
    # new_winner_indices = collect(partialsortperm!(collect(1:length(all_potential_winners)), all_potential_winners, 1:to_area.k, rev=true, initialized=true))
    new_winner_indices = sortperm(all_potential_winners, rev=true)[1:to_area.k]
    @info new_winner_indices
    num_first_winners = 0
    if (!to_area.explicit)
        first_winner_inputs = []
        for i in 1:to_area.k
            if (new_winner_indices[i] > to_area.w)
                push!(first_winner_inputs, potential_new_winners[new_winner_indices[i] - to_area.w])
                num_first_winners = num_first_winners + 1
                new_winner_indices[i] = to_area.w + num_first_winners
            end
        end
    end
    @info num_first_winners
    to_area.new_winners = new_winner_indices
    to_area.new_w = to_area.w + num_first_winners
    if (to_area.fixed_assembly)
        to_area.new_winners = to_area.winners
        to_area.new_w = to_area.w
        first_winner_inputs = []
        num_first_winners = 0
    end
    @info to_area.new_winners to_area.new_w
    first_winner_to_inputs = Dict{Int64,Vector{Int64}}()
    for i in 1:num_first_winners
        input_indices = sample(1:total_k, Int64(first_winner_inputs[i]), replace=false)
        inputs = zeros(num_inputs)
        total_so_far = 0
        for j in 1:num_inputs
            inputs[j] = 0
            for w in input_indices
                if ((total_so_far + input_sizes[j]) > (w - 1) >= total_so_far)
                    inputs[j] = inputs[j] + 1
                end
            end
            # inputs[j] = sum([((total_so_far + input_sizes[j]) > (w - 1) >= total_so_far) for w in input_indices])
            total_so_far = total_so_far + input_sizes[j]
        end
        first_winner_to_inputs[i] = inputs
    end
    m = 1
    for stim in from_stimuli
        if (num_first_winners > 0)
            append!(brain.stimuli_connectomes[stim][name], zeros(num_first_winners))
            for i in 1:num_first_winners
                brain.stimuli_connectomes[stim][name][to_area.w + i] = first_winner_to_inputs[i][m]
            end
            # brain.stimuli_connectomes[stim][name] = append!(brain.stimuli_connectomes[stim][name], first_winner_to_inputs[:][m])
        end
        stim_to_area_beta = to_area.stimulus_beta[stim]
        if (brain.no_plasticity)
            stim_to_area_beta = 0.0
        end
        for i in to_area.new_winners
            brain.stimuli_connectomes[stim][name][i] = brain.stimuli_connectomes[stim][name][i] * (1 + stim_to_area_beta)
        end
        @info stim name brain.stimuli_connectomes[stim][name]
        m = m + 1
    end
    for from_area in from_areas
        from_area_w = brain.areas[from_area].w
        from_area_winners = brain.areas[from_area].winners
        nr = size(brain.connectomes[from_area][name])[1]
        nc = size(brain.connectomes[from_area][name])[2]
        # It seems to me that the third solution is the fastest one
        @info "Executing padding of connectomes from $from_area to $name by adding $num_first_winners columns" from_area name num_first_winners 
        # brain.connectomes[from_area][name] = collect(PaddedView(0, brain.connectomes[from_area][name], (nr, nc + num_first_winners), (1, 1)))
        # brain.connectomes[from_area][name] = hcat(brain.connectomes[from_area][name], zeros(nr, num_first_winners))
        brain.connectomes[from_area][name] = reduce(hcat, [brain.connectomes[from_area][name], zeros(nr, num_first_winners)])

        for i in 1:num_first_winners
            total_in = first_winner_to_inputs[i][m]
            sample_indices = sample(from_area_winners, Int64(total_in), replace=false)
            for j in 1:from_area_w
                if (j in sample_indices)
                    brain.connectomes[from_area][name][j,to_area.w + i] = 1.0
                end
                if (!(j in from_area_winners))
                    brain.connectomes[from_area][name][j,to_area.w + i] = rand(Binomial(1, brain.p))
                end
            end
        end
        area_to_area_beta = to_area.area_beta[from_area]
        if brain.no_plasticity
            area_to_area_beta = 0.0
        end
        @info "Plasticity in projecting $from_area -> $name is $area_to_area_beta" from_area name area_to_area_beta
        for i in to_area.new_winners
            for j in from_area_winners
                brain.connectomes[from_area][name][j,i] = brain.connectomes[from_area][name][j,i] * (1.0 + area_to_area_beta)
            end
        end
        @info from_area name brain.connectomes[from_area][name]
        m = m + 1
    end
    for other_area in keys(brain.areas)
        if (!(other_area in from_areas))
            brain.connectomes[other_area] = get(brain.connectomes, other_area, Dict{String,Matrix{Float64}}())
            con = get(brain.connectomes[other_area], name, Array{Float64,2}(undef, 0, 0))
            nr = size(con)[1]
            nc = size(con)[2]
            # It seems to me that the third solution is the fastest one
            @info "Executing padding of connectomes from $other_area to $name by adding $num_first_winners columns" other_area name num_first_winners
            # brain.connectomes[other_area][name] = collect(PaddedView(0, con, (nr, nc + num_first_winners), (1, 1)))
            # brain.connectomes[other_area][name] = hcat(con, zeros(nr, num_first_winners)
            brain.connectomes[other_area][name] = reduce(hcat, [con, zeros(nr, num_first_winners)])
            for j in 1:brain.areas[other_area].w
                for i in (to_area.w + 1):to_area.new_w
                    brain.connectomes[other_area][name][j,i] = rand(Binomial(1, brain.p))
                end
            end
        end
        brain.connectomes[name] = get(brain.connectomes, name, Dict{String,Matrix{Float64}}())
        cno = get(brain.connectomes[name], other_area, Array{Float64,2}(undef, 0, 0))
        nr = size(cno)[1]
        nc = size(cno)[2]
        # It seems to me that the third solution is the fastest one
        @info "Executing padding of connectomes from $name to $other_area by adding $num_first_winners rows" name other_area num_first_winners
        # brain.connectomes[name][other_area] = collect(PaddedView(0, cno, (nr + num_first_winners, nc), (1, 1)))
        # brain.connectomes[name][other_area] = vcat(cno, zeros(num_first_winners, nc))
        brain.connectomes[name][other_area] = reduce(vcat, [cno, zeros(num_first_winners, nc)])
        @info name other_area size(brain.connectomes[name][other_area])
        columns = size(brain.connectomes[name][other_area])[2]
        @info to_area.w to_area.new_w
        for i in (to_area.w + 1):to_area.new_w
            for j in 1:columns
                brain.connectomes[name][other_area][i,j] = rand(Binomial(1, brain.p))
            end
        end
        @info name other_area brain.connectomes[name][other_area]
    end
    return num_first_winners
end

# TODO: docstringify the following and add assumptions corresponding to thrown errors

# Execute the project from stimuli and/or areas to areas. For each stimulus (key) in the first dictionary,
# the list (value) of areas to which the stimulus has the project is specified. For each area (key),
# in the second dictionary, the list (value) of areas to which the area has the project is specified.
# The function collects, for each area, the list of stimuli and areas that project to it (basically, it
# computes the inverse of the input mappings). Then, for each area which has to "receive" a projection
# (from either stimuli or areas), it invokes the function which actually performs the projection (this
# function returns the number of winners in the destination area). If the new winners have to be later
# analysed, then their list is appended to the the list of lists of winners of the area. When everything
# has been done, the function updates the destination areas.
function project(b::Brain, stim_to_area::Dict{String,Set{String}}, area_to_area::Dict{String,Set{String}}) 
    @info stim_to_area area_to_area
    stim_in = Dict{String,Set{String}}()
    area_in = Dict{String,Set{String}}()
    for (stim, areas) in stim_to_area
        (stim in keys(b.stimuli)) || throw(ArgumentError(stim * " is not in the stimuli of the brain"))
        for area in areas
            (area in keys(b.areas)) || throw(ArgumentError(area * " is not in the areas of the brain"))
            if (!(area in keys(stim_in)))
                stim_in[area] = Set([stim])
            else
                push!(stim_in[area], stim)
            end
        end
    end
    # TODO: Perhaps error should be as above (ArgumentError), and test should be added 
    for (from_area, to_areas) in area_to_area
        if (!(from_area in keys(b.areas)))
            error(from_area, " is not in the areas of the brain")
            return    
        end
        for to_area in to_areas
            if (!(to_area in keys(b.areas)))
                error(area, " is not in the areas of the brain")
                return    
            end
            if (!(to_area in keys(area_in)))
                area_in[to_area] = Set([from_area])
            else
                push!(area_in[to_area], from_area)
            end
        end
    end
    to_update = union(Set(keys(stim_in)), Set(keys(area_in)))
    @info to_update
    for area in to_update
        num_first_winners = project_into(b, b.areas[area], get(stim_in, area, Set{String}()), get(area_in, area, Set{String}()))
        b.areas[area].num_first_winners = num_first_winners
        if (b.save_winners)
            push!(b.areas[area].saved_winners, b.areas[area].new_winners)
        end
    end
    for area in to_update
        update_winners(b.areas[area])
        if (b.save_size)
            push!(b.areas[area].saved_w, b.areas[area].w)
        end
    end
end
