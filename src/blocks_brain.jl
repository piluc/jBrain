global const BLOCKS = "BLOCKS"
global const FIRST = "FIRST"
global const SECOND = "SECOND"
global const ABOVE = "ABOVE"
global const RELATION = "RELATION"

# TODO: test if using vectors instead of sets can be faster
mutable struct BlocksBrain
    brain::Brain
    num_blocks::Int64
    all_areas::Vector{String}
    explicit_area_k::Int64
    non_explicit_area_n::Int64
    non_explicit_area_k::Int64
    fiber_states::Dict{String,Dict{String,Set{Int64}}}
    # In general, inhibition and disinibition behave like mutex locks, and are not just booleans.
    # For example, it is possible for two different assemblies to inhibit an area A, and two separate disinhibitions
    # are necessary to fully disinhibit A. 
    area_states::Dict{String,Set{Int64}}
    function BlocksBrain(blocks_number::Int64, other_areas, p=0.1, eak::Int64=10, nean::Int64=10000, neak::Int64=100, db=0.2) #TODO: can we make this constructor *outer*?
        b = Brain(p)
        all_areas = vcat([BLOCKS], other_areas)
        ean = blocks_number * eak
        add_explicit_area(b, BLOCKS, ean, eak, db)
        for area_name in other_areas
            add_area(b, area_name, nean, neak, db)
        end
        fs = Dict{String,Dict{String,Set{Int64}}}()
        for from_area in all_areas
            fs[from_area] = Dict{String,Set{Int64}}()
            for to_area in all_areas
                fs[from_area][to_area] = push!(Set{Int64}(), 0)
            end
        end
        as = Dict{String,Set{Int64}}()
        for area in all_areas
            as[area] = push!(Set{Int64}(), 0)
        end
        return new(b, blocks_number, all_areas, eak, nean, neak, fs, as)
    end
end

function inhibit_area(blocks_brain::BlocksBrain, area_name::String, lock::Int64)
    push!(blocks_brain.area_states[area_name], lock)
end

function inhibit_fiber(blocks_brain::BlocksBrain, area1, area2, lock)
    push!(blocks_brain.fiber_states[area1][area2], lock)
    push!(blocks_brain.fiber_states[area2][area1], lock)
end

function disinhibit_area(blocks_brain::BlocksBrain, area_name, lock)
    delete!(blocks_brain.area_states[area_name], lock)
end

function disinhibit_fiber(blocks_brain::BlocksBrain, area1, area2, lock)
    delete!(blocks_brain.fiber_states[area1][area2], lock)
    delete!(blocks_brain.fiber_states[area2][area1], lock)
end

function activate_block(blocks_brain::BlocksBrain, index)
    area = blocks_brain.brain.areas[BLOCKS]
    k = area.k
    assembly_start = ((index - 1) * k) + 1
    area.winners = collect(assembly_start:(assembly_start + k - 1))
    fix_assembly(area)
end

function activate_assembly(blocks_brain::BlocksBrain, index, activation_area)
    area = blocks_brain.brain.areas[activation_area]
    k = area.k
    assembly_start = ((index - 1) * k) + 1
    area.winners = collect(assembly_start:(assembly_start + k - 1))
    fix_assembly(area)
end

function get_project_map(blocks_brain::BlocksBrain; verbose=false)::Dict{String,Set{String}}
    project_map = Dict{String,Set{String}}()
    for area1 in blocks_brain.all_areas
        if (verbose) println(area1) end
        as1 = get(blocks_brain.area_states, area1, Set{Int64}())
        if (length(as1) == 0)
            if (verbose) println("In get_project_map as1 empty") end
            for area2 in blocks_brain.all_areas
                if (area1 != BLOCKS || area2 != BLOCKS)
                    if (verbose) println("In get_project_map area1=", area1, " and area2=", area2) end
                    as2 = get(blocks_brain.area_states, area2, Set{Int64}())
                    if (length(as2) == 0)
                        fs1 = get(blocks_brain.fiber_states, area1, Dict{String,Set{Int64}}())
                        if (length(get(fs1, area2, Set{Int64}())) == 0)
                            if (length(blocks_brain.brain.areas[area1].winners) > 0)
                                pma1 = get(project_map, area1, Set{String}())
                                if (length(pma1) == 0)
                                    project_map[area1] = Set([area2])
                                else
                                    push!(project_map[area1], area2)
                                end
                            end
                            if (length(blocks_brain.brain.areas[area2].winners) > 0)
                                pma2 = get(project_map, area2, Set{String}())
                                if (length(pma2) == 0)
                                    project_map[area2] = Set([area2])
                                else
                                    push!(project_map[area2], area2)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return project_map
end

function fix_assemblies_for_project(blocks_brain::BlocksBrain, project_map)
    for area in keys(project_map)
        if ((!(BLOCKS in keys(project_map)) || !(area in project_map[BLOCKS])) )
            fix_assembly(blocks_brain.brain.areas[area])
        elseif (area != BLOCKS)
            unfix_assembly(blocks_brain.brain.areas[area])
            blocks_brain.brain.areas[area].winners = []
        end
    end
end

function project_star(blocks_brain; project_rounds=30, verbose=false)
    project_map = get_project_map(blocks_brain)
    fix_assemblies_for_project(blocks_brain, project_map)
    # fix_assemblies_for_project might reset the winners of an area to empty, so recompute project_map.
    project_map = get_project_map(blocks_brain)
    if (verbose)
        println("Got project map: ", project_map)
    end
    for _ in 1:project_rounds
        project_map = get_project_map(blocks_brain)
        project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)
        if (verbose)
            project_map = get_project_map(blocks_brain)
            println("Got project map: ", project_map)
        end
    end
end

function test_project(blocks_brain; verbose=false)
    # Temporarily disable plasticity so this projection doesn't change connectomes.
    blocks_brain.brain.no_plasticity = true
    project_map = get_project_map(blocks_brain)
    fix_assemblies_for_project(blocks_brain, project_map)
    # fix_assemblies_for_project might reset the winners of an area to empty, so recompute project_map.
    project_map = get_project_map(blocks_brain)
    for area in keys(project_map)
        if (area in project_map[area])
            delete!(project_map[area], area)
        end
    end
    if (verbose)
        println("        Projecting ", project_map, " (in test_project)")
    end
    project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)
    blocks_brain.brain.no_plasticity = false
end

## next function converts a Set{Int64} object into a Vector{Int64} object

function set_to_vector(S::Set{Int64})::Vector{Int64}
    n = length(S)
    vec = Vector{Int64}(undef,n)
    count = 1
    for x in S
        vec[count] = x
        count = count + 1
    end
    return vec
end

function is_assembly(blocks_brain::BlocksBrain, area_name, min_overlap=0.75; verbose=false)
    if (verbose)
        println("        Executing is_assembly on ", area_name)
        println("        States of ", area_name, ": ", get(blocks_brain.area_states, area_name, Set{Int64}()))
    end
    if (length(get(blocks_brain.area_states, area_name, Set{Int64}())) > 0)
        return false
    end
    area = blocks_brain.brain.areas[area_name]
    if (length(area.winners) == 0)
        return false
    end
    blocks_brain.brain.no_plasticity = true
    winners_before = Set{Int64}(area.winners)
    threshold = min_overlap * area.k
    project_map = Dict{String,Set{String}}()
    project_map[area_name] = Set([area_name])
    project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)
    winners_after = Set{Int64}(area.winners)
    blocks_brain.brain.no_plasticity = false
    # Restore previous winners (so testing for assembly does not affect system).
    area.winners = set_to_vector(winners_before)
    if (length(intersect(winners_before, winners_after)) >= threshold)
        return true
    end
    return false
end

# TODO: get_block_index returns 0 when it does not find the block, so we should enforce that the input block list are strictly positive `Int`s. A natural way could be to have a struct with a constructor that enfoces that.
"""
    get_block_index(blocks, area...)

Returns the index of the active block assembly in `area_name`. 
If none of the `blocks` is (sufficiently) active, returns 0. 
"""
function get_block_index(blocks_brain::BlocksBrain, area_name, min_overlap=0.75)
    area = blocks_brain.brain.areas[area_name]
    (length(area.winners) > 0) || throw(ArgumentError(string("Cannot get index (block) because no assembly in ", area_name)))
    winners = Set{Int64}(area.winners)
    area_k = area.k
    threshold = min_overlap * area_k
    for block_index in 1:blocks_brain.num_blocks
        block_assembly_start = ((block_index - 1) * area_k) + 1
        block_assembly = Set(block_assembly_start:(block_assembly_start + area_k - 1))
        if (length(intersect(winners, block_assembly)) >= threshold)
            return block_index
        end
    end
    return 0
end