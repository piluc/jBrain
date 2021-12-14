import Base.intersect

global TMP = "TEMPORARY"

# Look for the common part of two stacks. Starting from the bottom block of
# each stack (which is determined by using the function find_last_active_area_*),
# the function scans the two stacks until the top of one of them is reached
# or two different blocks are found. While doing so, it maintains the top
# block of the common part found so far. If at the end, this block exists,
# the function stores this block in a temporary area and return true, otherwise it 
# returns false. In order to check whether the top of a stack is reached,
# the function checks whether the correpsonding head area is activated.
## activates the assembly in NODE corresponding to the last block in the stack
# Parameters
# - blocks_brain: type BlocksBrain, after having parsed the stacks
# - source_stack_index: the index of the source stack (integer)
# - prefix_source: string "I", "G", or "B" whether the areas we work on are the input areas, the goal areas or the building areas.
# - target_stack_index: the index of the target stack (integer)
# - prefix_target: string "I", "G", or "B" whether the areas we work on are the input areas, the goal areas or the building areas.
# Example: intersect!(blocks_brain, 1, "I", 3, "G")

function intersect(blocks_brain::BlocksBrain, source_stack_index::Int, prefix_source::String, target_stack_index::Int, prefix_target::String; verbose::Bool=false)::Tuple{Bool,Int64}

    # adding prefix to the regions we work with
    source_regions = add_prefix(vcat(REGIONS...), prefix_source)
    target_regions = add_prefix(vcat(REGIONS...), prefix_target)

    # defining locally the regions involved with stack j
    source_head = add_prefix([HEADS[source_stack_index]], prefix_source)[1]
    source_nodes = add_prefix(NODES[source_stack_index], prefix_source)
    target_head = add_prefix([HEADS[target_stack_index]], prefix_target)[1]
    target_nodes = add_prefix(NODES[target_stack_index], prefix_target)

    # during this process each HEAD will ``forget'' the last active assembly: we save it and then re-project it at the end
    # this is a workaround
    save_heads_to_TMP!(blocks_brain, [source_head,target_head])

    top_index_source, top_block_source = top(blocks_brain, source_stack_index, prefix_source)
    top_index_target, top_block_target = top(blocks_brain, target_stack_index, prefix_target)

    if (verbose)
        println("The top of the source stack is block ", top_block_source, ".")
        println("The top of the target stack is block ", top_block_target, ".")
    end

    source_stack_length, last_index_source = fire_last_assembly(blocks_brain, source_head, source_nodes, top_index_source)
    target_stack_length, last_index_target = fire_last_assembly(blocks_brain, target_head, target_nodes, top_index_target)

    blocks_brain.brain.no_plasticity = true
    for area in vcat([BLOCKS], source_regions, target_regions)
        unfix_assembly(blocks_brain.brain.areas[area])
    end

    is_source_top = false
    is_target_top = false
    blocks_equal = true
    current_index_source = last_index_source
    current_index_target = last_index_target

    last_equal_block = 0
    while (!is_source_top && !is_target_top && blocks_equal)
        next_index_source = mod(current_index_source - 1, 1:MAX_NODES_AREAS)
        next_index_target = mod(current_index_target - 1, 1:MAX_NODES_AREAS)

        disinhibit_area(blocks_brain, BLOCKS, 0)
        disinhibit_areas(blocks_brain, [source_nodes[current_index_source],source_nodes[next_index_source], source_head], 0)
        disinhibit_areas(blocks_brain, [target_nodes[current_index_target],target_nodes[next_index_target], target_head], 0)

        disinhibit_fibers(blocks_brain, [[source_nodes[current_index_source],source_nodes[next_index_source]],[source_nodes[current_index_source],BLOCKS], [source_nodes[current_index_source],source_head] ], 0)
        disinhibit_fibers(blocks_brain, [[target_nodes[current_index_target],target_nodes[next_index_target]],[target_nodes[current_index_target],BLOCKS], [target_nodes[current_index_target],target_head] ], 0)

        project_map = Dict{String,Set{String}}()
        project_map[source_nodes[current_index_source]] = Set{String}([source_nodes[next_index_source],source_head,BLOCKS])
        project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)

        source_block = get_block_index(blocks_brain, BLOCKS)
        is_source_top = is_assembly(blocks_brain, source_head)

        project_map = Dict{String,Set{String}}()
        project_map[target_nodes[current_index_target]] = Set{String}([target_nodes[next_index_target],target_head,BLOCKS])
        project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)

        target_block = get_block_index(blocks_brain, BLOCKS)
        is_target_top = is_assembly(blocks_brain, target_head)

        if (verbose)
            println("Coomparing source block ", source_block, " and target block ", target_block, ".")
        end
        blocks_equal = source_block == target_block

        if (blocks_equal)
            last_equal_block = source_block
            if (verbose)
                println("Block ", source_block, " is shared.")
            end
        end

        current_index_source = next_index_source
        current_index_target = next_index_target
        inhibit_areas(blocks_brain, [source_nodes[current_index_source],source_nodes[next_index_source], source_head], 0)
        inhibit_areas(blocks_brain, [target_nodes[current_index_target],target_nodes[next_index_target], target_head], 0)
        inhibit_fibers(blocks_brain, [[source_nodes[current_index_source],source_nodes[next_index_source]],[source_nodes[current_index_source],BLOCKS], [source_nodes[current_index_source],source_head] ], 0)
        inhibit_fibers(blocks_brain, [[target_nodes[current_index_target],target_nodes[next_index_target]],[target_nodes[current_index_target],BLOCKS], [target_nodes[current_index_target],target_head] ], 0)
    end

    # we re-proejct the last active assembly in each HEAD
    for head in [source_head,target_head]
        reactivate_head_from_TMP!(blocks_brain, head)
    end

    blocks_brain.brain.no_plasticity = false
    inhibit_areas(blocks_brain, vcat([BLOCKS], source_regions, target_regions), 0)
    inhibit_all_fibers(blocks_brain, vcat([BLOCKS], source_regions, target_regions), 0)


    if (last_equal_block > 0)
        disinhibit_areas(blocks_brain, [BLOCKS,TMP], 0)
        disinhibit_fiber(blocks_brain, BLOCKS, TMP, 0)
        activate_block(blocks_brain, last_equal_block)
        proj_star_rounds = 50
        project_star(blocks_brain; project_rounds=proj_star_rounds)
        inhibit_areas(blocks_brain, [BLOCKS,TMP], 0)
        inhibit_fiber(blocks_brain, BLOCKS, TMP, 0)
        return (true, last_equal_block)
    end
    return (false, last_equal_block)
end


### saves the assembly in head to TMP
# Parameters
# - blocks_brain: type BlocksBrain, after having parsed the stacks
# - heads: vector of strings, the head regions to save in TMP
# Example: save_heads_to_TMP!(blocks_brain, [HEADS[1],HEADS[2]])
function save_heads_to_TMP!(blocks_brain::BlocksBrain, heads::Vector{String}; verbose::Bool=false)::Nothing
    blocks_brain.brain.no_plasticity = false
    for area in vcat(heads, [TMP])
        unfix_assembly(blocks_brain.brain.areas[area])
    end
    disinhibit_areas(blocks_brain, vcat(heads, [TMP]), 0)
    for head in heads
        disinhibit_fiber(blocks_brain, TMP, head, 0)
    end
    project_map = Dict{String,Set{String}}()
    for head in heads
        project_map[head] = Set{String}([head, TMP])
    end
    project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)
    proj_rounds = 50
    project_star(blocks_brain;project_rounds=proj_rounds)
    inhibit_areas(blocks_brain, vcat(heads, [TMP]), 0)
    for head in heads
        inhibit_fiber(blocks_brain, TMP, head, 0)
    end
    return    
end

# only after save_heads_to_TMP!
# reactivates the assembly in head from TMP
# Parameters
# - blocks_brain: type BlocksBrain, after having parsed the stacks
# - head: a string, the head region where we want the assembly to be active
# Example: save_heads_to_TMP!(blocks_brain, HEADS[1])
function reactivate_head_from_TMP!(blocks_brain::BlocksBrain, head::String; verbose::Bool=false)::Nothing
    blocks_brain.brain.no_plasticity = true
    for area in [head,TMP]
        unfix_assembly(blocks_brain.brain.areas[area])
    end
    disinhibit_areas(blocks_brain, [head,TMP], 0)
    disinhibit_fiber(blocks_brain, TMP, head, 0)
    project_map = Dict(TMP => Set([head, TMP]))

    project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)
    inhibit_areas(blocks_brain, [head,TMP], 0)
    inhibit_fiber(blocks_brain, TMP, head, 0)
    return
end


function planning_optimized(blocks_brain::BlocksBrain, input_stacks_number::Int, input_prefix::String, goal_stacks_number::Int, goal_prefix::String, building_prefix::String; to_readout::Bool = false, verbose::Bool=false)::Tuple{Vector{Vector{String}},Vector{Vector{String}}}
    dismantling_actions = Vector{Vector{String}}(undef, input_stacks_number)
    reassembling_actions = Vector{Vector{String}}(undef, goal_stacks_number)
    starting_blocks = zeros(Int, goal_stacks_number) ## this area is a workaround. To be correct, we should save each block in a temporary area from which to start reassembling.
    for i in 1:input_stacks_number
        input_stack_paired = false
        for j in 1:goal_stacks_number
            non_empty_intersect, block = intersect(blocks_brain, i, input_prefix, j, goal_prefix; verbose=false)
            input_stack_paired = input_stack_paired || non_empty_intersect
            if (non_empty_intersect)
                starting_blocks[j] = block
                if (verbose)
                    println("Input stack ", i, " and goal stack ", j, " intersect until block ", block, " included.")
                end
                dismantling_actions[i] = dismantle_stack_up_to_block!(blocks_brain, i, input_prefix, block,input_stacks_number; to_readout, verbose=false)
                # reassembling_actions[j] = [string("Fix the remaining of input stack ", i, " as final stack ", j, ".")]
                reassembling_actions[j] = [string(block)]
                save_up_to_block!(blocks_brain, j, goal_prefix, block, building_prefix; verbose=false)
                break
            end
        end
        if !(input_stack_paired)
            if (verbose)
                println("Input stack ", i, " does NOT intersect with ANY goal stack.")
            end
            dismantling_actions[i] = dismantle_stack!(blocks_brain, i, input_prefix, input_stacks_number; to_readout, verbose=false)
        end
    end
    if (verbose)
        println("Input stacks dismantled.")
    end
    for j in 1:goal_stacks_number
        block = starting_blocks[j]
        if (block != 0)
            if (verbose)
                println("Goal stack ", j, " is built until block ", block, " included.")
            end
            append!(reassembling_actions[j], reassemble_stack_from_block!(blocks_brain, j, goal_prefix, building_prefix, block, input_stacks_number,goal_stacks_number; to_readout, verbose=false))
        else
            if (verbose)
                println("Goal stack ", j, " is NOT partially built.")
            end
            reassembling_actions[j] = reassemble_stack!(blocks_brain, j, goal_prefix, building_prefix,input_stacks_number, goal_stacks_number; to_readout, verbose=false)
        end
    end
    if (verbose)
        println("Blocks reassembled. The actions to perform follow.")
        for j in 1:goal_stacks_number
            println("Stack number ", j, ":")
            println(reassembling_actions[j])
        end
    end
    return (dismantling_actions, reassembling_actions)
end



## dismantle stack and put each block on the table up to a given block
# Parameters
# - blocks_brain: type BlocksBrain, after having parsed the stacks
# - stack_index: an integer, it is the index of the stack
# - prefix_source: string "I", "G", "B", or "T" whether the areas we work on are the input areas, the goal areas, the building areas, or the table areas.
# - block_to_keep: Int, the block where to stop (we stop by deleting the one before)
# - stacks_number: an integer, it is the number of input stacks (only for readout purposes)

# Example: dismantle_stack_up_to_block!(blocks_brain, 5, "I", 4)
function dismantle_stack_up_to_block!(blocks_brain::BlocksBrain, stack_index::Int, prefix_source::String, block_to_keep::Int, stacks_number::Int; to_readout::Bool = false, verbose::Bool=false)::Vector{String}
    stack_actions = Vector{String}()
    node_index, block = top(blocks_brain, stack_index, prefix_source)
    if block == block_to_keep
        if is_last_block(blocks_brain, stack_index, node_index, block, prefix_source)
            #pop!(blocks_brain,stack_index,prefix_source)
            put!(blocks_brain, stack_index, block, "T")
            if (verbose)
                println("Saving block ", block, " on the table.")
            end
            if (to_readout)
                read_block_world_instance(blocks_brain,stacks_number,prefix_source)
                read_block_world_instance(blocks_brain,stacks_number,"T")
            end
        end
        return stack_actions
    end

    while block != block_to_keep
        if (verbose)
            println("Now popping block ", block, " in area NODE", node_index, ".")
        end
        old_block = block
        # append!(stack_actions, [string("Remove block ", old_block, " from top of stack ", stack_index, ".")])
        append!(stack_actions, [string(old_block)])
        pop!(blocks_brain, stack_index, prefix_source)
        put!(blocks_brain, stack_index, block, "T")
        if (verbose)
            println("Saving block ", block, " on the table.")
        end
        if (to_readout)
            read_block_world_instance(blocks_brain,stacks_number,prefix_source)
            read_block_world_instance(blocks_brain,stacks_number,"T")
        end
        node_index, block = top(blocks_brain, stack_index, prefix_source)
    end
    return stack_actions
end


## save the dismantled stack as final stack
# Parameters
# - blocks_brain: type BlocksBrain, after having parsed the stacks
# - stack_index: an integer, it is the index of the stack
# - prefix_source: string "I", "G", "B", or "T" whether the areas we work on are the input areas, the goal areas, the building areas, or the table areas.
# - block_to_keep: Int, the block where to stop (we stop by deleting the one before)
# Example: save_up_to_block!(blocks_brain, 5, "B", 4)
function save_up_to_block!(blocks_brain::BlocksBrain, stack_index::Int, prefix_source::String, block_to_keep::Int, prefix_working::String; verbose::Bool=false)::Nothing
    if (verbose) 
        println("Working on stack number ", stack_index, ".")
    end
    stack_actions = Vector{String}()
    blocks_brain.brain.no_plasticity = true

    # defining locally the regions involved with stack stack_index
    source_head = add_prefix([HEADS[stack_index]], prefix_source)[1]
    source_nodes = add_prefix(NODES[stack_index], prefix_source)
    # working_head = add_prefix([HEADS[stack_index]],prefix_working)[1]
    working_nodes = add_prefix(NODES[stack_index], prefix_working)

    for area in vcat([BLOCKS,source_head], source_nodes, working_nodes)
        unfix_assembly(blocks_brain.brain.areas[area])
    end

    top_index = 1
    source_stack_length, last_index = fire_last_assembly(blocks_brain, source_head, source_nodes, top_index)
    if (verbose)
        println("Last active node in stack ", stack_index, ": ", source_nodes[last_index], ".")
    end

    # preactions
    disinhibit_area(blocks_brain, source_nodes[last_index], 0)

    current_index = last_index
    count = 1
    while true
        if (verbose)
            println("Working on area ", source_nodes[current_index], ".")
        end
        blocks_brain.brain.no_plasticity = true
        next_index = mod(current_index - 1, 1:MAX_NODES_AREAS)
        # preactions
        for area in vcat([BLOCKS], source_nodes, working_nodes)
            unfix_assembly(blocks_brain.brain.areas[area])
        end
        disinhibit_areas(blocks_brain, [source_nodes[next_index],BLOCKS,source_nodes[current_index]], 0)
        disinhibit_fibers(blocks_brain, [[source_nodes[current_index], source_nodes[next_index]],[source_nodes[current_index],BLOCKS]], 0)

        project_map = Dict{String,Set{String}}()
        project_map[source_nodes[current_index]] = Set{String}([source_nodes[current_index],BLOCKS,source_nodes[next_index]])
        project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)

        block = get_block_index(blocks_brain, BLOCKS)

        if count == 1
            # append!(stack_actions,[string("Fix block ", block, " on the table as stack ", stack_index, ".")])
            if (verbose) 
                println(string("Fix block ", block, " on the table as stack ", stack_index, "."))
            end
        else
            append!(stack_actions, [string("Fix block ", block, " from the table to the top of stack ", stack_index, ".")])
            if (verbose) 
                println(string("Fix block ", block, " from the table to the top of stack ", stack_index, "."))
            end
        end

        ## Postactions
        inhibit_areas(blocks_brain, [source_nodes[current_index],source_nodes[next_index]], 0)
        inhibit_fibers(blocks_brain, [[source_nodes[current_index],source_nodes[next_index]],[source_nodes[current_index],BLOCKS]], 0)

        put!(blocks_brain, stack_index, block, prefix_working)

        if (block == block_to_keep)
            break
        end
        current_index = next_index
        count = count + 1
    end
    return 
end

## reassemble blocks for a given goal stack assuming it is already built until starting_block
## WARNING: in this function we loose track of the last active assembly in HEAD for the goal stack. We should project it into a temporary region
# Parameters
# - blocks_brain: type BlocksBrain, after having parsed the stacks
# - stacks_index: an integer, it is the index of the goal stack to reassemble
# - prefix_source: string "I", "G", or "B" whether the areas we work on are the input areas, the goal areas or the building areas.
#                  usually "G", since it is the source to reassemble the blocks
# - prefix_working: string "I", "G", or "B" whether the areas we work on are the input areas, the goal areas or the building areas.
#                  usually "B", since it is the final set of regions where to reassemble the blocks
# - starting_block: Integer, it is the last block already saved
# - input_stacks_number: an integer, it is the number of input stacks
# - input_stacks_number: an integer, it is the number of goal stacks (only for readout pursposes)

# Example: reassemble_stack_from_block!(blocks_brain, 5, "G", "B", 2)
function reassemble_stack_from_block!(blocks_brain::BlocksBrain, stack_index::Int, prefix_source::String, prefix_working::String, starting_block::Int, input_stacks_number::Int, goal_stacks_number::Int; to_readout::Bool = false, verbose::Bool=false)::Vector{String}
    if (verbose) 
        println("Working on stack number ", stack_index, ".")
    end
    stack_actions = Vector{String}()
    blocks_brain.brain.no_plasticity = true

    # defining locally the regions involved with stack stack_index
    source_head = add_prefix([HEADS[stack_index]], prefix_source)[1]
    source_nodes = add_prefix(NODES[stack_index], prefix_source)
    working_head = add_prefix([HEADS[stack_index]], prefix_working)[1]
    working_nodes = add_prefix(NODES[stack_index], prefix_working)

    for area in vcat([BLOCKS,source_head,working_head], source_nodes, working_nodes)
        unfix_assembly(blocks_brain.brain.areas[area])
    end

    top_index = 1
    last_index = activate_node_with_block!(blocks_brain, source_head, source_nodes, top_index, starting_block)
    if (verbose)
        println("Last active node in stack ", stack_index, ": ", source_nodes[last_index], ".")
    end

    # preactions
    disinhibit_area(blocks_brain, source_nodes[last_index], 0)

    current_index = last_index
    next_index = mod(current_index - 1, 1:MAX_NODES_AREAS)
    disinhibit_areas(blocks_brain, [source_nodes[next_index],source_head,source_nodes[current_index]], 0)
    disinhibit_fibers(blocks_brain, [[source_nodes[current_index], source_nodes[next_index]],[source_nodes[current_index],source_head]], 0)
    project_map = Dict{String,Set{String}}()
    project_map[source_nodes[current_index]] = Set{String}([source_nodes[current_index],BLOCKS,source_nodes[next_index],source_head])
    project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)

    stable_assembly = is_assembly(blocks_brain, source_head)
    inhibit_areas(blocks_brain, [source_nodes[current_index],source_nodes[next_index]], 0)
    inhibit_fibers(blocks_brain, [[source_nodes[current_index],source_nodes[next_index]],[source_head,source_nodes[current_index]]], 0)


    current_index = next_index
    count = 1
    while !(stable_assembly)
        if (verbose)
            println("Working on area ", source_nodes[current_index], ".")
        end
        blocks_brain.brain.no_plasticity = true
        next_index = mod(current_index - 1, 1:MAX_NODES_AREAS)
        # preactions
        for area in vcat([BLOCKS,source_head,working_head], source_nodes, working_nodes)
            unfix_assembly(blocks_brain.brain.areas[area])
        end
        disinhibit_areas(blocks_brain, [source_nodes[next_index],BLOCKS,source_head,source_nodes[current_index]], 0)
        disinhibit_fibers(blocks_brain, [[source_nodes[current_index], source_nodes[next_index]],[source_nodes[current_index],BLOCKS],[source_nodes[current_index],source_head]], 0)

        project_map = Dict{String,Set{String}}()
        project_map[source_nodes[current_index]] = Set{String}([source_nodes[current_index],BLOCKS,source_nodes[next_index],source_head])
        project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)

        stable_assembly = is_assembly(blocks_brain, source_head)
        block = get_block_index(blocks_brain, BLOCKS)

        # if count == 1
        #     append!(stack_actions,[string("Fix block ", block, " on the table as stack ", stack_index, ".")])
        #     if (verbose) 
        #         println(string("Fix block ", block, " on the table as stack ", stack_index, "."))
        #     end
        # else
        #     append!(stack_actions,[string("Move block ", block, " from the table to the top of stack ", stack_index, ".")])
        #     if (verbose) 
        #         println(string("Move block ", block, " from the table to the top of stack ", stack_index, "."))
        #     end
        # end
        # append!(stack_actions,[string("Move block ", block, " from the table to the top of stack ", stack_index, ".")])
        append!(stack_actions, [string(block)])
        if (verbose) 
            println(string("Move block ", block, " from the table to the top of stack ", stack_index, "."))
        end

        ## Postactions
        inhibit_areas(blocks_brain, [source_nodes[current_index],source_nodes[next_index]], 0)
        inhibit_fibers(blocks_brain, [[source_nodes[current_index],source_nodes[next_index]],[source_nodes[current_index],BLOCKS],[source_head,source_nodes[current_index]]], 0)

        put!(blocks_brain, stack_index, block, prefix_working)

        set_as_relocated!(blocks_brain, input_stacks_number, block, "T")

        if (verbose)
            println("Block ", block, " has been recollocated.")
            println("There's a stable assembly in ", source_head, ": ", stable_assembly, ".")
        end

        if (to_readout)
            read_block_world_instance(blocks_brain,goal_stacks_number,prefix_working)
            read_block_world_instance(blocks_brain,input_stacks_number,"T")
        end

        current_index = next_index
        count = count + 1
    end
    return stack_actions
end


## activates the assembly in NODE corresponding to the given block in the stack
# Parameters
# - blocks_brain: type BlocksBrain, after having parsed the stacks
# - head: a string representing the HEAD area of the interested stack
# - nodes: a vector of String objects representing the NODE areas of the interested stack
# - top_index: an integer representing the index of the NODE area where the first block of the stack is saved
# - block: an integer - the block we have to find
# Example: activate_node_with_block!(blocks_brain, "T2_H", ["T2_N1","T2_N2","T2_N3"], 2,3)
function activate_node_with_block!(blocks_brain::BlocksBrain, head::String, nodes::Vector{String}, top_index::Int, block::Int; verbose::Bool=false)::Int
    
    blocks_brain.brain.no_plasticity = true
    for area in vcat([head], nodes)
        unfix_assembly(blocks_brain.brain.areas[area])
    end

    # preactions
    disinhibit_areas(blocks_brain, [head,nodes[top_index],BLOCKS], 0)
    disinhibit_fiber(blocks_brain, head, nodes[top_index], 0)

    project_map = Dict{String,Set{String}}()
    project_map[head] = Set{String}([head,nodes[top_index]])
    project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)

    # postactions
    inhibit_area(blocks_brain, head, 0)
    inhibit_fiber(blocks_brain, head, nodes[top_index], 0)

    count = 1
    current_index = top_index
    last_index = 0

    block_found = false

    while !(block_found)
        next_index = mod(current_index + 1, 1:MAX_NODES_AREAS)
        # preactions
        disinhibit_area(blocks_brain, nodes[next_index], 0)
        disinhibit_fibers(blocks_brain, [[nodes[current_index],nodes[next_index]],[nodes[current_index],BLOCKS]], 0)

        project_map = Dict{String,Set{String}}()
        project_map[nodes[current_index]] = Set{String}([nodes[current_index],nodes[next_index],BLOCKS])
        project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)

        last_block = get_block_index(blocks_brain, BLOCKS)
        if last_block == block
            block_found = true
        end

        # postactions
        inhibit_area(blocks_brain, nodes[current_index], 0)
        inhibit_fibers(blocks_brain, [[nodes[current_index],nodes[next_index]],[nodes[current_index],BLOCKS]], 0)

        last_index = current_index
        current_index = next_index
        count = count + 1     
    end

    inhibit_areas(blocks_brain, [nodes[current_index],BLOCKS], 0)

    current_index = last_index
    if current_index == 0
        error("No blcok ", block, " in the goal stack.")
    end

    if (verbose)
        println("Block ", block, " in the goal stack is saved in ", nodes[current_index])
    end

    blocks_brain.brain.no_plasticity = false
    return current_index
end