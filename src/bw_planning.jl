# this region saves which blocks have been moved from the TABLE to some stack
global RELOCATED = "RELOCATED"

## fires the assembly in NODE corresponding to the last block in the stack: return the length of the stack and the node index of the last block
## WARNING: the return length is only for readout purposes. The AC should not be able to do this with the current model
## for the TABLE areas, it DOES COUNT also the relocated blocks
# Parameters
# - blocks_brain: type BlocksBrain, after having parsed the stacks
# - head: a string representing the HEAD area of the interested stack
# - nodes: a vector of String objects representing the NODE areas of the interested stack
# - top_index: an integer representing the index of the NODE area where the first block of the stack is saved
# Example: fire_last_assembly(blocks_brain, "T2_H", ["T2_N1","T2_N2","T2_N3"], 2)

function fire_last_assembly(blocks_brain::BlocksBrain, head::String, nodes::Vector{String}, top_index::Int; verbose::Bool = false)::Tuple{Int,Int}
    
    blocks_brain.brain.no_plasticity = true
    for area in vcat([head],nodes)
        unfix_assembly(blocks_brain.brain.areas[area])
    end

    # preactions
    disinhibit_areas(blocks_brain,[head,nodes[top_index]],0)
    disinhibit_fiber(blocks_brain,head,nodes[top_index],0)

    project_map = Dict{String,Set{String}}()
    project_map[head] = Set{String}([head,nodes[top_index]])
    project(blocks_brain.brain,Dict{String,Set{String}}(),project_map)

    # postactions
    inhibit_area(blocks_brain,head,0)
    inhibit_fiber(blocks_brain,head,nodes[top_index],0)

    count = 1
    current_index = top_index
    last_index = 0

    while (is_assembly(blocks_brain,nodes[current_index]))
        next_index = mod(current_index + 1, 1:MAX_NODES_AREAS)
        #preactions
        disinhibit_area(blocks_brain,nodes[next_index],0)
        disinhibit_fibers(blocks_brain,[[nodes[current_index],nodes[next_index]]],0)

        project_map = Dict{String,Set{String}}()
        project_map[nodes[current_index]] = Set{String}([nodes[current_index],nodes[next_index]])
        project(blocks_brain.brain,Dict{String,Set{String}}(),project_map)

        inhibit_area(blocks_brain,nodes[current_index],0)
        inhibit_fiber(blocks_brain,nodes[current_index],nodes[next_index],0)

        last_index = current_index
        current_index = next_index
        count = count + 1     
    end

    inhibit_area(blocks_brain,nodes[current_index],0)

    current_index = last_index
    if current_index == 0
        error("No blocks in the stack. Top index was ", top_index, ".")
    end

    if (verbose)
        println("Last block in the stack is saved in ", nodes[current_index])
    end

    blocks_brain.brain.no_plasticity = false
    return (count-1,current_index)
end


## reads the block world instance and prints it to standard output, in particular TABLE and the stacks with the given prefix
## Parameters
# - blocks_brain:: BlocksBrain
# - stacks_number:: how many stacks there are
# - prefix:: a string, the preifx of the stacks

function read_block_world_instance(blocks_brain::BlocksBrain,stacks_number::Int, prefix::String)::Nothing
    top_areas = Vector{Int}(undef,stacks_number)
    stack_lengths = Vector{Int}(undef,stacks_number)
    for j in 1:stacks_number
        head = add_prefix([HEADS[j]],prefix)[1]
        nodes = add_prefix(NODES[j],prefix)
        top_areas[j],block = top(blocks_brain,j,prefix)
        if top_areas[j] == 0
            stack_lengths[j] = 0
        else
            stack_lengths[j], last_index = fire_last_assembly(blocks_brain,head,nodes,top_areas[j]) 
        end
    end
    input_stacks = readout(blocks_brain,stacks_number,stack_lengths,top_areas,prefix)
    kind_of_stack = ""
    if prefix == "I"
        kind_of_stack = "Input"
    elseif prefix == "B"
        kind_of_stack = "Final"
    elseif prefix == "G"
        kind_of_stack = "Goal"
    elseif prefix == "T"
        kind_of_stack = "Table"
    end

    if prefix == "T"
        all_blocks = Vector{Int}()
        for j in 1:stacks_number
            append!(all_blocks, input_stacks[j])
        end
        println("Blocks on the table: ", all_blocks, ".")
    else 
        for j in 1:stacks_number
            println(kind_of_stack," stack ", j,": ", input_stacks[j])
        end
    end
    flush(stdout)
    # head = add_prefix([HEADS[1]],"T")[1]
    # nodes = add_prefix(NODES[1],"T")
    # top_area_table,block = top(blocks_brain,1,"T")
    # if top_area_table == 0
    #     stack_length = 0
    # else
    #     stack_length, last_index = fire_last_assembly(blocks_brain,head,nodes,top_area_table)
    # end
    #println("Blocks on the table: ", readout(blocks_brain,1,[stack_length],[top_area_table],"T")[1])
    return
end

## dismantle stacks and put each block on the table
# Parameters
# - blocks_brain: type BlocksBrain, after having parsed the stacks
# - stacks_number: an integer, it is the number of stacks
# - prefix: string "I", "G", "B", or "T" whether the areas we work on are the input areas, the goal areas, the building areas, or the table areas.
# Example: dismantle!(blocks_brain, 5, "I")
function dismantle!(blocks_brain::BlocksBrain, stacks_number::Int, prefix_source::String; to_readout::Bool = false, verbose::Bool = false)::Vector{Vector{String}}
    myactions = Vector{Vector{String}}(undef,stacks_number)
    # we dismantle each stack
    for j in 1:stacks_number
        myactions[j] = dismantle_stack!(blocks_brain,j,prefix_source,stacks_number; to_readout, verbose)
        # if (to_readout)
        #     read_block_world_instance(blocks_brain,stacks_number,prefix_source)
        #     read_block_world_instance(blocks_brain,stacks_number,"T")
        # end
    end
    return myactions
end

## dismantle stack and put each block on the table
# Parameters
# - blocks_brain: type BlocksBrain, after having parsed the stacks
# - stacks_index: an integer, it is the index of the stack
# - prefix: string "I", "G", "B", or "T" whether the areas we work on are the input areas, the goal areas, the building areas, or the table areas.
# - stacks_number: integer, the number of stacks (only for readout purposes)
# Example: dismantle_stack!(blocks_brain, 5, "I",3)
function dismantle_stack!(blocks_brain::BlocksBrain, stack_index::Int, prefix_source::String, stacks_number::Int; to_readout::Bool = false, verbose::Bool = false)::Vector{String}
    stack_actions = Vector{String}()
    node_index, block = top(blocks_brain,stack_index,prefix_source; verbose = false)

    while !(is_last_block(blocks_brain,stack_index,node_index,block,prefix_source))
        verbose && println("Now popping block ", block, " in area", NODES[stack_index][node_index], ".")
        old_block = block
        #append!(stack_actions, [string("Remove block ", old_block, " from top of stack ", stack_index, ".")])
        append!(stack_actions, [string(old_block)])
        pop!(blocks_brain,stack_index,prefix_source; verbose = false)

        verbose && println("Saving block ", block, " on the table.")
        put!(blocks_brain,stack_index,block,"T"; verbose = false)
        if (to_readout)
            read_block_world_instance(blocks_brain,stacks_number,prefix_source)
            read_block_world_instance(blocks_brain,stacks_number,"T")
        end
        node_index, block = top(blocks_brain,stack_index,prefix_source; verbose= false)
        verbose && println("New top area ", NODES[stack_index][node_index])
        verbose && println("New top block ", block)

    end
    verbose && println("Block ", block, " is last block.")
    pop!(blocks_brain,stack_index,prefix_source; verbose = false)
    verbose && println("Saving block ", block, " on the table.")
    put!(blocks_brain,stack_index,block,"T"; verbose= false )
    if (to_readout)
        read_block_world_instance(blocks_brain,stacks_number,prefix_source)
        read_block_world_instance(blocks_brain,stacks_number,"T")
    end
    return stack_actions
end

## reassemble blocks 
## WARNING: in this function we loose track of the last active assembly in HEAD for the goal stack. We should project it into a temporary region
# Parameters
# - blocks_brain: type BlocksBrain, after having parsed the stacks
# - input_stacks_number: an integer, it is the number of input stacks
# - goal_stacks_number: an integer, it is the number of goal stacks
# - prefix_source: string "I", "G", or "B" whether the areas we work on are the input areas, the goal areas or the building areas.
#                  usually "G", since it is the source to reassemble the blocks
# - prefix_working: string "I", "G", or "B" whether the areas we work on are the input areas, the goal areas or the building areas.
#                  usually "B", since it is the final set of regions where to reassemble the blocks
# Example: reassemble!(blocks_brain, 5, "G", "B")
function reassemble!(blocks_brain::BlocksBrain,input_stacks_number::Int,goal_stacks_number::Int,prefix_source::String,prefix_working::String; to_readout::Bool = false, verbose::Bool = false)::Vector{Vector{String}}
    
    myactions = Vector{Vector{String}}(undef,goal_stacks_number)
    # adding the prefix to the regions
    source_regions = add_prefix(vcat(REGIONS...),prefix_source)
    working_regions = add_prefix(vcat(REGIONS...),prefix_working)
    for j in 1:goal_stacks_number
        myactions[j] = reassemble_stack!(blocks_brain,j,prefix_source,prefix_working,input_stacks_number,goal_stacks_number;to_readout,verbose)
    end
    blocks_brain.brain.no_plasticity = true
    inhibit_areas(blocks_brain, vcat([BLOCKS],source_regions,working_regions),0)
    inhibit_all_fibers(blocks_brain, vcat([BLOCKS],source_regions,working_regions),0)
    return myactions
end

## reassemble blocks for a given goal stack
## WARNING: in this function we loose track of the last active assembly in HEAD for the goal stack. We should project it into a temporary region
# Parameters
# - blocks_brain: type BlocksBrain, after having parsed the stacks
# - stacks_index: an integer, it is the index of the goal stack to reassemble
# - prefix_source: string "I", "G", or "B" whether the areas we work on are the input areas, the goal areas or the building areas.
#                  usually "G", since it is the source to reassemble the blocks
# - prefix_working: string "I", "G", or "B" whether the areas we work on are the input areas, the goal areas or the building areas.
#                  usually "B", since it is the final set of regions where to reassemble the blocks
# - input_stacks_number: integer, the number of input stacks
# - goal_stacks_number: integer, the number of goal stacks (only for readout purposes)
# Example: reassemble_stack!(blocks_brain, 5, "G", "B")
function reassemble_stack!(blocks_brain::BlocksBrain,stack_index::Int,prefix_source::String,prefix_working::String,input_stacks_number::Int,goal_stacks_number::Int; to_readout::Bool = false, verbose::Bool = false)::Vector{String}
    if (verbose) 
        println("Working on stack number ", stack_index, ".")
    end
    stack_actions = Vector{String}()
    blocks_brain.brain.no_plasticity = true

    # defining locally the regions involved with stack stack_index
    source_head = add_prefix([HEADS[stack_index]],prefix_source)[1]
    source_nodes = add_prefix(NODES[stack_index],prefix_source)
    working_head = add_prefix([HEADS[stack_index]],prefix_working)[1]
    working_nodes = add_prefix(NODES[stack_index],prefix_working)

    for area in vcat([BLOCKS,source_head,working_head],source_nodes,working_nodes)
        unfix_assembly(blocks_brain.brain.areas[area])
    end

    top_index = 1
    stack_length,last_index = fire_last_assembly(blocks_brain,source_head,source_nodes,top_index)
    if (verbose)
        println("Last active node in stack ", stack_index, ": ", source_nodes[last_index],".")
    end

    # preactions
    disinhibit_area(blocks_brain,source_nodes[last_index],0)

    current_index = last_index
    count = 1
    while true
        if (verbose)
            println("Working on area ", source_nodes[current_index], ".")
        end
        blocks_brain.brain.no_plasticity = true
        next_index = mod(current_index - 1, 1:MAX_NODES_AREAS)
        # preactions
        for area in vcat([BLOCKS,source_head],source_nodes)
            unfix_assembly(blocks_brain.brain.areas[area])
        end
        disinhibit_areas(blocks_brain,[source_nodes[next_index],BLOCKS,source_head,source_nodes[current_index]],0)
        disinhibit_fibers(blocks_brain, [[source_nodes[current_index], source_nodes[next_index]],[source_nodes[current_index],BLOCKS],[source_nodes[current_index],source_head]],0)

        project_map = Dict{String,Set{String}}()
        project_map[source_nodes[current_index]] = Set{String}([source_nodes[current_index],BLOCKS,source_nodes[next_index],source_head])
        project(blocks_brain.brain,Dict{String,Set{String}}(),project_map)

        stable_assembly = is_assembly(blocks_brain,source_head)
        block = get_block_index(blocks_brain,BLOCKS)

        if count == 1
            #append!(stack_actions,[string("Fix block ", block, " on the table as stack ", stack_index, ".")])
            append!(stack_actions,[string(block)])
            if (verbose) 
                println(string("Fix block ", block, " on the table as stack ", stack_index, "."))
            end
        else
            #append!(stack_actions,[string("Move block ", block, " from the table to the top of stack ", stack_index, ".")])
            append!(stack_actions,[string(block)])
            if (verbose) 
                println(string("Move block ", block, " from the table to the top of stack ", stack_index, "."))
            end
        end

        ## Postactions
        inhibit_areas(blocks_brain,[source_nodes[current_index],source_nodes[next_index]],0)
        inhibit_fibers(blocks_brain,[[source_nodes[current_index],source_nodes[next_index]],[source_nodes[current_index],BLOCKS],[source_head,source_nodes[current_index]]],0)

        put!(blocks_brain,stack_index,block,prefix_working)
        set_as_relocated!(blocks_brain,input_stacks_number,block,"T"; verbose = false)

        if (verbose)
            println("Block ", block, " has been recollocated.")
            println("There's a stable assembly in ", source_head, ": ", stable_assembly, ".")
        end

        if (to_readout)
            read_block_world_instance(blocks_brain,goal_stacks_number,prefix_working)
            read_block_world_instance(blocks_brain,input_stacks_number,"T")
        end

        if (stable_assembly)
            break
        end
        current_index = next_index
        count = count + 1
    end
    return stack_actions
end

function set_as_relocated!(blocks_brain::BlocksBrain,stacks_number::Int,block::Int,prefix::String; verbose::Bool = false)::Nothing
    blocks_brain.brain.no_plasticity = true
    # adding prefix to regions
    # defining locally the regions involved

    for stack_index in 1:stacks_number
        nodes = add_prefix(NODES[stack_index],prefix)

        disinhibit_areas(blocks_brain,vcat([BLOCKS,RELOCATED],nodes),0)
        for area in vcat([BLOCKS,RELOCATED],nodes)
            unfix_assembly(blocks_brain.brain.areas[area])
        end
        activate_block(blocks_brain,block)
        for area in nodes
            disinhibit_fibers(blocks_brain,[[BLOCKS,area],[area,RELOCATED]],0)
            project_map = Dict(BLOCKS => Set([BLOCKS,area]))
            project(blocks_brain.brain, Dict{String,Set{String}}(),project_map)

            if (is_assembly(blocks_brain,area))
                blocks_brain.brain.no_plasticity = false

                inhibit_areas(blocks_brain,vcat([BLOCKS,RELOCATED],nodes), 0)
                inhibit_all_fibers(blocks_brain,vcat([BLOCKS,RELOCATED],nodes),0)

                project_map = Dict(area => Set([RELOCATED,area]))
                project(blocks_brain.brain, Dict{String,Set{String}}(),project_map)
                disinhibit_areas(blocks_brain,[area,RELOCATED],0)
                disinhibit_fiber(blocks_brain,area,RELOCATED,0)
                proj_rounds = 50
                project_star(blocks_brain; project_rounds = proj_rounds)
                if (verbose)
                    println("Block ", block, " has been relocated:", is_assembly(blocks_brain,RELOCATED))
                end
                inhibit_all_fibers(blocks_brain,vcat([BLOCKS,RELOCATED],nodes),0)
                inhibit_areas(blocks_brain,vcat([BLOCKS,RELOCATED],nodes),0)
                for area in vcat([BLOCKS,RELOCATED],nodes)
                    unfix_assembly(blocks_brain.brain.areas[area])
                end
                return
            end
            inhibit_fibers(blocks_brain,[[BLOCKS,area],[area,RELOCATED]],0)

        end
        inhibit_all_fibers(blocks_brain,vcat([BLOCKS,RELOCATED],nodes),0)
        inhibit_areas(blocks_brain,vcat([BLOCKS,RELOCATED],nodes),0)
    end
    error("Block ", block, " is NOT on the TABLE.")
    return
end