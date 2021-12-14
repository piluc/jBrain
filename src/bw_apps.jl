import Base.pop!
import Base.put!

# The number of node areas to be used to represent a block world instance. The default number is three.
global MAX_NODES_AREAS = 3

# The maximum number of stacks in the block rowld instance.
global MAX_STACKS = 5

# The "input" node areas. They are of the form " I2_N3 ", which means "input stack number 2, node area number 3"
# or " G1_N2 ", which means "goal stack number 2, node area number 3". The prefix will be given later and can be either "I" ("input"), "G" ("goal"), "B" ("building").
node_areas = Vector{Vector{String}}()
for j in 1:MAX_STACKS
    node_areas_stack_j = Vector{String}()
    for k in 1:MAX_NODES_AREAS
        append!(node_areas_stack_j, [string(j, "_N", k)])
    end
    append!(node_areas, [node_areas_stack_j])
end

global NODES = node_areas

# H stands for "heads" and stores the top block of each stack. The prefix will be given later and can be either "I" ("input"), "G" ("goal"), "B" ("building").
heads = Vector{String}()
for j in 1:MAX_STACKS
    append!(heads, [string(j, "_H")])
end

global HEADS = heads

# REGIONS collects all regions. It is a vector of vectors of strings, one for each stack. The prefix will be given later and can be either "I" ("input"), "G" ("goal"), "B" ("building").
regions = Vector{Vector{String}}()
for j in 1:MAX_STACKS
    regions_stack_j = vcat(NODES[j], [HEADS[j]])
    append!(regions, [regions_stack_j])
end
global REGIONS = regions

# some preliminary functions

function disinhibit_areas(blocks_brain::BlocksBrain, area_names, lock)
    for area_name in area_names
        disinhibit_area(blocks_brain, area_name, lock)
    end
end

function disinhibit_fibers(blocks_brain::BlocksBrain, area_pairs, lock)
    for pair in area_pairs
        disinhibit_fiber(blocks_brain, pair[1], pair[2], lock)
    end
end

function disinhibit_all_fibers(blocks_brain::BlocksBrain, area_names, lock)
    for area1 in area_names
        for area2 in area_names
            if area1 != area2
                disinhibit_fiber(blocks_brain, area1, area2, lock)
            end
        end
    end
end

function inhibit_areas(blocks_brain::BlocksBrain, area_names, lock)
    for area_name in area_names
        inhibit_area(blocks_brain, area_name, lock)
    end
end

function inhibit_fibers(blocks_brain::BlocksBrain, area_pairs, lock)
    for pair in area_pairs
        inhibit_fiber(blocks_brain, pair[1], pair[2], lock)
    end
end

function inhibit_all_fibers(blocks_brain::BlocksBrain, area_names, lock)
    for area1 in area_names
        for area2 in area_names
            if area1 != area2
                inhibit_fiber(blocks_brain, area1, area2, lock)
            end
        end
    end
end


function add_prefix(regions::Vector{String}, prefix::String)::Vector{String}
    new_regions = Vector{String}()
    for area in regions
        append!(new_regions, [string(prefix, area)])
    end
    return new_regions
end


# Parameters
# - blocks: array of blocks, that is, permutation of number between 1 and number of blocks
# - query_a: first block of the query (has to be one of the blocks in the stack)
# - query_b: second block of the query (has to be one of the blocks in the stack)
# - p: erdos-renyi parameter
# - eak: k for explicit areas
# - nean: n for non explicit areas
# - neak: k for non explicit areas
# - db: defaul plasticity
# Example: is_above([1,2,3,4],1,4,0.1,10,10000,100,0.2)
# This is the hello world program in the paper
function is_above(blocks::Vector{Int}, query_a::Int, query_b::Int, p::Float64, eak::Int, nean::Int, neak::Int, db::Float64)::Bool
    blocks_number = length(blocks)
    blocks_brain = BlocksBrain(blocks_number, [FIRST, SECOND, ABOVE, RELATION], p, eak, nean, neak, db)

    disinhibit_area(blocks_brain, BLOCKS, 0)
    disinhibit_area(blocks_brain, FIRST, 0)
    disinhibit_fiber(blocks_brain, BLOCKS, FIRST, 0)

    activate_block(blocks_brain, query_a)
    project_star(blocks_brain)

    inhibit_area(blocks_brain, FIRST, 0)
    disinhibit_area(blocks_brain, SECOND, 0)
    disinhibit_fiber(blocks_brain, BLOCKS, SECOND, 0)

    activate_block(blocks_brain, query_b)
    project_star(blocks_brain)

    for block in blocks
        inhibit_area(blocks_brain, SECOND, 0)
        disinhibit_area(blocks_brain, FIRST, 0)
        activate_block(blocks_brain, block)
        project_map = Dict{String,Set{String}}()
        project_map[BLOCKS] = Set{String}([BLOCKS,FIRST])        
        project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)
        if (is_assembly(blocks_brain, FIRST))
            return true
        end
        inhibit_area(blocks_brain, FIRST, 0)
        disinhibit_area(blocks_brain, SECOND, 0)
        activate_block(blocks_brain, block)
        project_map = Dict{String,Set{String}}()
        project_map[BLOCKS] = Set{String}([BLOCKS,SECOND])        
        project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)
        if (is_assembly(blocks_brain, SECOND))
            return false
        end
    end
end

# Creates and returns a brain

# Parameters
# - blocks_number: number of blocks involved
# - working_regions: Vector of area names involved
# - p: erdos-renyi parameter
# - eak: k for explicit areas
# - nean: n for non explicit areas
# - neak: k for non explicit areas
# - db: defaul plasticity
# Example: brain_creation(4,[FIRST,SECOND,ABOVE],0.1,50,1000000,50,0.2)

function brain_creation(blocks_number::Int, working_regions::Vector{String}, p::Float64, eak::Int, nean::Int, neak::Int, db::Float64)::BlocksBrain
    return BlocksBrain(blocks_number, working_regions, p, eak, nean, neak, db)
end


# Parser 

# Parameters
# - blocks_brain:: a brain (type BlocksBrain) 
#                  created with a number of blocks which is at least as big as the number of blocks we want to parse
#                  and with the regions with the given prefix
# - stacks: array (max length MAX_STACKS) of arrays of blocks, the lasts being integer from 1 to number of blocks
# - p: erdos-renyi parameter
# - eak: k for explicit areas
# - nean: n for non explicit areas
# - neak: k for non explicit areas
# - db: defaul plasticity
# - prefix: string "I", "G", "B", or "T" whether the areas we work on are the input areas, the goal areas, the building areas, or the table areas.
# Example: parse!(blocks_brain,[[1,2,3,4],[5,6]],"I")
function parse!(blocks_brain::BlocksBrain, stacks::Vector{Vector{Int}}, prefix::String; verbose::Bool=false, project_rounds=50)::Nothing
    stacks_number = length(stacks)
    blocks_number = 0
    for j in 1:stacks_number
        blocks_number = blocks_number + length(stacks[j])
    end

    # adding the prefix to the regions
    working_regions = add_prefix(vcat(REGIONS...), prefix)

    # we parse each stack
    for j in 1:stacks_number
        blocks = stacks[j]

        # defining locally the regions involved to parse stack number j
        head = add_prefix([HEADS[j]], prefix)[1]
        nodes = add_prefix(NODES[j], prefix)

        # Global preactions
        disinhibit_area(blocks_brain, BLOCKS, 0)
        disinhibit_area(blocks_brain, head, 0)
        disinhibit_area(blocks_brain, nodes[1], 0)
        disinhibit_fiber(blocks_brain, head, nodes[1], 0)
        disinhibit_fiber(blocks_brain, BLOCKS, nodes[1], 0)
        count = 1

        for block in blocks
            previous_to_area_index = ((count - 2) % MAX_NODES_AREAS) + 1
            to_area_index = ((count - 1) % MAX_NODES_AREAS) + 1
            next_to_area_index = (count % MAX_NODES_AREAS) + 1
            activate_block(blocks_brain, block)

            if count == 1
                # there is no difference between the "if" statement and the "else" statement
                project_star(blocks_brain; project_rounds = 100, verbose=false)
            else
                project_star(blocks_brain; project_rounds, verbose=false)
            end

            for area in working_regions
                unfix_assembly(blocks_brain.brain.areas[area])
            end

            # Postactions
            if count == 1
                inhibit_area(blocks_brain, head, 0)
                inhibit_fiber(blocks_brain, head, nodes[1], 0)
                inhibit_fiber(blocks_brain, BLOCKS, nodes[1], 0)
            else
                inhibit_area(blocks_brain, nodes[previous_to_area_index], 0)
                inhibit_fiber(blocks_brain, BLOCKS, nodes[to_area_index], 0)
                inhibit_fiber(blocks_brain, nodes[previous_to_area_index], nodes[to_area_index], 0)
            end
            # Preactions
            disinhibit_area(blocks_brain, nodes[next_to_area_index], 0)
            disinhibit_fiber(blocks_brain, BLOCKS, nodes[next_to_area_index], 0)
            disinhibit_fiber(blocks_brain, nodes[to_area_index], nodes[next_to_area_index], 0)
            count = count + 1
        end
        inhibit_areas(blocks_brain, vcat([BLOCKS], working_regions), 0)
        inhibit_all_fibers(blocks_brain, vcat([BLOCKS], working_regions), 0)
    end
    return
end

# Parameters
# - blocks_brain: a BlocksBrain represting a stack of blocks
# - stacks_number: number of stacks in the block world
# - stacks_lengths: array of numbers of blocks in each stack
# - top_areas: array of indeces of the working areas corresponding to the top block for each stack
# - prefix: string "I", "G", "B", or "T" whether the areas we work on are the input areas, the goal areas, the building areas, or the table areas.
# Example: readout(bb,4,[2,3,1,3],[1,1],"G")
function readout(blocks_brain::BlocksBrain, stacks_number::Int, stacks_lengths::Vector{Int}, top_areas::Vector{Int}, prefix::String; verbose::Bool=false)::Vector{Vector{Int}}
    stacks = Vector{Vector{Int}}()
    blocks_brain.brain.no_plasticity = true

    # adding prefix to regions
    working_regions = add_prefix(vcat(REGIONS...), prefix)
    if prefix == "T"
        append!(working_regions, [RELOCATED])
    end
        

    # unfix all assemblies for readout
	for area in vcat([BLOCKS], working_regions)
        unfix_assembly(blocks_brain.brain.areas[area])
    end

    for j in 1:stacks_number
        
        # the blocks of stack j
        blocks = Vector{Int}()

        if top_areas[j] == 0
            append!(stacks, [blocks])
            continue
        end

        # defining locally the regions involved to parse stack number j
        head = add_prefix([HEADS[j]], prefix)[1]
        nodes = add_prefix(NODES[j], prefix)

        # Preactions
        disinhibit_area(blocks_brain, head, 0)
        disinhibit_area(blocks_brain, nodes[top_areas[j]], 0)
        disinhibit_fiber(blocks_brain, head, nodes[top_areas[j]], 0)
        project(blocks_brain.brain, Dict{String,Set{String}}(), Dict{String,Set{String}}(head => Set([head,nodes[top_areas[j]]])))
        # Postactions
        inhibit_area(blocks_brain, head, 0)
        inhibit_fiber(blocks_brain, head, nodes[top_areas[j]], 0)
        # Preactions
        disinhibit_area(blocks_brain, BLOCKS, 0)
        disinhibit_fiber(blocks_brain, nodes[top_areas[j]], BLOCKS, 0)

        project_map = Dict{String,Set{String}}()
        if prefix == "T" # if the block was on the table, we must check if it has been RELOCATED
            disinhibit_area(blocks_brain, RELOCATED, 0)
            disinhibit_fiber(blocks_brain, nodes[top_areas[j]], RELOCATED, 0)
            project_map[nodes[top_areas[j]]] = Set{String}([BLOCKS,nodes[top_areas[j]],RELOCATED])
        else
            project_map[nodes[top_areas[j]]] = Set{String}([BLOCKS,nodes[top_areas[j]]])
        end

        project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)

        if !(is_assembly(blocks_brain, RELOCATED)) 
        # if the prefix is not TABLE, RELOCATED is inhibited and no assembly will be found
        # if the prefix is TABLE, it only appends the block if the block has not been relocated
            append!(blocks, get_block_index(blocks_brain, BLOCKS, 0.75))
        end
        
        # Postactions
        inhibit_area(blocks_brain, BLOCKS, 0)
        inhibit_fiber(blocks_brain, nodes[top_areas[j]], BLOCKS, 0)
        current_area = top_areas[j]
        next_area = mod(top_areas[j] + 1, 1:MAX_NODES_AREAS)
        for _ in 2:stacks_lengths[j]
            # Preactions
            disinhibit_area(blocks_brain, nodes[next_area], 0)
            disinhibit_fiber(blocks_brain, nodes[current_area], nodes[next_area], 0)
            project(blocks_brain.brain, Dict{String,Set{String}}(), Dict{String,Set{String}}(nodes[current_area] => Set([nodes[next_area]])))
            # Postactions
            inhibit_area(blocks_brain, nodes[current_area], 0)
            # Preactions
            disinhibit_area(blocks_brain, BLOCKS, 0)
            disinhibit_fiber(blocks_brain, nodes[next_area], BLOCKS, 0)
            project_map = Dict{String,Set{String}}()
            if prefix == "T" # if the block was on the table, we must check if it has been RELOCATED
                disinhibit_area(blocks_brain, RELOCATED, 0)
                disinhibit_fiber(blocks_brain, nodes[next_area], RELOCATED, 0)
                project_map[nodes[next_area]] = Set{String}([BLOCKS,nodes[next_area],RELOCATED])
            else
                project_map[nodes[next_area]] = Set{String}([BLOCKS,nodes[next_area]])
            end

            project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)

            if !(is_assembly(blocks_brain, RELOCATED)) 
            # if the prefix is not TABLE, RELOCATED is inhibited and no assembly will be found
            # if the prefix is TABLE, it only appends the block if the block has not been relocated
                append!(blocks, get_block_index(blocks_brain, BLOCKS, 0.75))
            end
            # Postactions
            inhibit_area(blocks_brain, RELOCATED, 0)
            inhibit_fiber(blocks_brain, nodes[next_area], RELOCATED, 0)
            inhibit_area(blocks_brain, BLOCKS, 0)
            inhibit_fiber(blocks_brain, nodes[next_area], BLOCKS, 0)
            current_area = next_area
            next_area = mod(next_area + 1, 1:MAX_NODES_AREAS)
        end
        append!(stacks, [blocks])
        inhibit_areas(blocks_brain, vcat([BLOCKS], working_regions), 0)
        inhibit_all_fibers(blocks_brain, vcat([BLOCKS], working_regions), 0)
    end
    blocks_brain.brain.no_plasticity = false
    return stacks
end

## Find the first areas in which an assembly is activated by firing from
# HEADS, and return corresponding node area index and block index (for each stack) in a dictionary
# Parameters
# - blocks_brain: a BlocksBrain representing stacks of blocks
# - stack_index: the index of the stack we look at
# - prefix: string "I", "G", "B", or "T" whether the areas we work on are the input areas, the goal areas, the building areas, or the table areas.
# Example: top(bb,4,"G")
function top(blocks_brain::BlocksBrain, stack_index::Int, prefix::String; verbose::Bool=false)::Tuple{Int,Int}
    blocks_brain.brain.no_plasticity = true

    # adding prefix to regions
    working_regions = add_prefix(vcat(REGIONS...), prefix)

    for area in vcat([BLOCKS], working_regions)
        unfix_assembly(blocks_brain.brain.areas[area])
    end

    # defining locally the regions involved to parse stack number stack_index
    head = add_prefix([HEADS[stack_index]], prefix)[1]
    nodes = add_prefix(NODES[stack_index], prefix)
    
    # Global preactions
    disinhibit_area(blocks_brain, head, 0)
    if (is_assembly(blocks_brain, head))
        for index in 1:MAX_NODES_AREAS
            # Preactions
            disinhibit_area(blocks_brain, nodes[index], 0)
            disinhibit_fiber(blocks_brain, head, nodes[index], 0)
            project_map = Dict(head => Set([head, nodes[index]]))
            project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)
            
            if (is_assembly(blocks_brain, nodes[index], 0.75))
                verbose && println("Area ", nodes[index], " is top area.")
                block = 0
                # Postactions
                inhibit_area(blocks_brain, head, 0)
                inhibit_fiber(blocks_brain, head, nodes[index], 0)

                disinhibit_areas(blocks_brain, [nodes[index],BLOCKS], 0)
                disinhibit_fiber(blocks_brain, nodes[index], BLOCKS, 0)

                project_map = Dict(nodes[index] => Set([nodes[index],BLOCKS]))
                project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)

                block = get_block_index(blocks_brain, BLOCKS, 0.75)

                verbose && println("Block ", block, " is top block.")
                blocks_brain.brain.no_plasticity = false
                # Postactions
                inhibit_area(blocks_brain, BLOCKS, 0)
                inhibit_area(blocks_brain, nodes[index], 0)
                inhibit_fiber(blocks_brain, BLOCKS, nodes[index], 0)
                return (index, block)
            end
            # Postactions
            inhibit_area(blocks_brain, nodes[index], 0)
            inhibit_fiber(blocks_brain, head, nodes[index], 0)
            
        end
    end
    inhibit_areas(blocks_brain, vcat([BLOCKS,RELOCATED], working_regions), 0)
    inhibit_all_fibers(blocks_brain, vcat([BLOCKS,RELOCATED], working_regions), 0)
    blocks_brain.brain.no_plasticity = false
    return (0, 0)
end


# Delete pointer from HEAD to old area by creating pointer from
# HEAD to next area (if it exists). It first projects the new
# area to HEAD and then reinforce their connection

# Parameters
# - blocks_brain: a BlocksBrain representing stacks of blocks
# - stack_index: the index of the stack we look at
# - prefix: string "I", "G", "B", or "T" whether the areas we work on are the input areas, the goal areas, the building areas, or the table areas.
# Example: pop!(bb,4,"G")
function pop!(blocks_brain::BlocksBrain, stack_index::Int, prefix::String; verbose::Bool=false)::Nothing
    top_area, top_block = top(blocks_brain, stack_index, prefix)
    
    new_top_area = mod(top_area + 1, 1:MAX_NODES_AREAS)

    # defining locally the regions involved to parse stack number stack_index
    head = add_prefix([HEADS[stack_index]], prefix)[1]
    nodes = add_prefix(NODES[stack_index], prefix)

    if (verbose)
        println("Top area ", nodes[top_area])
        println("Top block ", top_block)
    end

    if is_last_block(blocks_brain,stack_index,top_area,top_block,prefix)
        blocks_brain.brain.areas[head].winners = Vector{Int}() 
        verbose && println("Top block is last block.")
        return
    end

    verbose && println("New top area ", nodes[new_top_area])

    blocks_brain.brain.no_plasticity = true
    for area in vcat([BLOCKS,head], nodes)
        if is_assembly(blocks_brain, area)
            unfix_assembly(blocks_brain.brain.areas[area])
        end
    end


    # LOOKING FOR THE CORRECT ASSEMBLY IN NEW AREA
    # Preactions
    disinhibit_area(blocks_brain, head, 0)
    disinhibit_area(blocks_brain, nodes[top_area], 0)
    disinhibit_fiber(blocks_brain, head, nodes[top_area], 0)
    project_map = Dict{String,Set{String}}()
    project_map[head] = Set{String}([head,nodes[top_area]])
    project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)
    # Postactions
    inhibit_area(blocks_brain, head, 0)
    inhibit_fiber(blocks_brain, head, nodes[top_area], 0)
    # Preactions
    disinhibit_area(blocks_brain, nodes[new_top_area], 0)
    disinhibit_fiber(blocks_brain, nodes[top_area], nodes[new_top_area], 0)
    project_map = Dict{String,Set{String}}()
    project_map[nodes[top_area]] = Set{String}([nodes[top_area],nodes[new_top_area]])
    project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)

    # Postactions
    inhibit_area(blocks_brain, nodes[top_area], 0)
    inhibit_fiber(blocks_brain, nodes[top_area], nodes[new_top_area], 0)

    # CREATING NEW ASSOCIATION BETWEEN NEW AREA AND head
    blocks_brain.brain.no_plasticity = false
    disinhibit_areas(blocks_brain, [head,BLOCKS], 0)
    disinhibit_fiber(blocks_brain, head, nodes[new_top_area], 0)
    disinhibit_fiber(blocks_brain, BLOCKS, nodes[new_top_area], 0)

    project_map = Dict{String,Set{String}}()
    project_map[nodes[new_top_area]] = Set{String}([nodes[new_top_area], head, BLOCKS])
    project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)
    #println("After popping, we will have top block ", get_block_index(blocks_brain,BLOCKS))

    # REINFORCING THE ASSOCIATION BETWEEN NEW AREA AND head
    proj_rounds = 40
    # Preactions
    project_star(blocks_brain; project_rounds=proj_rounds, verbose=false)
    # Postactions
    inhibit_areas(blocks_brain, [BLOCKS,head], 0)
    inhibit_area(blocks_brain, nodes[new_top_area], 0)
    inhibit_fibers(blocks_brain, [[nodes[new_top_area], head],[nodes[new_top_area], BLOCKS]], 0)
    for area in vcat([BLOCKS,head], nodes)
        if is_assembly(blocks_brain, area)
            unfix_assembly(blocks_brain.brain.areas[area])
        end
    end
    return 
end

# Adding pointer head from new area by projecting new block first to new area and then to head

# Parameters
# - blocks_brain: a BlocksBrain representing stacks of blocks
# - stack_index: the index of the stack we look at
# - block: (Int) number of block to put on top of the stack
# - prefix: string "I", "G", "B", or "T" whether the areas we work on are the input areas, the goal areas, the building areas, or the table areas.
# Example: put!(bb,4,2,"G")
function put!(blocks_brain::BlocksBrain, stack_index::Int, block::Int, prefix::String; verbose::Bool=false)::Nothing
    top_node_index, top_block = top(blocks_brain, stack_index, prefix)
    

    blocks_brain.brain.no_plasticity = false

    # adding prefix to regions
    working_regions = add_prefix(vcat(REGIONS...), prefix)

    for area in vcat([BLOCKS], working_regions)
        unfix_assembly(blocks_brain.brain.areas[area])
    end

    # defining locally the regions involved to parse stack number stack_index
    head = add_prefix([HEADS[stack_index]], prefix)[1]
    nodes = add_prefix(NODES[stack_index], prefix)


    if (top_node_index == 0 || top_block == 0)
        verbose && println("We begin with a new block!")
        new_node_index = MAX_NODES_AREAS
    else
        new_node_index = mod(top_node_index - 1, 1:MAX_NODES_AREAS)
        if (verbose)
            println("Top area ", nodes[top_node_index])
            println("Top block ", top_block)
            println("New top area ", nodes[new_node_index])
        end
    end
    disinhibit_areas(blocks_brain, [nodes[new_node_index],BLOCKS, head], 0)
    disinhibit_fibers(blocks_brain, [[BLOCKS,nodes[new_node_index]],[nodes[new_node_index],head]], 0)

    # create new assembly in area nodes[new_node_index] from BLOCKS
    activate_block(blocks_brain, block)
    project_map = Dict{String,Set{String}}()
    project_map[BLOCKS] = Set{String}([BLOCKS,nodes[new_node_index]])
    project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)

    # create new assembly in head from nodes[new_node_index]
    project_map = Dict{String,Set{String}}()
    project_map[nodes[new_node_index]] = Set{String}([head,nodes[new_node_index]])
    project(blocks_brain.brain, Dict{String,Set{String}}(), project_map)

    proj_rounds = 50
    if (top_node_index != 0 && top_block != 0)
        disinhibit_area(blocks_brain, nodes[top_node_index], 0)
        disinhibit_fiber(blocks_brain, nodes[new_node_index], nodes[top_node_index], 0)
    end

    activate_block(blocks_brain, block)
    project_star(blocks_brain; project_rounds=proj_rounds)

    for area in vcat([BLOCKS,head], nodes)
        if is_assembly(blocks_brain, area)
            unfix_assembly(blocks_brain.brain.areas[area])
        end
    end

    inhibit_areas(blocks_brain, vcat([BLOCKS], working_regions), 0)
    inhibit_all_fibers(blocks_brain, vcat([BLOCKS], working_regions), 0)
    return
end

# this function is necessary for dismantle! and pop! and checks if a given block is the last one
## is_last_block
# Parameters
# - blocks_brain: type BlocksBrain, after having parsed the stacks
# - stacks_index: an integer, it is the index of stack we investigate
# - node_index: an integer, the NODE area the block is projected into
# - block: an integer, the block we check
# - prefix: string "I", "G", "B", or "T" whether the areas we work on are the input areas, the goal areas, the building areas, or the table areas.
# Example: is_last_block(blocks_brain, 5, 2, 3, "I")
function is_last_block(blocks_brain::BlocksBrain, stack_index::Int, node_index:: Int, block::Int, prefix::String; verbose::Bool = false)::Bool
    # defining locally the regions involved
    nodes = add_prefix(NODES[stack_index],prefix)
    blocks_brain.brain.no_plasticity = true

    #preactions
    next_index = mod(node_index + 1, 1:MAX_NODES_AREAS)
    disinhibit_areas(blocks_brain, [nodes[node_index],nodes[next_index],BLOCKS],0)
    disinhibit_fibers(blocks_brain, [[nodes[node_index],nodes[next_index]],[BLOCKS,nodes[node_index]]], 0 )

    for area in vcat([BLOCKS],nodes)
        unfix_assembly(blocks_brain.brain.areas[area])
    end

    activate_block(blocks_brain,block)
    project_map = Dict(BLOCKS => Set([nodes[node_index],BLOCKS]))
    project(blocks_brain.brain,Dict{String,Set{String}}(),project_map)

    project_map = Dict(nodes[node_index] => Set([nodes[node_index],nodes[next_index]]))
    project(blocks_brain.brain,Dict{String,Set{String}}(),project_map)

    is_last = !(is_assembly(blocks_brain,nodes[next_index]))

    unfix_assembly(blocks_brain.brain.areas[BLOCKS])
    inhibit_areas(blocks_brain, [nodes[node_index],nodes[next_index]],0)
    inhibit_fiber(blocks_brain, nodes[node_index],nodes[next_index], 0 )
    blocks_brain.brain.no_plasticity = false

    return is_last

end
