include("../src/brain.jl")
include("../src/blocks_brain.jl")
include("../src/bw_apps.jl")
include("../src/bw_planning.jl")
include("../src/bw_planning_optimized.jl")
include("../src/utils.jl")


## example
## test_parse_and_read([[1,3,5,8,9,11],[2,4],[6,7,10,12,13]],0.1,50,1000000,50,0.1)
function test_parse_and_read(stacks::Vector{Vector{Int}}, p::Float64, eak::Int, nean::Int, neak::Int, db::Float64)::Nothing
    stacks_number = length(stacks)
    stacks_lengths = Vector{Int}()
    blocks_number = 0
    top_areas = ones(Int, stacks_number)

    for stack in stacks
        stack_length = length(stack)
        append!(stacks_lengths, [stack_length])
        blocks_number = blocks_number + stack_length
    end

    prefix = "I"
    working_regions = add_prefix(vcat(REGIONS...), prefix)
    append!(working_regions, [RELOCATED])

    println("Creating the brain.")
    blocks_brain = brain_creation(blocks_number, working_regions, p, eak, nean, neak, db)
    println("Creation completed.")
    println("Parsing the stacks.")
    parse!(blocks_brain, stacks, prefix)
    println("Parsing completed.")
    println("Reading the brain.")
    stacks_read = readout(blocks_brain, stacks_number, stacks_lengths, top_areas, prefix)
    count = 1
    for stack in stacks_read
        println("Stacks number ", count, ":")
        count = count + 1
        println(stack)
    end
    return
end

## example
## test_top([[1,3,5,8,9,11],[2,4],[6,7,10,12,13]],0.1,50,1000000,50,0.1,"T")
function test_top(stacks::Vector{Vector{Int}}, p::Float64, eak::Int, nean::Int, neak::Int, db::Float64, prefix::String)::Nothing
    stacks_number = length(stacks)
    stacks_lengths = Vector{Int}()
    blocks_number = 0

    for stack in stacks
        stack_length = length(stack)
        append!(stacks_lengths, [stack_length])
        blocks_number = blocks_number + stack_length
    end

    working_regions = add_prefix(vcat(REGIONS...), prefix)
    append!(working_regions, [RELOCATED])

    println("Creating the brain.")
    blocks_brain = brain_creation(blocks_number, working_regions, p, eak, nean, neak, db)
    println("Creation completed.")

    println("Building the stacks.")
    for j in 1:stacks_number
        for block_index in stacks_lengths[j]:-1:1
            put!(blocks_brain, j, stacks[j][block_index], prefix)
        end
    end

    println("You want to relocate some blocks?")
    to_relocate = readline()
    if to_relocate == "yes"

        blocks_to_relocate = [1,5,8]
        stack_index_to_relocate = 1
        for block_to_relocate in blocks_to_relocate
            set_as_relocated!(blocks_brain, stack_index_to_relocate, block_to_relocate, prefix)
            println("We now relocate block ", block_to_relocate, " in stack ", stack_index_to_relocate, ".")
        end
    end

    println("Finding the tops.")
    top_areas = Vector{Int}(undef, stacks_number)
    for j in 1:stacks_number
        top_areas[j], block = top(blocks_brain, j, prefix)
    end
    println("Tops are ", top_areas, ".")
    println("Now reading the stacks.")

    stacks_read = readout(blocks_brain, stacks_number, stacks_lengths, top_areas, prefix)
    count = 1
    for stack in stacks_read
        println("Stacks number ", count, ":")
        count = count + 1
        println(stack)
    end
    return
end

## example
## test_pop([[1,3,5,8,9,11],[2,4],[6,7,10,12,13]],0.1,50,1000000,50,0.1; verbose = true)
function test_pop(stacks::Vector{Vector{Int}},p::Float64, eak::Int, nean::Int, neak::Int, db::Float64; verbose::Bool = false)::Nothing
    stacks_number = length(stacks)
    stacks_lengths = Vector{Int}()
    blocks_number = 0
    for stack in stacks
        stack_length = length(stack)
        append!(stacks_lengths, [stack_length])
        blocks_number = blocks_number + stack_length
    end
    source_prefix = "I"
    working_regions = add_prefix(vcat(REGIONS...), source_prefix)
    table_prefix = "T"
    append!(working_regions, add_prefix(vcat(REGIONS...), table_prefix))
    append!(working_regions, [RELOCATED])

    verbose && println("Creating the brain.")
    blocks_brain = brain_creation(blocks_number, working_regions, p, eak, nean, neak, db)
    verbose && println("Creation completed.")
    verbose && println("Parsing the stacks.")
    parse!(blocks_brain, stacks, source_prefix)
    verbose && println("Parsing completed.")

    verbose && println("Initial configuration: ")
    read_block_world_instance(blocks_brain,stacks_number,"I")
    read_block_world_instance(blocks_brain,stacks_number,"T")

    verbose && println("Popping from the stacks.")
    for i in 1:stacks_number
        top_node_index, top_block = top(blocks_brain,i,source_prefix)
        pop!(blocks_brain, i, source_prefix)
        put!(blocks_brain, i, top_block, table_prefix)
    end

    verbose && println("Final configuration: ")
    read_block_world_instance(blocks_brain,stacks_number,source_prefix)
    read_block_world_instance(blocks_brain,stacks_number,table_prefix)
    return
end

## example
## test_put([[11],[14],[16]],[[1,3,5,8,9,15],[2,4,17],[6,7,10,12,13]],0.1,50,1000000,50,0.1; verbose = true)
function test_put(table::Vector{Vector{Int}},stacks::Vector{Vector{Int}},p::Float64, eak::Int, nean::Int, neak::Int, db::Float64; verbose::Bool = false)::Nothing
    stacks_number = length(stacks)
    stacks_lengths = Vector{Int}()
    blocks_number = 0
    for stack in stacks
        stack_length = length(stack)
        append!(stacks_lengths, [stack_length])
        blocks_number = blocks_number + stack_length
    end
    blocks_number = blocks_number + length(table)
    source_prefix = "I"
    working_regions = add_prefix(vcat(REGIONS...), source_prefix)
    table_prefix = "T"
    append!(working_regions, add_prefix(vcat(REGIONS...), table_prefix))
    append!(working_regions, [RELOCATED])

    verbose && println("Creating the brain.")
    blocks_brain = brain_creation(blocks_number, working_regions, p, eak, nean, neak, db)
    verbose && println("Creation completed.")
    verbose && println("Parsing the stacks.")
    parse!(blocks_brain, stacks, source_prefix)
    verbose && println("Parsing completed.")
    verbose && println("Parsing the table.")
    parse!(blocks_brain, table, table_prefix)
    verbose && println("Parsing completed.")


    verbose && println("Initial configuration: ")
    read_block_world_instance(blocks_brain,stacks_number,source_prefix)
    read_block_world_instance(blocks_brain,stacks_number,table_prefix)
    
    verbose && println("Putting the blocks on top.")
    for i in 1:stacks_number
        put!(blocks_brain, i, table[i][1], source_prefix)
        set_as_relocated!(blocks_brain,i,table[i][1],table_prefix)
    end

    verbose && println("Final configuration: ")
    read_block_world_instance(blocks_brain,stacks_number,source_prefix)
    read_block_world_instance(blocks_brain,stacks_number,table_prefix)
    return
end

# ## example
# test_intersect([[1,2],[3,4]],[[2,3],[1,4]],0.1,50,4*10^6,50,0.1; verbose = true)
function test_intersect(first_stacks::Vector{Vector{Int}}, second_stacks::Vector{Vector{Int}},p::Float64, eak::Int, nean::Int, neak::Int, db::Float64; to_readout::Bool = false, verbose::Bool = false)::Nothing
    input_stacks_number = length(first_stacks)
    input_stacks_lengths = Vector{Int}()
    input_blocks_number = 0
    input_top_areas = ones(Int, input_stacks_number)

    for input_stack in first_stacks
        input_stack_length = length(input_stack)
        append!(input_stacks_lengths, [input_stack_length])
        input_blocks_number = input_blocks_number + input_stack_length
    end

    goal_stacks_number = length(second_stacks)
    goal_stacks_lengths = Vector{Int}()
    goal_blocks_number = 0
    goal_top_areas = ones(Int, goal_stacks_number)
    building_top_areas = Vector{Int}()

    for goal_stack in second_stacks
        goal_stack_length = length(goal_stack)
        append!(goal_stacks_lengths, [goal_stack_length])
        append!(building_top_areas, mod(MAX_NODES_AREAS + 1 - goal_stack_length, 1:MAX_NODES_AREAS))
        goal_blocks_number = goal_blocks_number + goal_stack_length
    end

    input_prefix = "I"
    working_regions = add_prefix(vcat(REGIONS...), input_prefix)

    goal_prefix = "G"
    append!(working_regions, add_prefix(vcat(REGIONS...), goal_prefix))

    building_prefix = "B"
    append!(working_regions, add_prefix(vcat(REGIONS...), building_prefix))

    table_prefix = "T"
    append!(working_regions, add_prefix(vcat(REGIONS...), table_prefix))
    append!(working_regions, [RELOCATED,TMP])

    verbose && println("Creating the brain.")
    blocks_brain = brain_creation(input_blocks_number, working_regions, p, eak, nean, neak, db)
    verbose && println("Brain creation completed.")

    verbose && println("Parsing the fist stacks.")
    parse!(blocks_brain, first_stacks, input_prefix)
    verbose && println("Parsing completed.")

    verbose && println("Parsing the second stacks.")
    parse!(blocks_brain, second_stacks, goal_prefix)
    verbose && println("Parsing completed.")


    verbose && println("Intersecting the stacks pairwise.")
    for i in 1:min(input_stacks_number,goal_stacks_number)
        verbose && println("Intersecting the stacks in position ", i, ".")
        intersection, block = intersect(blocks_brain,i,input_prefix,i,goal_prefix)
        verbose && println("There is intersection:", intersection, ".")
        if (intersection)
            verbose && println("Intersection up to block ", block, ".")
        end
    end
    return
end

## example
## test_dismantle([[1,3,5,8,9,11],[2,4],[6,7,10,12,13]],0.1,50,1000000,50,0.1)
function test_dismantle(stacks::Vector{Vector{Int}}, p::Float64, eak::Int, nean::Int, neak::Int, db::Float64; to_readout::Bool = false, verbose::Bool = false)::Nothing
    stacks_number = length(stacks)
    stacks_lengths = Vector{Int}()
    blocks_number = 0

    for stack in stacks
        stack_length = length(stack)
        append!(stacks_lengths, [stack_length])
        blocks_number = blocks_number + stack_length
    end

    source_prefix = "I"
    working_regions = add_prefix(vcat(REGIONS...), source_prefix)
    table_prefix = "T"
    append!(working_regions, add_prefix(vcat(REGIONS...), table_prefix))
    append!(working_regions, [RELOCATED])

    println("Creating the brain.")
    blocks_brain = brain_creation(blocks_number, working_regions, p, eak, nean, neak, db)
    println("Creation completed.")
    println("Parsing the stacks.")
    parse!(blocks_brain, stacks, source_prefix)
    println("Parsing completed.")

    myactions = dismantle!(blocks_brain, stacks_number, source_prefix; to_readout, verbose)

    count = 1
    for stack_actions in myactions
        println("Actions for stack ", count, ":")
        println(stack_actions)
        print("\n")
        count = count + 1
    end
    return
end

## example
## test_planning([[1,3,5,8,9,11],[4,2],[6,7,10,12,13]],[[1,5,3,8,9,11],[2,4],[6,10,7,12,13]],0.1,50,1000000,50,0.1)
## test_planning([[1,3,5,8,6,7],[2,4]],[[1,2,3,4,5],[8,7,6]],0.1,50,10000000,50,0.1)
## test_planning([[1,3,5],[2,4]],[[1,2,3],[4,5]],0.1,50,10000000,50,0.1)
function test_planning(input_stacks::Vector{Vector{Int}}, goal_stacks::Vector{Vector{Int}}, p::Float64, eak::Int, nean::Int, neak::Int, db::Float64)::Nothing
    input_stacks_number = length(input_stacks)
    input_stacks_lengths = Vector{Int}()
    input_blocks_number = 0
    input_top_areas = ones(Int, input_stacks_number)

    for input_stack in input_stacks
        input_stack_length = length(input_stack)
        append!(input_stacks_lengths, [input_stack_length])
        input_blocks_number = input_blocks_number + input_stack_length
    end

    goal_stacks_number = length(goal_stacks)
    goal_stacks_lengths = Vector{Int}()
    goal_blocks_number = 0
    goal_top_areas = ones(Int, goal_stacks_number)
    building_top_areas = Vector{Int}()

    for goal_stack in goal_stacks
        goal_stack_length = length(goal_stack)
        append!(goal_stacks_lengths, [goal_stack_length])
        append!(building_top_areas, mod(MAX_NODES_AREAS + 1 - goal_stack_length, 1:MAX_NODES_AREAS))
        goal_blocks_number = goal_blocks_number + goal_stack_length
    end

    input_prefix = "I"
    working_regions = add_prefix(vcat(REGIONS...), input_prefix)

    goal_prefix = "G"
    append!(working_regions, add_prefix(vcat(REGIONS...), goal_prefix))

    building_prefix = "B"
    append!(working_regions, add_prefix(vcat(REGIONS...), building_prefix))

    table_prefix = "T"
    append!(working_regions, add_prefix(vcat(REGIONS...), table_prefix))
    append!(working_regions, [RELOCATED,TMP])

    println("Creating the brain.")
    blocks_brain = brain_creation(input_blocks_number, working_regions, p, eak, nean, neak, db)
    println("Brain creation completed.")

    println("Parsing the input stacks.")
    parse!(blocks_brain, input_stacks, input_prefix)
    println("Parsing completed. To double check, we run a readout.")
    input_stacks_read = readout(blocks_brain, input_stacks_number, input_stacks_lengths, input_top_areas, input_prefix)
    println("Readout completed.")
    for j in 1:input_stacks_number
        println("Stack number ", j, ":")
        println(input_stacks_read[j])
    end

    println("Parsing the goal stacks.")
    parse!(blocks_brain, goal_stacks, goal_prefix)
    println("Parsing completed. To double check, we run a readout.")
    goal_stacks_read = readout(blocks_brain, goal_stacks_number, goal_stacks_lengths, goal_top_areas, goal_prefix)
    println("Readout completed.")
    for j in 1:goal_stacks_number
        println("Stack number ", j, ":")
        println(goal_stacks_read[j])
    end

    println("Dismantling the input stacks.")
    dismantle_actions = dismantle!(blocks_brain, input_stacks_number, input_prefix; to_readout=false, verbose=false)
    println("Input stacks dismantled. The actions to perform follow.")
    for j in 1:input_stacks_number
        println("Stack number ", j, ":")
        println(dismantle_actions[j])
    end

    println("Reassembling the blocks correctly.")
    reassemble_actions = reassemble!(blocks_brain, input_stacks_number, goal_stacks_number, goal_prefix, building_prefix; to_readout=false, verbose=false)
    println("Blocks reassembled. The actions to perform follow.")
    for j in 1:goal_stacks_number
        println("Stack number ", j, ":")
        println(reassemble_actions[j])
    end

    println("Planning completed.")    

    println("To double check the correctness, a readout is run.")
    built_stacks_read = readout(blocks_brain, goal_stacks_number, goal_stacks_lengths, building_top_areas, building_prefix)
    println("Readout completed.")
    for j in 1:goal_stacks_number
        println("Stack number ", j, ":")
        println(built_stacks_read[j])
    end
    return
end


## example
# planning_experiment([[1,3,5],[4,2]],[[2,3,5],[1,4]],0.1,50,1000000,50,0.1; verbose = true)
# planning_experiment([[1,3,5,8,9,11],[4,2],[6,7,10,12,13]],[[1,5,3,8,9,11],[2,4],[6,10,7,12,13]],0.1,50,1000000,50,0.1; verbose = true)
function planning_experiment(input_stacks::Vector{Vector{Int}}, goal_stacks::Vector{Vector{Int}}, p::Float64, eak::Int, nean::Int, neak::Int, db::Float64; to_readout::Bool = false, verbose::Bool=false)::Tuple{Bool,Vector{Vector{Int}},Vector{Vector{Int}},Vector{Vector{String}},Vector{Vector{String}}}
    input_stacks_bkp = deepcopy(input_stacks)
    goal_stacks_bkp = deepcopy(goal_stacks)
    input_stacks_number = length(input_stacks)
    input_stacks_lengths = Vector{Int}()
    input_blocks_number = 0

    for input_stack in input_stacks
        input_stack_length = length(input_stack)
        append!(input_stacks_lengths, [input_stack_length])
        input_blocks_number = input_blocks_number + input_stack_length
    end

    goal_stacks_number = length(goal_stacks)
    goal_stacks_lengths = Vector{Int}()
    goal_blocks_number = 0
    building_top_areas = Vector{Int}()

    for goal_stack in goal_stacks
        goal_stack_length = length(goal_stack)
        append!(goal_stacks_lengths, [goal_stack_length])
        append!(building_top_areas, mod(MAX_NODES_AREAS + 1 - goal_stack_length, 1:MAX_NODES_AREAS))
        goal_blocks_number = goal_blocks_number + goal_stack_length
    end

    input_prefix = "I"
    working_regions = add_prefix(vcat(REGIONS...), input_prefix)
    goal_prefix = "G"
    append!(working_regions, add_prefix(vcat(REGIONS...), goal_prefix))
    building_prefix = "B"
    append!(working_regions, add_prefix(vcat(REGIONS...), building_prefix))
    table_prefix = "T"
    append!(working_regions, add_prefix(vcat(REGIONS...), table_prefix))
    append!(working_regions, [RELOCATED,TMP])
    verbose && println(input_stacks, " => ", goal_stacks)

    blocks_brain = brain_creation(input_blocks_number, working_regions, p, eak, nean, neak, db)
    verbose && println("Parsing the input stacks.")
    parse!(blocks_brain, input_stacks, input_prefix)
    verbose && println("Input stacks parsed. Parsing the goal stack.")
    parse!(blocks_brain, goal_stacks, goal_prefix)
    verbose && println("Goal stacks parsed.")
    verbose && println("Dismantling input stacks.")
    dismantle_actions = dismantle!(blocks_brain, input_stacks_number, input_prefix; to_readout, verbose=false)
    verbose && println("Input stacks dismantled.")
    verbose && println("Reassembling.")
    reassemble_actions = reassemble!(blocks_brain,input_stacks_number ,goal_stacks_number, goal_prefix, building_prefix; to_readout, verbose=false)
    verbose && println("Planning completed.")
    verbose && println(": ")
    verbose && println("The actions to perform follow.")
    after_pops_stacks = execute_pops(input_stacks_bkp, dismantle_actions; verbose)
    after_puts_stacks = execute_puts(after_pops_stacks, reassemble_actions; verbose)
    success = (after_puts_stacks == goal_stacks)
    if (verbose)
        if (success)
            print("true (pops: ", length(reduce(vcat, dismantle_actions)))
            println(", puts: ", (length(reduce(vcat, reassemble_actions)) - length(reassemble_actions)), ")")
        else
            println("false (pops: ", dismantle_actions, ", puts: ", reassemble_actions, ")")
            println("constructed stacks: ", after_puts_stacks)
        end
    end
    return (success, input_stacks, goal_stacks, dismantle_actions, reassemble_actions)
end

# optimized_planning_experiment([[1,3,5,8,9,11],[4,2],[6,7,10,12,13]],[[1,5,3,8,9,11],[2,4],[6,10,7,12,13]],0.1,50,1000000,50,0.1; verbose = true)
function optimized_planning_experiment(input_stacks::Vector{Vector{Int}}, goal_stacks::Vector{Vector{Int}}, p::Float64, eak::Int, nean::Int, neak::Int, db::Float64; to_readout::Bool = false, verbose::Bool=false)::Tuple{Bool,Vector{Vector{Int}},Vector{Vector{Int}},Vector{Vector{String}},Vector{Vector{String}}}
    input_stacks_bkp = deepcopy(input_stacks)
    goal_stacks_bkp = deepcopy(goal_stacks)
    input_stacks_number = length(input_stacks)
    input_stacks_lengths = Vector{Int}()
    input_blocks_number = 0

    for input_stack in input_stacks
        input_stack_length = length(input_stack)
        append!(input_stacks_lengths, [input_stack_length])
        input_blocks_number = input_blocks_number + input_stack_length
    end

    goal_stacks_number = length(goal_stacks)
    goal_stacks_lengths = Vector{Int}()
    goal_blocks_number = 0
    building_top_areas = Vector{Int}()

    for goal_stack in goal_stacks
        goal_stack_length = length(goal_stack)
        append!(goal_stacks_lengths, [goal_stack_length])
        append!(building_top_areas, mod(MAX_NODES_AREAS + 1 - goal_stack_length, 1:MAX_NODES_AREAS))
        goal_blocks_number = goal_blocks_number + goal_stack_length
    end

    input_prefix = "I"
    working_regions = add_prefix(vcat(REGIONS...), input_prefix)
    goal_prefix = "G"
    append!(working_regions, add_prefix(vcat(REGIONS...), goal_prefix))
    building_prefix = "B"
    append!(working_regions, add_prefix(vcat(REGIONS...), building_prefix))
    table_prefix = "T"
    append!(working_regions, add_prefix(vcat(REGIONS...), table_prefix))
    append!(working_regions, [RELOCATED,TMP])
    verbose && println(input_stacks, " => ", goal_stacks)
    blocks_brain = brain_creation(input_blocks_number, working_regions, p, eak, nean, neak, db)
    verbose && println("Parsing the input stacks.")
    parse!(blocks_brain, input_stacks, input_prefix)
    verbose && println("Input stacks parsed. Parsing the goal stack.")
    parse!(blocks_brain, goal_stacks, goal_prefix)
    verbose && println("Goal stacks parsed.")
    dismantle_actions, reassemble_actions = planning_optimized(blocks_brain, input_stacks_number, input_prefix, goal_stacks_number, goal_prefix, building_prefix; to_readout, verbose)
    verbose && println("Planning completed.")
    verbose && println(": ")
    verbose && println("The actions to perform follow.")
    after_pops_stacks = execute_pops(input_stacks_bkp, dismantle_actions; verbose)
    after_puts_stacks = execute_puts(after_pops_stacks, reassemble_actions; verbose)
    success = (after_puts_stacks == goal_stacks)
    if (verbose)
        if (success)        
            print("true (pops: ", length(reduce(vcat, dismantle_actions)))
            println(", puts: ", (length(reduce(vcat, reassemble_actions)) - length(reassemble_actions)), ")")
        else
            println("false (pop actions: ", dismantle_actions, ", ", reassemble_actions, ")")
            println("constructed stacks: ", after_puts_stacks)
        end
    end
    return (success, input_stacks, goal_stacks, dismantle_actions, reassemble_actions)
end

## spiking_experiment([[1,3,5,8,9,11],[2,4],[6,7,10,12,13]],0.1,50,1000000,50,0.1, 50)
function spiking_experiment(stacks::Vector{Vector{Int}}, p::Float64, eak::Int, nean::Int, neak::Int, db::Float64, pr::Int64)::Bool
    stacks_number = length(stacks)
    stacks_lengths = Vector{Int}()
    blocks_number = 0
    top_areas = ones(Int, stacks_number)

    for stack in stacks
        stack_length = length(stack)
        append!(stacks_lengths, [stack_length])
        blocks_number = blocks_number + stack_length
    end

    prefix = "I"
    working_regions = add_prefix(vcat(REGIONS...), prefix)
    append!(working_regions, [RELOCATED])

    blocks_brain = brain_creation(blocks_number, working_regions, p, eak, nean, neak, db)
    parse!(blocks_brain, stacks, prefix; project_rounds=pr)
    stacks_read = readout(blocks_brain, stacks_number, stacks_lengths, top_areas, prefix)
    return (stacks_read == stacks)
end

## examples
# spiking_binary_search([[1,3,5,8,9,11],[4,2],[6,7,10,12,13]], 0.1, 50, 1000000, 50, 0.1, 10)
function spiking_binary_search(stacks::Vector{Vector{Int}}, p::Float64, eak::Int, nean::Int, neak::Int, db::Float64, pr::Int64)::Vector{Int}
    lpr = 0
    upr = pr
    while (!spiking_experiment(stacks, p, eak, nean, neak, db, upr))
        lpr = upr
        upr = 2 * upr
    end
    while (lpr < upr)
        mpr = (lpr + upr) ÷ 2
        if (spiking_experiment(stacks, p, eak, nean, neak, db, mpr))
            upr = mpr
        else
            lpr = mpr + 1
        end
    end
    return [lpr, upr]
end