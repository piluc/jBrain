using BSON
using DelimitedFiles
using Random

savebrain(filename::String, bb::BlocksBrain) = bson(filename, Dict(:bb => bb))

loadbrain(filename::String) = BSON.load(filename)

function execute_pops(before_pops::Vector{Vector{Int}}, pops::Vector{Vector{String}}; verbose::Bool = false)
    number_of_stacks = length(before_pops)
    for s in 1:number_of_stacks
        for b in 1:length(pops[s])
            block_to_be_popped = parse(Int64, pops[s][b])
            if (block_to_be_popped == before_pops[s][1])
                x = [block_to_be_popped]
                append!(before_pops, [x])
                verbose && println("Remove block ", block_to_be_popped, " from input stack ", s, ".")
                deleteat!(before_pops[s], 1)
            else
                println("Error!")
            end
        end
    end
    return before_pops
end

function execute_puts(before_puts::Vector{Vector{Int}}, puts::Vector{Vector{String}}; verbose::Bool = true)::Vector{Vector{Int}}
    number_of_stacks = length(before_puts)
    number_of_put_vectors = length(puts)
    after_puts = Vector{Vector{Int}}(undef, number_of_put_vectors)
    for p in 1:number_of_put_vectors
        bottom_block = parse(Int64, puts[p][1])
        verbose && println("Fix the stack with top block ", bottom_block, " as final stack ", p, ".")
        # indeces_to_delete = []
        stack_index = 0
        for s in 1:number_of_stacks
            if before_puts[s][1] == bottom_block
                stack_index = s
                after_puts[p] = deepcopy(before_puts[s])
                deleteat!(before_puts, s)
                number_of_stacks = number_of_stacks - 1
                break
            end
        end
        if (stack_index == 0)
            println("Error stack_index: ", bottom_block)
            println("Current configuration: ", before_puts)
            return
        end
        for b in 2:length(puts[p])
            block = parse(Int64, puts[p][b])
            verbose && println("Put block ", block , " on top of final stack ", p, ".")
            stack_to_be_deleted = 0
            for s in 1:number_of_stacks
                if before_puts[s][1] == block
                    stack_to_be_deleted = s
                    # append!(indeces_to_delete, [stack_to_be_deleted])
                    break
                end
            end    
            if (stack_to_be_deleted == 0)
                println("Error stack_to_be_deleted: ", block)
                println(before_puts)
                return
            end
            pushfirst!(after_puts[p], block)
            # pushfirst!(before_puts[stack_index], block)
            deleteat!(before_puts, stack_to_be_deleted)
            number_of_stacks = number_of_stacks - 1
        end
        # indeces_to_keep = [index for index in 1:(number_of_stacks) if !(index in indeces_to_delete) ]
        # number_of_stacks = number_of_stacks - length(indeces_to_delete)
        # before_puts = before_puts[indeces_to_keep]
    end
    if length(before_puts) > 0
        error("Something has been kept on the table have been deleted: ", before_puts)
    end
    return after_puts
end

## we take an array (d_i) 1 ≤ i ≤ n which represents a bw configuration in the following way
## d_i = 0 ---> block i on the table
## d_i = j > 0 ----> block i is right on block j
## we write a convert function which gives as output the bw config in our formattation
function Base.convert(input_config::Vector{Int}, stacks_number::Int, first_blocks::Vector{Int})::Vector{Vector{Int}}
    output_config = Vector{Vector{Int}}(undef, stacks_number)
    for i in 1:stacks_number
        output_config[i] = [first_blocks[i]]
        next_block = findall(x -> x == first_blocks[i], input_config)
        while (length(next_block) == 1)
            output_config[i] = append!([next_block[1]], output_config[i])
            next_block = findall(x -> x == next_block[1], input_config)
        end
    end
    return output_config
end


## we take in input a file having as lines vectors of configuration to be converted
## converts only those with 5 stacks max
function converter(input_file::String, max_stacks::Int)::Vector{Vector{Vector{Int}}}
    bw_instances = Vector{Vector{Vector{Int}}}()
    open(input_file, "r") do io
        input_configs = readdlm(io, ',', Int)
        num_raws, num_cols = size(input_configs)
        for i in 1:num_raws
            input_config = input_configs[i,:]
            first_blocks = findall(x -> x == 0, input_config)
            stacks_number = length(first_blocks)
            if stacks_number < max_stacks + 1
                append!(bw_instances, [convert(input_config, stacks_number, first_blocks)])
            end
        end 
    end
    return bw_instances
end

function prune_configuration(bw_instances::Vector{Vector{Vector{Int}}}, max_blocks)::Vector{Vector{Vector{Int}}}
    pruned_bw_instances = Vector{Vector{Vector{Int}}}()
    for c in 1:length(bw_instances)
        to_be_included = true
        for s in 1:length(bw_instances[c])
            if (length(bw_instances[c][s]) > max_blocks)
                to_be_included = false
            end
        end
        if (to_be_included)
            append!(pruned_bw_instances, [bw_instances[c]])
        end
    end
    return pruned_bw_instances
end

function param_tests(blocks_num ::Int = 10, optimized::Bool = false, max_stacks::Int = 5, max_inputs::Int = 10, max_blocks::Int = 7)::Nothing
    to_add = ""
    if !(optimized)
       to_add = "non_"
    end
    input_file::String = string("./bw_instances/data_to_create_instances/bw_",blocks_num,".txt")
    output_file =string("./bw_instances/data_to_create_instances/params_",blocks_num,"blocks_",to_add,"optimized.txt")
    
    inputs = prune_configuration( converter(input_file, max_stacks),max_blocks)
    num_lines = length(inputs)
    open(output_file, "w") do io
        mystring = ""
        for _ in 1:max_inputs
            for _ in 1:max_inputs
                tmp_i = rand(1:num_lines)
                indeces = Vector{Int}(1:num_lines)
                deleteat!(indeces, tmp_i)
                tmp_j = rand(indeces)
                #tmp_inputs = deepcopy(inputs[indeces])
                if (optimized)
                    mystring = string(mystring, "\"optimized\" \"", inputs[tmp_i], "\" \"", inputs[tmp_j], "\"\n")
                else
                    mystring = string(mystring, "\"non_optimized\" \"", inputs[tmp_i], "\" \"", inputs[tmp_j], "\"\n")
                end
            end
        end
        write(io, mystring)
    end
    return 
end

maxmem = 0
function updatemaxmem()
    s = readchomp(`ps -p $(getpid()) -o vsize`)
    ms = split(s, "\n")[2]
    mem = parse(Int, ms)
    global maxmem = max(maxmem, mem)
end