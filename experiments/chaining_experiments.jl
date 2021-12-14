# run with `julia -t auto src/large_chaining_experiments.jl`

#include("../src/brain.jl")
#include("../src/blocks_brain.jl")
include("experiment_functions.jl")

using BSON, Dates, ProgressMeter

# Test limits on chaining with 2 areas using unidirectional fibers
# first block: blocks -> A
# second block: A, blocks -> B
# third block: B, blocks -> A
# i-th block: {i even: A ; odd: B} , blocks -> {i even: B ; odd: A}
# in the end, given first block assembly, sucesss = can read out entire stack going BLOCKS->A->B->A etc. (and detect end)


function two_area_chaining(chain_length, p, eak, nean, neak, db, proj_rounds)
	blocks_brain = BlocksBrain(blocks_number, [FIRST, SECOND], p, eak, nean, neak, db)
	disinhibit_area(blocks_brain, BLOCKS, 0)
end

# Test limits on chaining with 3 areas using bidirectional fibers
# first block: blocks -> A
# second block: A, blocks -> B
# third block: B, blocks -> C
# fourth block: C, blocks -> A
# i-th block: {i mod 3 = 2: A ; 0: B ; 1: C}, blocks-> {i mod 3 = 2: B ; 0: C ; 1: A}

global const THIRD = "THIRD"

# three_area_chaining(3, 0.1, 50, 1000000, 50, 0.1, 40), succeeds
# three_area_chaining(5, 0.1, 50, 1000000, 50, 0.1, 40), doesn't succeed
function three_area_chaining(chain_length::Int64, p, eak::Int64, nean::Int64, neak::Int64, db, proj_rounds::Int64; verbose=false, blocks_recurse=false)
	blocks_brain = BlocksBrain(chain_length, [FIRST, SECOND, THIRD], p, eak, nean, neak, db)

	# Set up chain by projecting first block into FIRST
	disinhibit_area(blocks_brain, BLOCKS, 0)
	disinhibit_area(blocks_brain, FIRST, 0)
	disinhibit_fiber(blocks_brain, BLOCKS, FIRST, 0)
	disinhibit_fiber(blocks_brain, FIRST, FIRST, 0)
	inhibit_fiber(blocks_brain, BLOCKS, BLOCKS, 0)
	activate_block(blocks_brain, 1)
	project_star(blocks_brain; project_rounds=proj_rounds, verbose=verbose)

	for i in 2:chain_length
		activate_block(blocks_brain, i)
		current_area = FIRST
		next_area = SECOND
		other_area = THIRD
		if mod(i, 3) == 0
			current_area = SECOND
			next_area = THIRD
			other_area = FIRST
		elseif mod(i, 3) == 1
			current_area = THIRD
			next_area = FIRST
			other_area = SECOND
		end
		inhibit_area(blocks_brain, other_area, 0)
		disinhibit_area(blocks_brain, current_area, 0)
		disinhibit_area(blocks_brain, next_area, 0)
		disinhibit_fiber(blocks_brain, current_area, next_area, 0)
		disinhibit_fiber(blocks_brain, BLOCKS, next_area, 0)
		# disinhibit_fiber(blocks_brain, next_area, next_area, 0)
		# inhibit_fiber(blocks_brain, current_area, current_area, 0)  # prevent over-reinforcement?
		inhibit_fiber(blocks_brain, BLOCKS, current_area, 0)

		project_star(blocks_brain; project_rounds=proj_rounds, verbose=verbose)
	end

	# Readout to test success
	# How many assemblies out of chain_length are indeed assemblies
	success_assemblies = 0
	# How many out of chain_length blocks are correctly restored
	success_blocks = 0

	blocks_brain.brain.no_plasticity = true
	for area in [BLOCKS, FIRST, SECOND, THIRD]
        unfix_assembly(blocks_brain.brain.areas[area])
    end
	activate_block(blocks_brain, 1)
	project(blocks_brain.brain, Dict{String,Set{String}}(), Dict{String,Set{String}}(BLOCKS => Set([FIRST])))
	# project from 1 -> FIRST (test is_assembly ?)
	if is_assembly(blocks_brain, FIRST)
		success_assemblies += 1
	end
	unfix_assembly(blocks_brain.brain.areas[BLOCKS])
	project(blocks_brain.brain, Dict{String,Set{String}}(), Dict{String,Set{String}}(FIRST => Set([BLOCKS])))
	if blocks_recurse
		project(blocks_brain.brain, Dict{String,Set{String}}(), Dict{String,Set{String}}(BLOCKS => Set([BLOCKS])))
	end
	
	answer = get_block_index(blocks_brain, BLOCKS, 0.75)
	if answer == 1
		success_blocks += 1
	end

	current_area = FIRST
	next_area = SECOND

	for i in 2:chain_length
		project(blocks_brain.brain, Dict{String,Set{String}}(), Dict{String,Set{String}}(current_area => Set([next_area])))
		if is_assembly(blocks_brain, next_area, 0.5)
			success_assemblies += 1
		end
		project(blocks_brain.brain, Dict{String,Set{String}}(), Dict{String,Set{String}}(next_area => Set([BLOCKS])))
		if blocks_recurse
			project(blocks_brain.brain, Dict{String,Set{String}}(), Dict{String,Set{String}}(BLOCKS => Set([BLOCKS])))
		end
		answer = get_block_index(blocks_brain, BLOCKS, 0.75)
		if answer == i
			success_blocks += 1
		end
		if current_area == FIRST
			current_area = SECOND
			next_area = THIRD
		elseif current_area == SECOND
			current_area = THIRD
			next_area = FIRST
		elseif current_area == THIRD
			current_area = FIRST
			next_area = SECOND
		end
	end

	# Could add final test: now project CURRENT_AREA -> NEXT_AREA, test is_assembly, should be false
	verbose && println("Good assemblies obtained for ", success_assemblies, "/", chain_length)
	verbose && println("Correct blocks obtained for ", success_blocks, "/", chain_length)

	return success_assemblies, success_blocks
end


# 1) Fix parameters, try longer chains, get percentage success (or map against y=x to show drop off)
# exp_vary_chain_length(1, 40, 0.1, 50, 100000, 50, 0.1, 40, 1, 10)
# exp_vary_chain_length(1, 40, 0.1, 50, 500000, 50, 0.1, 40, 1, 10)
# exp_vary_chain_length(1, 40, 0.1, 50, 1000000, 50, 0.1, 40, 1, 10)
function exp_vary_chain_length(min::Int64, max::Int64, p, eak::Int64, nean::Int64, neak::Int64, db, proj_rounds::Int64, step::Int64, repeat::Int64)
	results = []
	results_percent = []
	chain_lengths = []
	assemblies = []
	chain_length = min
	while chain_length < max
		success_blocks_sum = 0
		success_assemblies_sum = 0
		for _ in 1:repeat
			success_assemblies, success_blocks = three_area_chaining(chain_length, p, eak, nean, neak, db, proj_rounds)
			success_blocks_sum += success_blocks
			success_assemblies_sum += success_assemblies
		end
		success_blocks_avg = success_blocks_sum / repeat
		success_assemblies_avg = success_assemblies_sum / repeat
		append!(chain_lengths, chain_length)
		append!(results, success_blocks_avg)
		append!(results_percent, (success_blocks_avg / chain_length) * 100.0)
		append!(assemblies, success_assemblies_avg)
		chain_length += step
	end

	return chain_lengths, results, results_percent, assemblies
end

# n=100000,k=50, from 1 to 31, steps of 2
# (Any[1, 3, 5, 7, 2, 11, 13, 0, 1, 2, 11, 5, 3, 5, 5], Any[100.0, 100.0, 100.0, 100.0, 22.22222222222222, 100.0, 100.0, 0.0, 5.88235294117647, 10.526315789473683, 52.38095238095239, 21.73913043478261, 12.0, 18.51851851851852, 17.24137931034483])
# n=1000000,k=50, starting from 30 (steps of 3??)
# (Any[30, 14, 32, 5, 34, 22, 1, 12, 12, 2], Any[100.0, 43.75, 94.11764705882352, 13.88888888888889, 89.47368421052632, 55.00000000000001, 2.380952380952381, 27.27272727272727, 26.08695652173913, 4.166666666666666])

# 2) Fix n, fix chain length (5,10?) vary k, find "window" of success
# exp_find_k_range(10, 0.1, 100000, 0.1, 40)
function exp_find_k_range(chain_length, p, n, db, proj_rounds; verbose = false)
	k = 10
	reached_lower_bound = false
	lower_bound = -1
	upper_bound = -1
	while k < (n/2)
		verbose && println("Testing k = ", k)
		success_assemblies, success_blocks = three_area_chaining(chain_length, p, k, n, k, db, proj_rounds)
		if !reached_lower_bound && success_blocks == chain_length
			reached_lower_bound = true
			lower_bound = k
		elseif reached_lower_bound && success_blocks < chain_length
			upper_bound = k 
			break
		end
		k *= 2
	end
	return lower_bound, upper_bound
end

# 3) Fix n, vary k, find max chain length with 100% success
# exp_vary_k_max_chain(10, 180, 0.1, 100000, 0.1, 40, 10)
# exp_vary_k_max_chain(10, 180, 0.1, 500000, 0.1, 40, 10)
# exp_vary_k_max_chain(10, 180, 0.1, 1000000, 0.1, 40, 10)
function exp_vary_k_max_chain(lower_k, higher_k, p, nean, db, proj_rounds, repeat; verbose = false)
	all_k = []
	max_chain_lengths = []
	k = lower_k
	while k < higher_k
		verbose && println("Testing k = ", k)
		append!(all_k, k)
		sum = 0
		for _ in 1:repeat
			chain_length = 1
			while true
				_, success_blocks = three_area_chaining(chain_length, p, k, nean, k, db, proj_rounds)
				if success_blocks < chain_length 
					chain_length -= 2  # succeeded at previous chain_length
					break
				end
				chain_length += 2
			end
			sum += chain_length
		end
		average_max_length = sum / repeat
		verbose && println("Got max chain=", average_max_length, " for k=", k)
		append!(max_chain_lengths, average_max_length)
		k += 10
	end
	return all_k, max_chain_lengths
end

const varychainargs = Tuple{Int64, Int64, Float64, Int64, Int64, Int64, Float64, Int64, Int64, Int64}
runexperim(args::varychainargs) = exp_vary_chain_length(args...)
const varykargs = Tuple{Int64, Int64 , Float64, Int64, Float64, Int64, Int64}
runexperim(args::varykargs) = exp_vary_k_max_chain(args...)

params = Vector{Union{varychainargs,varykargs}}()
append!(params, [
	(1, 40, 0.1, 50, 100000, 50, 0.1, 40, 1, 1),
	(1, 40, 0.1, 50, 500000, 50, 0.1, 40, 1, 1),
	(1, 40, 0.1, 50, 1000000, 50, 0.1, 40, 1, 1),
	(10, 180, 0.1, 100000, 0.1, 40, 1),
	(10, 180, 0.1, 500000, 0.1, 40, 1),
	(10, 180, 0.1, 1000000, 0.1, 40, 1)
	]
)
runs = [[params[j],i]  for j in 1:length(params) for i in 1:50]

results = Dict{Any,Any}()

numruns = length(runs)
p = Progress(numruns);
update!(p,0)
counter = Threads.Atomic{Int}(0)
l = Threads.SpinLock()

Threads.@threads for i in 1:numruns
	Threads.atomic_add!(counter, 1)
	Threads.lock(l)
	update!(p, counter[])
	Threads.unlock(l)  
	
	results[runs[i]] = runexperim(runs[i][1]);
end

bson("experiments_"*string(now())*".bson", Dict(:results => results))


