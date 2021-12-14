include("../src/brain.jl")
include("../src/blocks_brain.jl")
include("../src/utils.jl") 

maxstacks = (@isdefined MAXSTACKS) ? MAXSTACKS : 5 ## we want only a limited number of stacks
const limitconfig = 3 # in the end uses limitconfig^2
const pathname = "block_world_instances/"
const blocksizes = [10, 20, 30]

# aliases
const ConfigType = Vector{Vector{Int64}}
const RunType = NamedTuple{(:blocks, :initial, :goal, :p, :eak, :nean, :neak, :db, :stacks), Tuple{Int, ConfigType, ConfigType, Float64, Int, Int, Int, Float64, Int}}
 
"""
	readconfigs(blocksizes...)

	#TODO...
Read the configurations as provided by the Java generator.
It assumes that the latter files are located in `pathname` and named `bw_\$BLOCKSIZE.txt`, where `BLOCKSIZE`∈`blocksizes`. 
"""
function readconfigs(blocksizes::Vector{Int} = [10], pathname::String =  pathname, max_stacks::Int = maxstacks)::Dict{Int, Vector{ConfigType}}
	configurations = Dict{Int, Vector{ConfigType}}()
	for i in blocksizes
		filepath = pathname*"bw_"*string(i)*".txt"
		configurations[i] = converter(filepath, max_stacks)
	end
	return configurations
end

# TODO for parameters should provide ranges instead.. (also blocksizes)
"""
	generateruns(blocksizes...)

Assumes `limitconfig` is no more than the number of configurations available for each blocksize. 
"""
function generateruns(blocksizes = blocksizes, prange = 0.1:0.1, eakrange = 50:50, neanrange = 10^6:10^6, neakrange = 50:50, dbrange = 0.1:0.1, pathname = pathname, limitconfig = limitconfig, maxstacks = maxstacks)::Vector{RunType}
	runslist = Vector{RunType}()
	# importing instances 
	configurations = readconfigs(blocksizes, pathname, maxstacks)
	# creating runs  
	for blocksize in blocksizes,  p in prange, eak in eakrange, nean in neanrange, neak in neakrange, db in dbrange
		(limitconfig > length(configurations[blocksize])) && throw(error("Argument `limitconfig` is $(limitconfig), which is greater than the number of configurations $(length(configurations[blocksize]))!"))
		for i in 1:limitconfig, j in 1:limitconfig 
			append!(runslist, [
				( 
					blocks = blocksize, 
					initial = configurations[blocksize][i], 
					goal = configurations[blocksize][j], 
					p = p, 
					eak = eak, 
					nean = nean, 
					neak = neak, 
					db = db, 
					stacks = maxstacks
				)
			])
		end
	end
	return runslist
end

function printruns(runslist::Vector{RunType})
	for i in runslist
		for j in i
			print("\"",j,"\""," ")
		end
		println()
	end
end

# printruns(generateruns())