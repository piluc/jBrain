include("../experiments/experiment_functions.jl")

using ArgParse
function parse_commandline()
    settings = ArgParseSettings()
    settings.epilog = "Example:\n
        \$ julia src/bwACprogram.jl -i \"[[3,1],[4,2]]\" -g \"[[1,4],[3,2]]\""
    @add_arg_table! settings begin
        "--optimized", "-o"
        help = "Wheter to use the optimized algorithm"
        arg_type = Bool
        default = false
        "--initial", "-i"
        help = "Initial configuration"
        arg_type = String
        default = "[[1,3,5],[4,2]]"
        "--goal", "-g"
        help = "Goal configuration"
        arg_type = String
        default = "[[2,3,5],[1,4]]"
        "--edge-probability", "-p"
        help = "Edge probability of the underlying random graph"
        arg_type = Float64
        default = 0.1
        "--beta", "-b"
        help = "Plasticity parameter"
        arg_type = Float64
        default = 0.1
        "--area-neurons", "-n"
        help = "Number of neurons per area"
        arg_type = Int64
        default = 1000000
        "--area-winners", "-k"
        help = "Number of winners per area"
        arg_type = Int64
        default = 50
        # consider adding option for `explicit-k` (eak), which is currently set to 50
    end
    return parse_args(settings)
end

function getconfig(s::String)::Vector{Vector{Int}}
    c = try 
        eval(Meta.parse(s)) 
    catch e
        throw(e)
    end
    if !isa(c, Vector{Vector{Int}}) 
        throw(error("Input $(s) cannot be converted to type Vector{Vector{Int}"))
    else
        return c
    end
end

function main()
    args = parse_commandline()
    println("Running AC program for the block world with arguments")
    for (arg,val) in args
        println("  $arg : $val")
    end
    println()
    initial = getconfig(args["initial"])
    goal = getconfig(args["goal"])
    bwACprogram(args["optimized"], initial, goal, args["edge-probability"], args["area-neurons"], args["area-winners"], args["beta"], 50)
end

const ConfigType = Vector{Vector{Int}}
function bwACprogram(optimized::Bool, initial::ConfigType, goal::ConfigType, p::Float64, nean::Int, neak::Int, db::Float64, eak::Int)
        println("Initial configuration: ", initial)
        println("Goal configuration: ", goal)
        println("Edge probability of the underlying random graph: ", p)
        println("Number of neurons per area: ", nean)
        println("Number of winners per area: ", neak)
        println("Plasticity parameter: ", db)
        println("Using optimized algorithm: ", optimized)
        println()
        if optimized
		    planning_experiment(initial, goal, p, eak, nean, neak, db; verbose = true)
        else
            optimized_planning_experiment(initial, goal, p, eak, nean, neak, db; verbose = true)
        end
end

main()