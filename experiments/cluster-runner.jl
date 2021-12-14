include("experiment_functions.jl")
Random.seed!(0)


getarg(s::String) = eval(Meta.parse(s))
optimized = ARGS[1]
input_stacks = getarg(ARGS[2])
goal_stacks = getarg(ARGS[3])
# example
# julia experiments/runner.jl 3 "[[1,3],[2,4]]"

nean = 4*10^6
println(ARGS)
flush(stdout)

if (optimized == "optimized")
    s = optimized_planning_experiment(input_stacks, goal_stacks,  0.1, 50, nean, 50, 0.1; to_readout = false, verbose = true)
else
    s = planning_experiment(input_stacks, goal_stacks,  0.1, 50, nean, 50, 0.1; to_readout = false, verbose = true)
end
println(s)
flush(stdout)

