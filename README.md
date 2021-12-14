# jBrain
Software associated with the AAAI 2022 paper

*Francesco D'Amore, Daniel Mitropolsky, Pierluigi Crescenzi, Emanuele Natale, Christos H. Papadimitriou*, **Planning with Biological Neurons and Synapses**.

In order to get a first overview of the software, you are encouraged to issue the terminal command 

`julia src/bwACprogram.jl`

from the project folder, which will run a simple default experiment. The input and output configurations, and the main model parameters, can be chosen by providing suitable options, which can be displayed by adding the `--help` flag in the aforementioned command.

Some basic program simulations can be run by compiling the file `src/experiments.jl` and executing, for example, the following functions:

`test_parse_and_read([[1,3,5,7],[2,4,6]],0.1,50,1000000,50,0.1)`

or

`test_optimized_planning([[1,3,5,7,9],[2,4,6,8]],[[2,4,6,8,1,3,5,7,9]],0.1,50,1000000,50,0.1)`.
