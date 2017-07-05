#This script involves basic steps involved in reading-grouping data
# running tests for training a very simple MLP model for
# automatic singing assessment. The model is implemented using the Knet toolbox.
# The readers are invited to develop more complicated models using the example
# provided here which should be straight forward.
#
# For referring to this implementation in your work please use:
#   Bozkurt, B., Baysal, O., Yuret, D. A Dataset and Baseline System for Singing
#   Voice Assessment, 13th Int. Symposium on Computer Music Multidisciplinary
#   Research, Porto, Sept. 25-28, 2017.
#
# Dependencies: KNET, JSON, JLD, DSP,
# 'Dynamic Time Warping' by Joe Fowler and Galen O'Neil, NIST Boulder Laboratories
# https://github.com/joefowler/DynamicTimeWarp.jl

#Installing dependencies if not installed
for p in ("JSON","JLD","DSP","Knet")
    if Pkg.installed(p) == nothing; Pkg.add(p); end
end

# Setting database and Julia tools directories
# 'dbaDir' is the folder that contains all f0s.txt files
# 'toolsDir' is the folder that contains all julia scripts in this example
toolsDir=dirname(@__FILE__); #current directory assigned as the tools directory
dbaDir=replace(toolsDir,"baseline","f0data");

#Database processing from f0s.txt files:
#The following lines will load the functions in "gatherMelSegsInAPool.jl" and
# call the 'runDBAPrepProcess' function to read the txt files, pool them
# with respect to melodies, produce json files containing these melody sets
# and a jld file("groupedMelSegData.jld") contaning all data
include("gatherMelSegsInAPool.jl");
runDBAPrepProcess(dbaDir);

#Running automatic scoring tests:
# Each instance will be created by pairing two recordings from the same melody pool:
# a reference recording and a performance recording

#Loading data to 'data' (groupedMelSegData.jld file created by the previous line)
data = load(joinpath(dbaDir,"groupedMelSegData.jld"), "data")

#Running automatic grading learning and testing steps
#
# The model is a very simple linear model which takes in a histogram of
# differences computed by subtracting two melodic curves and has a binary output [pass or fail]
# The histogram is computed on a single octave, so the bin size is=1200cent/inputVectorSize
# One can easily test much more complicated models by altering the MLPsizes, the variable
# which defines number of neurons in layers of an MLP (see 'mlp.jl')
# and incorporating other features about the distance between two f0-series data
#
include("MASTmelSimTest.jl");#includes functions for data splitting, running learning tests
inputVectorSize=150
MLPsizes=[inputVectorSize,1]
runTests(data,MLPsizes);
