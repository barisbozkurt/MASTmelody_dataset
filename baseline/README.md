[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-ff69b4.svg)](http://creativecommons.org/licenses/by-nc-sa/4.0/)

# Baseline system for singing assessment from f0-series data

A baseline system for assessment (comparing a singing sample to a reference piano recording):
> Bozkurt, B., Baysal, O., Yuret, D. A Dataset and Baseline System for Singing Voice Assessment, 13th Int. Symposium on Computer Music Multidisciplinary Research, Porto, Sept. 25-28, 2017.
```
@inproceedings{inproceedings,
  author={Bozkurt, B., Baysal, O., Yuret, D.},
  title={A Dataset and Baseline System for Singing Voice Assessment},
  year={2017},
  booktitle={13th Int. Symposium on Computer Music Multidisciplinary Research, CMMR 2017}
}
```

Please cite the publication if you use this baseline system in your work.

For experimenting with the tools, the following database is required: <a href="https://github.com/barisbozkurt/MASTmelody_dataset#mastmelody_dataset">The dataset of f0-series data computed from signing and piano samples</a>. This dataset can be found under the f0data directory of the <a href="https://github.com/barisbozkurt/MASTmelody_dataset">MASTmelody_dataset repo</a>.

<a name="Introduction"></a>Introduction
--------------------
This repository contains a singing assessment system implemented in <a href="https://julialang.org">Julia language</a> using <a href="https://github.com/denizyuret/Knet.jl">the Knet toolbox</a>. Please use the version v0.5.2 of the Julia language, available at the <a href="https://julialang.org/downloads/oldreleases.html">Julia downloads</a> website.

The codes involve basic steps for reading and grouping data, running tests for training a very simple MLP model for automatic singing assessment. The model is implemented using the Knet toolbox. The readers are invited to develop more complicated models using the example provided here which should be straight forward.

Dependencies (Julia packages): <a href="https://github.com/denizyuret/Knet.jl">Knet</a>, <a href="https://github.com/JuliaIO/JSON.jl">JSON</a>, <a href="https://github.com/JuliaIO/JLD.jl">JLD</a>, <a href="https://github.com/JuliaDSP/DSP.jl">DSP</a> and <a href="https://github.com/joefowler/DynamicTimeWarp.jl">Dynamic Time Warping</a> by Joe Fowler and Galen O'Neil, NIST Boulder Laboratories.

In Julia, these packages can be easily installed using a line of command: Pkg.add("PackageName") [ex: Pkg.add("Knet")] except the last one. For this reason "DynamicTimeWarp.jl" and "WindowedMatrix.jl" files were included in this repo.

#### Where to start:
First, you would need to <a href="https://github.com/barisbozkurt/MASTmelody_dataset">clone or download the repository</a> that contains both the data and this actual codes repository. The tools are designed to run with two lines of code if you do not change the folder structure of the repo.

MASTbaselineProcess.jl is the batch file that contains the sequence of operations for running ONE learning experiment. It can be launched in Julia via first changing directory to your local directory where MASTbaselineProcess.jl exists:
```
julia> cd("MASTmelody_dataset/baseline")
```

... and running the batch:
```
julia> include("MASTbaselineProcess.jl")
```

Here are the steps of the process executed with this batch:
1) Installing dependencies,
2) Running data reading (of 3618 f0 files), grouping process (calling the 'runDBAPrepProcess' function implemented in "gatherMelSegsInAPool.jl") which will produce json files containing these data grouped in melody sets and a jld file ("groupedMelSegData.jld") contaning all data.
3) The data will be re-loaded from the .jld file (which in fact does not need to be re-created in each experiment. Step 2 can be omitted to repeat the experiment using the same data)
4) Defining MLP sizes and running the learning experiments
5) Saving the final model

These processes will create the following new files in your "f0data" directory:
1) "trainedModel.jld": contains the trainedModel
2) Json files for each melody set: f0-series data that is grouped with respect to melody (i.e. f0-series for all performances of one single melody) "RefSegsTrue": f0 data for the reference recordings, "PerSegsTrue": f0 data for singing performances labeled as true/pass, "PerSegsFalse": f0 data for singing performances labeled as false/fail.
3) groupedMelSegData.jld: all data packed in one single jld file that can be loaded with a single command: `data = load("groupedMelSegData.jld","data")`.  

For further details involved please refer to the comments in the codes and the scientific paper referred above.

<a name="License"></a>License
--------------------
<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

<a name="Acknowledgement"></a>Acknowledgement
--------------------
This dataset has been curated within the TUBITAK (The Scientific and Technological Research Council of Turkey) funded research project 1001-215K017 targeting development of automatic assessment tools for music performances.
