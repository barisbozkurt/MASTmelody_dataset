#function to gather all melodic segments (f0-sequences) by
# - reading all f0 files
# - converting hz to cents
# - grouping them with respect to melody (40 melodies exist in the dba)
# - resampling f0-series to the mean length of all
# - saving the groups to json files
# - normalizing f0-cent information in range fmin=5800, fmax=7300
# - saving normalized segments to a jld file

using JSON, JLD, DSP;

#melodyPool carries all performances for a melody
type melodyPool
  melodyIndex;#index for the pool in the whole dba
  RefSegsTrue;#data for piano recordings which serve as reference
  PerSegsTrue;#data for signing(performance) recordings labeled as pass/true
  PerSegsFalse;#data for signing(performance) recordings labeled as fail/false
end

#helper method to find melody index in pools
#used for assigning a file's data with the corresponding melody pool
function findIndexWithMelInd(pools,melInd)
  for k=1:length(pools)
    if pools[k].melodyIndex==melInd
      return k;#if found return it
    end
  end
  #if non-existant, pick one pool with index=-1
  for k=1:length(pools)
    if pools[k].melodyIndex==-1
      return k;
    end
  end
end

#Function to remove zero-valued elements, perform ocatve correction and resampling
function removeZerosResample(arr,targetLen)
  #to avoid re-sampling artifacts at boundaries,
  #we extend length by 2*offset points and further removed from the ends of the series
  offset=10;
  isPositive(x)=(x>0);#in-line function definition to be used in finding positive valued elements
  #modifying the array
  arr=arr[find(isPositive,arr)];#removing all zero valued elements
  #correcting octave jumps by finding jumps in the derivative and correcting them
  diffArr=diff(arr);
  for k=1:length(diffArr)
    while diffArr[k]>800
      diffArr[k]=diffArr[k]-1200;
    end
    while diffArr[k]<(-800)
      diffArr[k]=diffArr[k]+1200;
    end
  end
  arr=cumsum(diffArr)+arr[1];
  arr=resample(arr,convert(Float64,targetLen+offset*2)/length(arr));#resampling, we avoid integer division by convert
  return arr[offset:offset+targetLen-1];#removing ends to get rid of re-sampling artifacts
end

#Reshaping all f0-series in the pool:
# All zero values are removed from f0-series data (you may like to question and modify that :) )
# octave-correction within an f0-series applied
# f0-series resampled to have 'targetLen'
# further a final step of octave correction is applied that compares performance recordings with reference recordings
#   and shifts points at a distance of 800 cents larger than median of references by an octave
function reshapeMelSeg(pool,targetLen)
  meansRefSeg=[];
  for k=1:length(pool.RefSegsTrue)
    pool.RefSegsTrue[k]=removeZerosResample(pool.RefSegsTrue[k],targetLen);
    push!(meansRefSeg,mean(pool.RefSegsTrue[k]));
  end
  medianOfMeansRef=median(meansRefSeg);#this value will be used to apply to decide octave shift for performance segments
  for k=1:length(pool.PerSegsTrue)
    pool.PerSegsTrue[k]=removeZerosResample(pool.PerSegsTrue[k],targetLen);
    #if performance is at an octave difference compared to references, shift it
    meanPer=mean(pool.PerSegsTrue[k]);
    if (meanPer-medianOfMeansRef)>800
      pool.PerSegsTrue[k]=pool.PerSegsTrue[k]-1200;
    elseif (meanPer-medianOfMeansRef)<(-800)
      pool.PerSegsTrue[k]=pool.PerSegsTrue[k]+1200;
    end
  end
  for k=1:length(pool.PerSegsFalse)
    pool.PerSegsFalse[k]=removeZerosResample(pool.PerSegsFalse[k],targetLen);
    #if performance is at an octave difference compared to references, shift it
    meanPer=mean(pool.PerSegsFalse[k]);
    if (meanPer-medianOfMeansRef)>800
      pool.PerSegsFalse[k]=pool.PerSegsFalse[k]-1200;
    elseif (meanPer-medianOfMeansRef)<(-800)
      pool.PerSegsFalse[k]=pool.PerSegsFalse[k]+1200;
    end
  end
  return pool;
end

#Function to pool f0s to json packages for each melody (a total of 40 melodies exist in the dba)
# - reading all f0 files
# - converting hz to cents
# - grouping them with respect to melody (using melodyIndex deduced from filename)
# - resampling all f0-series to have the same length (f0SeriesLen=1000)
# - saving the groups(each containing reference and performance recordings
#   of a specific melody) to json files
function poolF0s_SaveInJson(dbaDir)
  #variable initialization
  c1_hz=float(8.17579891564371);#frequency of C1 for conversion from hz to cents
  f0SeriesLen=1000;#manually set value for number of f0 points in a melodic segment
  #creating and initializing melody pools (array of 'melodyPool' type defined above)
  pools=melodyPool[];
  for ind=1:40
    #melodyIndex assigned to -1 indicates that melodyIndex has not been assigned yet
    push!(pools,melodyPool(-1,Any[],Any[],Any[]));
  end
  arrayLengths=Any[];
  #assuming all f0 files are gathered in the 'dbaDir' folder, reading all
  files = readdir(dbaDir);
  numFiles=0

  for f in files#for each f0 file ('f0s.txt')
    if contains(f,".f0s.txt")
      numFiles=numFiles+1
      print(numFiles," ")
      f0cent=Float64[];#initialize f0 array for the current file
      #read f0 file
      fid = open(joinpath(dbaDir,f));
      for ln in eachline(fid)
        vals=split(ln);f0hz=float(vals[2]);#reading second column as f0-hz info
        if(f0hz>30)#frequencies below 30Hz will be considered as zero cent
          #push!(f0cent,convert(Int32,round(log(2,float(vals[2])/c1_hz)*1200)));#conversion to cents as Int32
          push!(f0cent,round(log(2,float(vals[2])/c1_hz)*1200));#conversion to cents
        else
          push!(f0cent,0.0);
        end
      end
      close(fid);
      push!(arrayLengths,length(f0cent));

      #putting f0-cent time series as melody segment to its corresponding pool
      # file name starts with melody index and mel1/mel2
      # melInd is formed by concatenating digits before the first two '_'s
      melInd=parse(Int,string(f[1:search(f,'_')-1],f[search(f[5:end],'_')+3]));
      # find index of the pool with that melodyIndex
      #   if not already exists, it will be assigned to an empty one (with -1)
      ind=findIndexWithMelInd(pools,melInd);
      pools[ind].melodyIndex = melInd;
      #using filename to decide type of the melodic segment
      if contains(f,"ref")#reference melody
        push!(pools[ind].RefSegsTrue,f0cent);
      else#student performance melody
        if contains(f,"pass")#performance-true/pass
          push!(pools[ind].PerSegsTrue,f0cent);
        else#performance-false/fail
          push!(pools[ind].PerSegsFalse,f0cent);
        end
      end
      #melody segment added to the pool, move to next f0s.txt file
    end
  end
  println(" data files read");
  #Re-shaping f0series data and writing to json files;
  println("Writing grouped melodic segments to json files")
  for p in pools
    if p.melodyIndex!=-1
      #removing zeros from f0-series, resampling, octave correction
      p=reshapeMelSeg(p,f0SeriesLen);
      #writing pool to a json file
      open(joinpath(dbaDir,string(p.melodyIndex,".json")), "w") do f
        write(f, JSON.json(p));
        print(string(p.melodyIndex,".json"),"\t")
      end
    end
  end
  println();
  return pools;
end

#loading and normalizing data in json files from the dba folder
function loaddata(dbaDir, fminCent=4800, octaveSize=1200)
  println("Loading json files to start training")
  data  = Any[]
  files = readdir(dbaDir)
  for f in files
    ismatch(r"\.json$", f) || continue
    dict = JSON.parsefile(joinpath(dbaDir,f))
    for (k,v) in dict
      if isa(v,Void) || isa(v,Int64)#skipping melodyIndex or empty arrays
        dict[k] = Any[]
      else
        for i=1:length(v)
          #normalisation
          v[i] = (convert(Array{Float32},v[i]) - fminCent) ./ octaveSize;
        end
      end
    end
    #saving each pool of melody with its reference and performance segments
    # as an element in data. Ex: to access data of the first melody, use data[1]
    push!(data, dict)
  end
  return data
end

# Running basic steps of the database preparation process
function runDBAPrepProcess(dbaDir)
  #Reading f0s.txt files, grouping with respect to melody index, saving to Json files
  println("Reading f0 files")
  poolF0s_SaveInJson(dbaDir);
  #loading and normalizing data in [0,2] range(2 octaves) and saving as a jld archive file
  JLD.save(joinpath(dbaDir,"groupedMelSegData.jld"),"data",loaddata(dbaDir));
end
