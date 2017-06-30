#Loading Dynamic Time Warping function
# Credits: Joe Fowler and Galen O'Neil, NIST Boulder Laboratories
# link: https://github.com/joefowler/DynamicTimeWarp.jl
include("DynamicTimeWarp.jl");
using StatsBase, DynamicTimeWarp

#Feature vector computation: vect=[histogram(delta_signal_after_dtw), dtw_cost, dtw_length_change]
# DTW applied to two signals
# The difference between two matched signals are found
# The histogram of the difference is returned as the feature vector
# We would expect the highest errors would be due to octave errors in pitch estimation
# If directly used in training the MLP, these will possibly get lower weights
# So, the last two elements are replaced with the cost computed from DTW matching and
# and the ratio of time-extention applied to match the signals
# which we think could be more discriminant features
function DTWCostHist(sig1,sig2,numHistBins)
  costScaleFactor=30.0;#manually set scaling/normalisation factor for the dtw cost
  #octave wrapping before operations
  sig1=mod(vec(sig1),1.0);sig2=mod(vec(sig2),1.0);
  #dtw matching of the two signals
  cost, match1, match2 = DynamicTimeWarp.fastdtw(sig1,sig2,16);
  #dist: absolute delta signal after dtw matching of two signals
  dist=abs(sig1[match1]-sig2[match2]);
  
  h = fit(Histogram, dist, 0:1/numHistBins:1.0);
  hist=h.weights/length(match1);
  #last two values altered to contain two other features
  hist[end]=cost/costScaleFactor;
  hist[end-1]=abs(length(match1)-length(sig1))*2/length(match1);
  return hist;
end

# Generate positive and negative instances(feature vector for pairs of melodic segments):
# Returns a pair of matrices: one for true one for false examples.
# The columns of matrices contain the feature vectors.
function pairdata(data,numHistBins)
  tout = Float32[];#tout: true-out: matching melodic pairs' (piano-singing) feature
  fout = Float32[];#fout: false-out: non-matching melodic pairs' (piano-singing) feature
  for d in data
    refT = d["RefSegsTrue"];#in julia string indexes can be safely used :)
    perT = d["PerSegsTrue"];
    perF = d["PerSegsFalse"];

    #Pair all reference segments with performance segments marked as true/pass
    # compute the feature vector and the add to true-pairs-data
    for i=1:length(refT)#all reference recordings
      for j=1:length(perT)#all true-performance recordings
        # Computation of the feature as histogram of differences after matching two signals by DTW
        # the vector is appended to true-output pair
        # If you would prefer another input for the MLP, consider modifying this line
        append!(tout, DTWCostHist(refT[i],perT[j],numHistBins));
      end
    end

    #pair all reference segments with performance segments marked as false/fail
    # compute the feature vector and the add to false-pairs-data
    for i=1:length(refT)#all reference recordings
      for j=1:length(perF)#all false-performance recordings
        #computation of the feature as histogram of differences after matching two signals by DTW
        # the vector is appended to false-output pair
        # If you would prefer another input for the MLP, consider modifying this line
        append!(fout, DTWCostHist(refT[i],perF[j],numHistBins));
      end
    end
  end

  #re-shape feature vectors computed from true and false pairs in matrix form and return both
  rows=numHistBins;
  cols = div(length(tout),rows);
  tout = reshape(tout, (rows,cols));
  cols = div(length(fout),rows);
  fout = reshape(fout, (rows,cols));
  return (tout, fout);

end


# Split the pair data into a balanced train, dev and test set with a given minibatch size.
#   pdata: tout and fout returned by pairdata(),
#   contains feature vectors computed from true and false pairs of melodic segments
#
#   While this function can be used to split a whole pool into three packages, we prefer to
#   handle this operation in two steps because we want to use different pools for train,
#   development and test sets. So, this function is called twice: one for
#   preparing the training and development sets
#   and one for the test set

function trntst(pdata; batch=100, splt=((50000,50000),(5000,5000),(5000,5000)))
  (tout,fout) = pdata;
  (nd,nt) = size(tout);#nt: number of true-pairs
  (nd,nf) = size(fout);#nf: number of false-pairs
  rt = randperm(nt);#forming random order of indexes to shuffle
  rf = randperm(nf);
  nt = nf = 0
  bdata = Any[];#batch data for all sets of splits
  for (t,f) in splt
    if(t>0 && f>0)#if there exists (0,0) in splits, skip it
      #forming input-vector for MLP
      # Left half contains true samples, right half contains false samples
      x = hcat(tout[:,rt[nt+1:nt+t]], fout[:,rf[nf+1:nf+f]]);
      #forming the output vector: 1:true, -1: false
      y = hcat(ones(Float32,1,t), -ones(Float32,1,f))

      nt += t; nf += f;
      r = randperm(t+f);
      batches = Any[];
      #forming a single batch in each loop, putting in 'batches'
      for i=1:batch:length(r)-batch
        xbatch = x[:,r[i:i+batch-1]];
        ybatch = y[:,r[i:i+batch-1]];
        push!(batches, (xbatch, ybatch))
      end
      #adding to 'batches' for this split in 'bdata'
      push!(bdata, batches);
    end
  end
  return bdata
end

# To load and save KnetArrays:
import JLD: writeas, readas
type _KnetArray; a::Array; end
writeas(c::KnetArray) = _KnetArray(Array(c))
readas(d::_KnetArray) = KnetArray(d.a)
