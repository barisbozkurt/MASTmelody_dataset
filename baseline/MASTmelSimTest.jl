using Knet, JLD;

#loading data processing and pooling functions
include("formdata.jl");
#loading multi-layer perceptron functions
include("mlp.jl")

#Function to decide split values starting from data size
# option=1->create splits/sizes for train and development data
# option=2->create splits/sizes for test data and a dummy (0,0)
#   as the minimum of the lowest size of true-pairs data and false pairs-data
#   to guarantee a balanced set of data for the tests.
#   i.e. the splits will contain equal number of true and false pairs
function decideSplits(pairData,option)
  if option==1
    (tempX,tempY)=pairData;
    minNumPairs=min(size(tempX,2),size(tempY,2));
    splitRatio=0.9;#split for train and development
    splts=((Int(floor(minNumPairs*splitRatio)),Int(floor(minNumPairs*splitRatio))),(Int(floor(minNumPairs*(1-splitRatio))),Int(floor(minNumPairs*(1-splitRatio)))));
    return splts;
  elseif option==2
    (tempX,tempY)=pairData;
    #minimum of number of true-pairs and false-pairs will be used for training
    minNumPairs=min(size(tempX,2),size(tempY,2));
    splts=((minNumPairs,minNumPairs),(0,0));
    return splts;
  end
end

#Function to run all tests for training
# For seperating the test set, data is splitted in melody sets so that
# samples of the same melody do not fall into both train and test sets
#
# runTests(data,sizes):
# sizes=[numHistBins(input vector size),...,1(output node)] (sizes of MLP)
function runTests(data,sizes)
  print("MLP sizes: ");println(sizes);
  numHistBins=sizes[1];#feature_vector/input size, first layer of MLP

  #Random separation for the test/train sets (%20, %80)
  # !!Since this separation is random, the results will vary each time
  # you run your tests. We left the decision of the test-train separation to
  # the user. In our cross-validation tests (performing this in a for loop and
  # averaging the accuracies obtained for the test set we observed an average
  # accuracy of 0.74)

  randInds=randperm(length(data));
  #preparing train data from 32 randomly selected melody pools
  dataTrain=data[randInds[1:32]];

  #Reading melodic segments, producing feature vector for
  # each piano-singing pair in a melodic set. The output is a tuple that carries
  # features of true-pairs(grade=pass) and features of false-pairs(grade=fail)
  println("Gathering features for training data...involves DTW computation, will take some time...(maybe 10 mins.)")
  trainPairData=pairdata(dataTrain,numHistBins);#see loaddata.jl for implementation

  # deciding how the trainPairData should be splitted
  # to have balanced set of true and false samples
  splts=decideSplits(trainPairData,1);
  #creating batch data for the train and development sets
  bdata = trntst(trainPairData,batch=100,splt=splts);
  print("Splits for train-dev sets (number of true-pairs and false-pairs): ");println(splts);

  #run training
  model=mlprun(bdata; epochs=100, sizes=sizes);

  #preparing test set from the other 8 melody pools
  println("Gathering features for test data...involves DTW computation, will also take some time...")
  dataTest=data[randInds[33:40]];
  testPairData=pairdata(dataTest,numHistBins);
  splts=decideSplits(testPairData,2);println("Number of true-pairs and false-pairs: ",splts[1]);
  #creating single batch data for the test set
  bTestData=trntst(testPairData,batch=100,splt=splts);
  println("Accuracy for the test set:")
  println(map(d->acc(model,d),bTestData));
  #save the model in a jld file (remember, the test-train separation was random)
  println("Saving the model to trainedModel.jld")
  save("trainedModel.jld", "model", model);

  println("FINISHED!")
end
