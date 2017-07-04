using Knet

# Defining the MLP
# sample usage: m128 = MLP([3644,128,1])
# mlprun(data; model=m128, epochs=20)

# Define a model type with weights and optimization params so we can
# use the Adam optimizer which works faster than SGD:

type MLP
    weights
    oparams
    function MLP(sizes; optimizer=Adam, winit=0.1, atype=Array{Float32})
        m = new(Any[],Any[])
        for i=2:length(sizes)
            w = convert(atype,winit*randn(sizes[i],sizes[i-1]))
            b = convert(atype,zeros(sizes[i],1))
            push!(m.weights, w)
            push!(m.oparams, optimizer())
            push!(m.weights, b)
            push!(m.oparams, optimizer())
        end
        return m
    end
end

#Prediction function: given weights and input, computes the output
function mlppred(w,x)
    for i=1:2:length(w)-2
        x = max(0, w[i]*x .+ w[i+1])
    end
    return w[end-1]*x .+ w[end]
end

# The logistic loss function
mlploss(w,x,y) = mean(log(1 .+ exp(-y .* mlppred(w,x))))
# Gradient of the loss function
mlpgrad = grad(mlploss)

# Training loop, does one pass over data, modifies mlp in place
function train!(m::MLP, data)
    for (x,y) in data
        dw = mlpgrad(m.weights, x, y)
        for i in 1:length(m.weights)
            update!(m.weights[i], dw[i], m.oparams[i])
        end
    end
end

# Computing logistic loss for model on data
function test(m::MLP, data)
    sumloss = numloss = 0
    for (x,y) in data
        sumloss += mlploss(m.weights, x, y)
        numloss += 1
    end
    sumloss / numloss
end

# Computing classification accuracy for model on data
function acc(m::MLP, data)
    sumloss = numloss = 0
    for (x,y) in data
        z = mlppred(m.weights,x)
        sumloss += mean((z .* y) .> 0)
        numloss += 1
    end
    sumloss / numloss
end

# Function for running training and testing on data with epochs and reporting the results
#   data[1]:train set, data[2]: dev. set
function mlprun(data; epochs=100, sizes=[150,24,1], model=MLP(sizes; atype = typeof(data[1][1][1])))
    #msg(e) = println((e,map(d->acc(model,d),data)...,map(d->test(model,d),data)...)); msg(0)
    println("Accuracies for: train-dev sets:")
    msg(e) = println((e,map(d->acc(model,d),data)...)); msg(0)
    for epoch = 1:epochs
      # Alternative: one could keep the model with highest accuracy in development set results
      # and return that one instead of the last model
        train!(model, data[1])#training on the train set (data[1])
        msg(epoch)
    end
    return model
end

# Adding basic missing functions in Knet: (these are part of Knet 0.8.3 now)
# Base.mean(a::KnetArray) = sum(a)/length(a)
# Base.mean(a::AutoGrad.Rec) = sum(a)/length(a)

# Model conversion for gpu:
function cpu2gpu(m::MLP)
    g = deepcopy(m)
    for i=1:length(g.weights)
        g.weights[i] = KnetArray(g.weights[i])
        g.oparams[i].fstm = KnetArray(g.oparams[i].fstm)
        g.oparams[i].scndm = KnetArray(g.oparams[i].scndm)
    end
    return g
end
