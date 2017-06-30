# Copied from: https://github.com/joefowler/DynamicTimeWarp.jl
# Joe Fowler and Galen O'Neil
# NIST Boulder Laboratories
# December 2014 - February 2015
#
# Copyright (c) 2014: Joseph Fowler.
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

module DynamicTimeWarp

include("WindowedMatrix.jl")


# Dynamic Time Warping with a user-specified distance function

function dtw(seq1::Vector, seq2::Vector, distance::Function=Distance.square)

    # Build the cost matrix
    const m=length(seq2)
    const n=length(seq1)
    cost11 = distance(seq1[1], seq2[1])
    cost = zeros(typeof(cost11), m, n)

    # Initialize first column and first row
    cost[1,1] = cost11
    for r=2:m
        cost[r,1] = cost[r-1,1] + distance(seq1[1], seq2[r])
    end
    for c=2:n
        cost[1,c] = cost[1,c-1] + distance(seq1[c], seq2[1])
    end

    # Complete the cost matrix
    for c=2:n
        for r=2:m
            best_neighbor_cost = min(cost[r-1,c], cost[r-1,c-1], cost[r,c-1])
            cost[r,c] = best_neighbor_cost + distance(seq1[c], seq2[r])
        end
    end

    trackcols, trackrows = trackback(cost)
    cost[end,end], trackcols, trackrows
end


# Compute the optimal track backwards through the cost matrix from end to beginning.
# Return (columns, rows) of the optimal track.

function trackback(D)
    r,c = size(D)
    rows,cols = [r],[c]
    while r > 1 && c > 1
        tb = indmin([D[r-1,c-1], D[r-1,c], D[r,c-1]])
        tb in [1,2] && (r-=1)
        tb in [1,3] && (c-=1)
        push!(rows,r)
        push!(cols,c)
    end
    # Possibly either r>1 or c>1 at this point (but not both).
    # Add the unfinished part of the track to reach [1,1]
    for r=r-1:-1:1
        push!(rows,r)
        push!(cols,1)
    end
    for c=c-1:-1:1
        push!(rows,1)
        push!(cols,c)
    end
    reverse(cols), reverse(rows)
end


# The FastDTW approximation to the DTW, described in "FastDTW: Toward Accurate
# Dynamic Time Warping in Linear Time and Space", S Salvador & P Chan, __Intelligent
# Data Analysis__ (2007).

function fastdtw(seq1::Vector, seq2::Vector, radius::Integer,
                 distance::Function=Distance.square)
    const MinSize = max(radius + 2, 10)
    const N1 = length(seq1)
    const N2 = length(seq2)
    if N1 <= MinSize || N2 <= MinSize
        return (dtw(seq1, seq2, distance))
    end

    # Call recursively on a pair of sequences half this length
    compressed1 = compress(seq1)
    compressed2 = compress(seq2)
    _cost, lowrescol, lowresrow = fastdtw(compressed1, compressed2, radius, distance)

    # Now resample that path to the finer resolution, find the correct
    # window around it, and get the DTW given that window.
    hirescol, hiresrow = expandpath(lowrescol, lowresrow, N1, N2)
    idx2min, idx2max = computewindow(hirescol, hiresrow, radius)
    cost1, newcol, newrow = dtwwindowed(seq1, seq2, idx2min, idx2max, distance)
end



# Given a path through low-res space, generate an approximate path
# through high-res space. It should have dimension Ncol x Nrow

function expandpath(lowrescol, lowresrow, Ncol, Nrow)
    @assert div(Ncol+1,2) == lowrescol[end]
    @assert div(Nrow+1,2) == lowresrow[end]
    const Np = length(lowrescol)
    @assert Np == length(lowresrow)

    hirescol = zeros(eltype(lowrescol), 2*Np)
    hiresrow = zeros(eltype(lowresrow), 2*Np)
    hirescol[1] = hiresrow[1] = c = r = 1
    for i=1:Np-1
        # Select plan according to the next move in lowres path.
        if lowrescol[i+1] == lowrescol[i]  # Next move is up
            r += 1
            hirescol[2*i] = c
            hiresrow[2*i] = r
            r += 1
            hirescol[2*i+1] = c
            hiresrow[2*i+1] = r

        elseif lowresrow[i+1] == lowresrow[i] # Next move is sideways
            c += 1
            hirescol[2*i] = c
            hiresrow[2*i] = r
            c += 1
            hirescol[2*i+1] = c
            hiresrow[2*i+1] = r

        else  # Next move is diagonal.
            c += 1; r += 1
            hirescol[2*i] = c
            hiresrow[2*i] = r
            c += 1; r += 1
            hirescol[2*i+1] = c
            hiresrow[2*i+1] = r
        end
    end
    hirescol[end] = Ncol
    hiresrow[end] = Nrow
    # When expanding to an odd numbered size, it's possible to repeat
    # the last step.  Fix that:
    if hirescol[end]==hirescol[end-1] && hiresrow[end]==hiresrow[end-1]
        hirescol = hirescol[1:end-1]
        hiresrow = hiresrow[1:end-1]
    end
    hirescol, hiresrow
end


# Given the lists of (col,row) indices for the optimal path, compute a "window"
# around that path of the given radius.
# Returns (rowmin, rowmax), each a vector of length pathcols[end], representing
# for each column, the minimum and maximum row numbers used in that column.

function computewindow(pathcols, pathrows, radius)
    const Np = length(pathcols)
    @assert Np == length(pathrows)
    const Ncol = pathcols[end]
    const Nrow = pathrows[end]

    # Find the min/max row at each column in the path.
    pathmin = zeros(Int, Ncol)
    pathmax = zeros(Int, Ncol)
    for i=1:Np
        c,r = pathcols[i], pathrows[i]
        pathmax[c] = r
        if pathmin[c] == 0
            pathmin[c] = r
        end
    end

    # The window in each column for "radius" r starts at the pathmin
    # of the rth-previous column and ends at the pathmax of the
    # rth-next column, plus (in each case) the radius.
    if radius < Ncol-1 && radius < Nrow-1
        rowmin = vcat(fill(1,radius), pathmin[1:end-radius]-radius)
        rowmax = vcat(pathmax[radius+1:end]+radius, fill(Nrow,radius))

        # Window values must be in the range [1:Nrow].
        for c=1:Ncol
            if rowmin[c]<1; rowmin[c]=1; end
            if rowmax[c]>Nrow; rowmax[c]=Nrow; end
        end
    else
        rowmin = fill(1,Ncol)
        rowmax = fill(Nrow,Ncol)
    end
    rowmin, rowmax
end



# Do DTW in a subset of the full space, the subset specified by
# a "window". Arguments [idx2min,idx2max] give the inclusive lowest
# and highest index in the seq2 direction, one element for each index along
# the seq1 direction. Thus seq1, idx2min, and idx2max should all be of
# equal length

function dtwwindowed(seq1::Vector, seq2::Vector,
                     idx2min::Vector, idx2max::Vector,
                     distance::Function=Distance.square)

    const m=length(seq2) # of rows  in cost matrix
    const n=length(seq1) # of columns in cost matrix
    @assert n==length(idx2min)
    @assert n==length(idx2max)
    @assert 1==minimum(idx2min)
    @assert m==maximum(idx2max)

    # Build the (n x m) cost matrix into a WindowedMatrix, because it's ragged.
    # That type gives efficient storage with convenient [r,c] indexing and returns
    # Inf when accessed outside the window.
    cost = WindowedMatrix(idx2min, idx2max, Inf)

    # First column first
    cost[1,1] = distance(seq1[1], seq2[1])
    for r=2:idx2max[1]
        cost[r,1] = cost[r-1,1]  + distance(seq1[1], seq2[r])
    end

    # Complete the cost matrix from columns 2 to m.
    for c=2:n
        for r=idx2min[c]:idx2max[c]
            best_neighbor_cost = min(cost[r-1,c], cost[r-1,c-1], cost[r,c-1])
            cost[r,c] = best_neighbor_cost + distance(seq1[c], seq2[r])
        end
    end
    trackcols, trackrows = trackback(cost)
    cost[end,end], trackcols, trackrows
end


function compress(seq::Vector)
    Navg = div(length(seq), 2)
    evenseq = 0.5*(seq[1:2:end-1]+seq[2:2:end])
    if length(seq)%2 == 1
        return vcat(evenseq, [seq[end]])
    end
    evenseq
end




using PyPlot

function dtwspeedtest(N::Integer, functype=1)
    @assert N>10
    if functype == 1
        stretch = div(N,10)
        stretch = max(stretch, 10)
        t = linspace(0,1,N-stretch)
        x = vcat(fill(0.0, stretch), sin((t.^1.3).*5))
        y = vcat(sin((t.^0.7).*5), fill(sin(5), stretch))
    else
        # Random path
        rows=Array(Int, N)
        cols=Array(Int, N)
        r=c=1
        const PSTRETCH=0.8
        for i=1:N
            rnum = rand()
            if rnum > PSTRETCH # go diag
                r+=1; c+=1
            elseif rnum > 0.5*PSTRETCH
                r+=1
            else
                c+=1
            end
            rows[i] = r
            cols[i] = c
        end
        t = linspace(0, N, max(r[end],c[end]))
        f = sin((N-t).^0.7)+ (t.*(3./N))
        x = f[rows]
        y = f[cols]
    end

    clf()
    plt.subplot(211)
    dstd = std(x)+std(y)
    plot(x-mean(x), "r")
    plot(y-mean(y)+dstd, "b")

    plt.subplot(224)
    radii=[40,20,10,5,2]
    colors=["red","orange","gold","green","cyan"]
    for (rad,color) in zip(radii, colors)
        print("Testing size $N radius $rad: ")
        @time _,r1,c1=fastdtw(x,y,rad)
        plot(r1,c1,color=color)
    end

    if N<=2500
        @time _,r,c=dtw(x,y)
        plot(r,c,"k")
    end
    nothing
end


function plotdtw(seq1::Vector, seq2::Vector, offset=0.0)
    cost, match1, match2 = dtw(seq1, seq2)
    if offset == 0.0
        offset = 2*(std(seq1) + std(seq2)) + mean(seq1) - mean(seq2)
    end
    clf()
    plot(seq1, "-r", seq2+offset, "-b")
    for i=1:length(match1)
        plot([match1[i],match2[i]]-1, [seq1[match1[i]], seq2[match2[i]]+offset],
             color="gray")
    end
end



#function dtwbaryavg_iteration(dbavg::Vector, sequences::Array{Array{1},1})
function dtwbaryavg_iteration(dbavg::Vector, sequences::Array)
    const Nseq = length(sequences)
    count = zeros(Int, length(dbavg))
    sumcoords = zeros(Float64, length(dbavg))

    for i=1:Nseq
        cost, match1, match2 = dtw(dbavg, sequences[i])
        for j=1:length(match2)
            count[match1[j]] += 1
            sumcoords[match1[j]] += sequences[i][match2[j]]
        end
        println("Compared $i to the standard")
    end
    sumcoords = sumcoords ./ count
end



#function dtwbaryavg{T<:Real}(sequences::Array{Array{T,N},1})
function dtwbaryavg(sequences)
    dbavg = copy(sequences[1])
    clf()
    for i = 1:5
        dbavg = dtwbaryavg_iteration(dbavg, sequences)
        plot(dbavg)
    end
    dbavg
end


# This sub-module contains various functions returning pointwise distances
# between elements from two sequences.

module Distance
square(x,y) = (x-y)^2
absval(x,y) = abs(x-y)

# And now a Poisson-sensitive distance.
# Call this with the 2 sequences, and it returns a closure, a function that
# measures Poisson-distance for elements of those two particular sequences.
# Alternately, call it with two integers, representing the number of elements
# in each sequence.
function poissonpenalty(seq1::Vector, seq2::Vector)
    n1=sum(seq1)
    n2=sum(seq2)
    poissonpenalty(n1, n2)
end


function poissonpenalty(n1::Integer, n2::Integer)
    function pdist(x,y)
        if x<=0 && y<=0; return 0; end
        mu = float(x+y)/float(n1+n2)
        return ((x-mu*n1)^2 + (y-mu*n2)^2) / (x+y)
    end
end

end # module Distance


export
dtw,
dtwwindowed,
fastdtw,
dtwbaryavg

end # module
