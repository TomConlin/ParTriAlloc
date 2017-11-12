ignore = """
given N side of (lower) triangle matrix

k              is number of processors

N^2 / 2        is the apx area of work to be done

N^2 / 2k       is each processors fairish share

a:x            is an interval in 1:N

a * (x-a) + (x-a)^2 / 2  is the area of work in the interval

Given:
    a = 1
    N = 100
    k = 4


N^2 / 2k = 10000 / 8  == 1250  units of work
(expect [12 +/- 1] as first/last interval)


N^2 / 2k  = a * (x-a) + (x-a)^2 / 2

N^2 / k  = 2a * (x-a) + (x-a)^2

N^2 / k  = 2ax - 2a^2 + (x-a)^2

N^2 / k + 2a^2  = 2ax + (x-a)^2

yuck, find a solver

x = (k*a^2 + N^2)^(1/2) / k^(1/2)


N=60000; k=4; a=1;
for z=1:k
    a = round((k*a^2 + N^2)^(1/2) / k^(1/2))
    println(a)
end

30000.0
42426.0
51961.0
60000.0

for a lower triangle, reverse and subtract from N

"""



"""
Given:
 - `N` the size of a triangular matrix (one dimension)
 - `k` the number of processes sharing the work
 - `lower`  [true]|false

    Defaults to lower triangular matrix.  
    Include false for upper triangular matrix  

    expects:  1 < k < N

Return Array of k intervals in 1:N balancing area covered per interval

"""
function triallocate(N::Int64, k::Int64 , lower::Bool=true)
    if N < k || k < 2
        return [1,N]
    end
    result = Array{Real}(2k)
    denom =  sqrt(k)
    N2 = N*N
    x = 1
    for i = 1:k-1
        w = x>1 ? x + 1 : 1
        x = (k*x^2 + N2)^(1/2) / denom
        result[2i-1] = w
        result[2i] = x
    end
    result[2k-1] = x+1
    result[2k] = N
    if lower == true
        reverse!(result)
        result =  N - result
        result[1] = 1
        result[2k] = N
    end
    # round down after chosing direction
    result = [Int(floor(result[i])) for i=1:2k]
    return [result[2i-1]:result[2i] for i=1:k]
end
