function consecutive(f, A::AbstractVector)
    [ f(A[i+1], A[i]) for i = 1:length(A)-1 ]
end

"""
    Cumulative integral of a function given as tabulate values in arrays xs and ys.
"""
function cumintegrate(xs, ys)
    # trapezoid rule
    xdif = diff(xs)
    yavg = consecutive((x, y) -> (x + y)/2, ys)
    cumint = cumsum(xdif .* yavg)
    # add zero at the start
    pushfirst!(cumint, 0.)
    return cumint
end

function cumintegrate(xs, ys, lower, upper)
    # assume upper < xs[end] and lower < xs[1]
    upper < xs[end] || error("upper bound in cumintegrate is outside the tabulated interval")
    lower < xs[1] || error("lower bound in cumintegrate should be outside the tabulated interval")
    # the data should be sorted by xs
    p = sortperm(xs)
    xs, ys = xs[p], ys[p]
    len = length(xs)
    sum = 0.
    cumint = zeros(len)
    for i in 2:len
        x0, y0 = xs[i-1], ys[i-1]
        x1, y1 = xs[i], ys[i]
        yavg = (y1 + y0) / 2
        dx = 0.
        if x0 < upper & x1 ≤ upper
            dx = x1 - x0
        elseif x0 < upper & x1 > upper
            dx = upper - x0
        else
            dx = 0.
        end
        sum += dx * yavg
        cumint[i] = sum
    end
    return cumint
end
