struct MatrixInterpolation
    size::Tuple{Int64, Int64}
    ys
    xs
end

# struct for matrix interpolations
# to be consistent with DataInterpolations ys are given before xs
"""
    LinearMatrixInterpolation(ys, xs; extrapolate=false)

## Arguments

    - `ys`: data points; each element is a matrix.
    - `xs`: time points.

## Keyword Arguments

    - `extrapolate`: boolean value to allow extrapolation. Defaults to `false`.
"""
function MatrixInterpolation(
    ys::Vector{Matrix{AbstractFloat}}
    xs::Vector{AbstractFloat},
)
    Nbands, _ = size(first(ys))
    ys_tensor = reshape(reduce(hcat, ys), Nbands, Nbands, :)

end

(interp::MatrixInterpolation)(t::Number) = _interpolate(interp, t)

function _interpolate(interp, t)
    
end

function interpolatematrix(ρs::Vector{Matrix{AbstractFloat}})
    
end