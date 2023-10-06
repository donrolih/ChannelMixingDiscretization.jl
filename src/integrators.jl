"""
    Given an array of K points ωs = [ω_1, ..., ω_K] and the values of the function f in an array values at these points, calculate the integral 
    ```math
    F(\omega) = \int_{\omega_1}^{\omega} f(\omega') \mathrm{d} \omega'
    ```
    using Gaussian quadrature from the QuadGK.jl package. If you know there
    is a singularity inside the interval (ω_1, ω_K), you should provide an
    array of singularities as an optional parameter.
"""
function integrate(ωs, values; singularities=[])
    func = interpolate(ωs, values)
    K = length(ωs)
    cumulint = zeros(Float64, K)
    errors = zeros(Float64, K)
    lb = ωs[1] # lower bound of integration
    if isempty(singularities)
        for (i, ω) in enumerate(ωs)
            if i > 1
                int, error = quadgk(x -> func(x), lb, ω)
                cumulint[i] = int
                errors[i] = error
            end
        end
    else
        error("integration with singularities not yet implemented")
    end
end

function interpolate(xs, ys; type="linear")
    if type == "linear"
        return linear_interpolation(xs, ys)
    elseif type == "cubic"
        return cubic_spline_interpolation(xs, ys)
    else
        error("type of interpolation unknown, try linear or cubic")
    end
end