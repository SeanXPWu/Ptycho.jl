"""
Convert array to the desinated backend and precision.
"""
function to_backend(
    backend::B,
    precision::T,
    arr::A,
) where {B<:Backend,T<:DataType,A<:AbstractArray}
    if backend == KernelAbstractions.get_backend(arr) && precision == eltype(arr)
        return arr
    else
        newarr = allocate(backend, precision, size(arr))
        copyto!(newarr, arr)
        return newarr
    end
end

@kernel function dot_kernel!(result, a, b)
    i = @index(Global, Linear)
    result[i] = a[i] * b[i]
end

function dot!(result, a, b)
    kernel! = dot_kernel!(KernelAbstractions.get_backend(a))
    kernel!(result, a, b, ndrange = size(result))
end
