module Ptycho

using FFTW
using LinearAlgebra
using KernelAbstractions

using JLD2, HDF5, MAT

using ImageView

import Base: size, length, iterate

"""
Convert array to the desinated backend and precision.
"""
function to_backend(
    backend::B,
    precision::T,
    arr::A,
) where {B<:Backend,T<:DataType,A<:Array{<:Any}}
    if backend == KernelAbstractions.get_backend(arr) && precision == eltype(arr)
        return arr
    else
        newarr = allocate(backend, precision, size(arr))
        copyto!(newarr, arr)
        return newarr
    end
end

include("Dataset.jl")
include("Reconstruction.jl")
include("IO.jl")
include("PIE.jl")

include("Iteration.jl")

end
