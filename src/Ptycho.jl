module Ptycho

using FFTW
using LinearAlgebra
using KernelAbstractions

using JLD2, HDF5, MAT

using ImageView

import Base: read, size, iterate

"""
Convert array to the desinated backend and precision.
"""
function to_backend(backend::B, precision::T, arr::A) where {B<:Backend, T<:DataType, A<:Array{<:Any}}
    if backend == CPU()
        return Array{precision}(arr)
    else
        gpuarr = allocate(backend, precision, size(arr))
        copyto!(gpuarr, arr)
        return gpuarr
    end
end

include("Data.jl")
include("Reconstruction.jl")
include("PIE.jl")

include("iteration.jl")

end
