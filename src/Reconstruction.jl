struct Probe{T<:Real}
    ProbeMatrix::Array{T,3}
    RecordStep::Integer
end

function init_probe()

end

struct Object{T<:Real}
    ObjectMatrix::Array{T,3}
    RecordStep::Integer
end

struct UpdateParameters

end

struct Reconstruction
    Probe::Probe
    Object::Object
    Iteration::Integer
    UpdateParameters::UpdateParameters
end

function read(filename::String)

    return Reconstruction{}()
end

function save(filepath::String, recon::Reconstruction)

end
