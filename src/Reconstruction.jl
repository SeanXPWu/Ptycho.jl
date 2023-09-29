export Probe, Object, UpdateParameters, Reconstruction
export save

struct Probe{T<:Real}
    ProbeMatrix::Array{T,3}
    RecordStep::Integer
end

function init_probe(params::Parameters, dps::DiffractionPatterns)
    x, y = size(dps)[1:2]
    return Probe(Array{ComplexF64,2}(undef, x, y), 1)
end

function init_probe(ds::DataSet)
    return init_probe(ds.Params, ds.DPs)
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

function save_recon(filepath::String, recon::Reconstruction, ext::String)
    if ext == "mat"
        file = matopen(filepath, "w")
    elseif ext == "jld2"
        file = jldopen(filepath, "w")
    end
    write(file, "Probe", recon.Probe)
    write(file, "Object", recon.Object)
    write(file, "Iteration", recon.Iteration)
    write(file, "UpdateParameters", recon.UpdateParameters)
    close(file)
end
