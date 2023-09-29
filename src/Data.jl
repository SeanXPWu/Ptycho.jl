export AberrationParameters, ScanTrajectory, Parameters, DiffractionPatterns, DataSet
export load_dp

struct AberrationParameters{T<:Real}
    Defocus::T
end

struct ScanTrajectory
    ScanStep::Real
    Angle::Real
end

struct Parameters
    Voltage::Real
    Semiangle::Real
    dx::Real
    Scan::ScanTrajectory
    Aberrations::AberrationParameters
end

function Parameters(Voltage, Semiangle, dx, ScanStep, Angle, Defocus)
    return Parameters(
        Voltage,
        Semiangle,
        dx,
        ScanTrajectory(ScanStep, Angle),
        AberrationParameters(Defocus),
    )
end

struct DiffractionPatterns{T<:Real}
    DPs::Array{T,4}
end

function size(dps::DiffractionPatterns)
    return size(dps.DPs)
end

struct DataSet
    Params::Parameters
    DPs::DiffractionPatterns
end

function load_dp(filename::String, varname::String)
    ext = split(filename, ".")[end]
    if ext == "mat"
        file = matopen(filename)
    else
        error("Unsupported file extension: $ext")
    end
    dp = read(file, varname)
    close(file)
    return dp
end

function load_dps(dirpath::String, varname::String = "dp")
    files = readdir(dirpath)
    tmp = load_dp(files[1], varname)
    x, y = Base.size(tmp)
    dps = Array{UInt8,3}(undef, x, y, length(files))
    for (i, filename) in enumerate(files)
        tmp = load_dp(filename, varname)
        dps[:, :, i] = tmp
    end
    return DiffractionPatterns(dps)
end

function save_dataset(filepath::String, data::DataSet, ext::String = "mat")
    if ext == "mat"
        file = matopen(filepath, "w")
    elseif ext == "jld2"
        file = jldopen(filepath, "w")
    end
    write(file, "Params", data.Params)
    write(file, "DPs", data.DPs)
    close(file)
end
