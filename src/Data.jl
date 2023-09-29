export AberrationParameters, ScanTrajectory, Parameters, DiffractionPatterns, DataSet
export load_dp, save_dataset

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
    DPs::Array{T}
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

function load_dps(
    dirpath::String,
    scansize::T,
    varname::String = "dp",
) where {T<:Tuple{Integer,Integer}}
    files = readdir(dirpath)
    if length(files) != prod(scansize)
        error("Number of files in $dirpath does not match scan dimensions.")
    end
    tmp = load_dp(files[1], varname)
    x, y = Base.size(tmp)
    dps = Array{UInt8,3}(undef, x, y, scansize[1], scansize[2])
    for (i, filename) in enumerate(files)
        tmp = load_dp(filename, varname)
        dps[:, :, i%scansize(1), i÷scansize(1)+1] = tmp
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
