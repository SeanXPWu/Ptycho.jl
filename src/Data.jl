export AberrationParameters, ScanTrajectory, Parameters, DiffractionPatterns, DataSet
export load_dp, load_dps, save_dataset, load_dps_3d

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
    Adjustment::Integer
end

function Parameters(Voltage, Semiangle, dx, ScanStep, Angle, Defocus)
    return Parameters(
        Voltage,
        Semiangle,
        dx,
        ScanTrajectory(ScanStep, Angle),
        AberrationParameters(Defocus),
        50
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

"""
Load diffraction patterns from a collection of data files in a directory.

The scansize is a tuple of the form (x, y) where x and y are the number of scan points in the x and y directions respectively.
"""
function load_dps(
    dirpath::String,
    scansize::T,
    varname::String = "dp",
) where {T<:Tuple{Integer,Integer}}
    files = readdir(dirpath, join=true, sort=false)
    if length(files) != prod(scansize)
        error("Number of files in $dirpath does not match scan dimensions.")
    end
    tmp = load_dp(files[1], varname)
    x, y = Base.size(tmp)
    dps = Array{UInt8,4}(undef, x, y, scansize[1], scansize[2])
    for i =1:length(files)
        filename="$dirpath/dp_$i.mat"
        tmp = load_dp(filename, varname)
        y = (i-1) ÷ scansize[1] + 1
        x = i - (y-1)*scansize[1]
        #println(x,y)
        dps[:, :, x, y] = tmp
    end
    return DiffractionPatterns(dps)
end

function load_dps_3d(
    dirpath::String,
    scansize::T,
    varname::String = "dp",
) where {T<:Tuple{Integer,Integer}}
    files = readdir(dirpath, join=true, sort=false)
    if length(files) != prod(scansize)
        error("Number of files in $dirpath does not match scan dimensions.")
    end
    tmp = load_dp(files[1], varname)
    x, y = Base.size(tmp)
    dps = Array{UInt8,3}(undef, x, y, scansize[1]*scansize[2])
    for i =1:length(files)
        filename="$dirpath/dp_$i.mat"
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
