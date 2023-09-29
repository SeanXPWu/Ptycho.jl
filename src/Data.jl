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

struct DataSet
    Params::Parameters
    DPs::DiffractionPatterns
end

function load_dp(filename::String, varname::String)
    ext = split(filename,".")[end]
    if ext == "mat"
        file = matopen(filename)
    end
    dp = read(file, varname)
    close(file)
    return dp
end

function save(filepath::String, data::DataSet)

end
