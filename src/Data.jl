struct AberrationParameters{T} where T<:Real
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
    Scan::ScanTrajectory{Real}
    Aberrations::AberrationParameters{Real}
end

struct DiffractionPatterns{T} where T<:Real
    DPs::Array{T,4}
end

struct DataSet
    Params::Parameters
    DPs::DiffractionPatterns
end

function read(filename::String)

    return DiffractionPatterns{}()
end

function save(filepath::String, data::DataSet)

end
