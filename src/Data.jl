export Parameters, DiffractionPatterns, DataSet

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
    Aberrations::AberrationParameters{Real}
end

function Parameters(Voltage,Semiangle,dx,ScanStep,Angle,Defocus)
    return Parameters(Voltage,
                      Semiangle,
                      dx,
                      ScanTrajectory(ScanStep,
                                     Angle),
                      AbberationParameters(Defocus)
                      )
end

struct DiffractionPatterns{T<:Real}
    DPs::Array{T,4}
end

struct DataSet
    Params::Parameters
    DPs::DiffractionPatterns
end

function save(filepath::String, data::DataSet)

end
