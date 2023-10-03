export AberrationParameters, ScanTrajectory, Parameters, DiffractionPatterns

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
    ObjUpdate::Float32
    ProbeUpdate::Float32
    ObjPadding::Integer
end

function Parameters(;
    Voltage::Real,
    Semiangle::Real,
    dx::Real,
    ScanStep::Real,
    ScanAngle::Real,
    Defocus::Real,
    ObjUpdate::T=0.1,
    ProbeUpdate::T=0.01,
    ObjPadding::Integer=50,
) where T<:AbstractFloat
    return Parameters(
        Voltage,
        Semiangle,
        dx,
        ScanTrajectory(ScanStep, ScanAngle),
        AberrationParameters(Defocus),
        ObjUpdate,
        ProbeUpdate,
        ObjPadding,
    )
end

struct DiffractionPatterns{T<:Real}
    DPs::AbstractArray{T, 4}
end

function size(dps::DiffractionPatterns)
    return size(dps.DPs)
end

function length(dps::DiffractionPatterns)
    return length(dps.DPs[1,1,:,:])
end

"""function getindex(dps::DiffractionPatterns, range::Vector{UnitRange{Int}})

end"""