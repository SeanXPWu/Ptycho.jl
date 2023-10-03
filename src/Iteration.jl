export ePIE
export run_iteration

abstract type Engine end

abstract type PIE <: Engine end

struct ePIE <: PIE end

function run_iteration(epie::ePIE, recon::Reconstruction, trans_exec::AbstractArray, params::Parameters, dps::DiffractionPatterns)
    rmse = 0.0
    obj = recon.Object.ObjectMatrix
    probe = recon.Probe.ProbeMatrix
    dps = dps.DPs
    ePIE_iteration!()
end
