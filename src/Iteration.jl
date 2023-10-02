abstract type Engine end

abstract type PIE <: Engine end

struct ePIE <: PIE end

function run_iteration(epie::ePIE)
    ePIE_iteration!()
end
