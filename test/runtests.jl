using Ptycho
using Test

@testset "Ptycho.jl" begin
    # Write your tests here.
    @test Parameters(300, 4.64, 1.03, 32.55, 0, -2300) ==
          Parameters(300, 4.64, 1.03, ScanTrajectory(32.55, 0), AberrationParameters(-2300))
    @test init_probe(Parameters(300, 4.64, 1.03, 32.55, 0, -2300),DiffractionPatterns(load_dp("../test_data/dp_1.mat","dp"))).ProbeMatrix == load_dp("../probe.mat","probe")
end
