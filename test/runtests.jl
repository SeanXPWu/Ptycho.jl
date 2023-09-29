using Ptycho
using Test

@testset "Ptycho.jl" begin
    # Write your tests here.
    @test Parameters(300, 1.03, 4.65, 31.25, -126, -13500) ==
          Parameters(300, 1.03, 4.65, ScanTrajectory(31.25, -126), AberrationParameters(-13500))
    @test init_probe(Parameters(300, 1.03, 4.65, 31.25, -126, -13500),DiffractionPatterns(load_dp("../test_data/dp_1.mat","dp"))).ProbeMatrix == load_dp("../probe.mat","probe")
end
