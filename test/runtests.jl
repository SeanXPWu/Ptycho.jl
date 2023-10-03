using Ptycho
using Test

@testset "Ptycho.jl" begin
    # Write your tests here.
    @test Parameters(Voltage=300, Semiangle=1.03, dx=4.65, ScanStep=31.25, ScanAngle=-126, Defocus=-13500) == Parameters(
        300,
        1.03,
        4.65,
        ScanTrajectory(31.25, -126),
        AberrationParameters(-13500),
        0.1,
        0.01,
        50,
    )

end
