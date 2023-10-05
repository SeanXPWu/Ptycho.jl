using Revise
using Ptycho, KernelAbstractions
using ImageView

function main()
    params = Parameters(
        Voltage = 300,
        Semiangle = 1.03,
        dx = 4.65,
        ScanStep = 31.25,
        ScanAngle = -126,
        Defocus = -13500,
    )

    dps = load_dps("./test_data/", (127, 127))

    println("Diffraction patterns loaded \nInitialsing...")

    recon, trans_exec = initialise_recon(params, dps)

    println("Initialisation finished \nStarting reconstruction...")

    iter_step = 10
    backend = CPU()
    precision = Float32

    rmse_list = Vector{Float64}(undef, iter_step)
    for i = 1:iter_step
        println("Current iteration : ", i)
        @time "Iteration $i finished in" recon, rmse =
            ePIE_iteration(recon, params, dps, trans_exec, backend, precision)
        rmse_list[i] = rmse
        println("Current RMSE : $rmse")
    end

    println("Reconstruction finished")

    obj = recon.Object.ObjectMatrix
    imshow(convert(Array{precision}, abs.(obj)))
    return params, recon, rmse_list, dps, trans_exec
end

params, recon, rmse_list, dps, trans_exec = main()
