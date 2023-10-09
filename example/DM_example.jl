using Revise
using Ptycho
using ImageView

function main()
    params = Parameters(
        Voltage = 300,
        Semiangle = 1.03,
        dx = 4.65,
        ScanStep = 31.25,
        ScanAngle = -126,
        Defocus = -13500,
        ObjUpdate = 1e-4
    )

    dps = load_dps("./test_data/", (127, 127))

    println("Diffraction patterns loaded")

    iter_step = 2
    #backend = CUDABackend()
    #precision = Float32

    recon, _ = DM_iterations(params, dps, iter_step)

    obj = recon.Object.ObjectMatrix

    obj = convert(Array{Float32}, obj)
    imshow(abs.(obj))
    imshow(angle.(obj))

    #return params, recon, dps, rmse
end

main()
