using Revise
using Ptycho
using ImageView

function main()
    params = Parameters(Voltage=300, Semiangle=1.03, dx=4.65, ScanStep=31.25, ScanAngle=-126, Defocus=-13500)

    dps =load_dps("./test_data/",(127,127))

    println("Diffraction patterns loaded \nInitialsing...")

    recon, trans_exec, _ = prestart(params, dps)    

    println("Initialisation finished \nStart iteration...")

    iter_step = 1
    rmse_list = Vector{Float64}(undef, iter_step)
    for i = 1:iter_step
        println("Current iteration : ", i)
        @time recon,rmse = iterate(recon, trans_exec, params, dps)
        rmse_list[i] = rmse
        println("Iteration $i finished, current RMSE : ", rmse)
    end

    println("Reconstruction finished")

    return params, recon, rmse_list, dps, trans_exec
end

params,recon, rmse_list,dps, trans_exec=main()
obj = recon.Object.ObjectMatrix
imshow(abs.(obj))
