using Revise
using Ptycho
using ImageView

function main()
params = Parameters(300, 1.03, 4.65, 31.25, -126, -13500)
#params = Parameters(300, 4.64, 1.03, 32.55, 0, -2300)
dps =load_dps("./test_data/",(127,127))
#dps =load_dps_3d("../001/",(128,128))

recon, trans_exec, dpList = prestart(params,dps,0.01,0.01)

iter_step = 1
rmse_list = Vector{Float64}(undef, iter_step)
for i = 1:iter_step
    @time recon,rmse = iterate(recon,trans_exec,dps,dpList)
    rmse_list[i] = rmse
end
return params, recon, rmse_list, dps, trans_exec
end

params,recon, rmse_list,dps, trans_exec=main()
obj = recon.Object.ObjectMatrix
imshow(abs.(obj))
