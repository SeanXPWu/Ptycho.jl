export initialise_recon, ePIE_iteration, get_rotational_offset

function get_rotational_offset(params::Parameters,dps::DiffractionPatterns)
    scanstep = params.Scan.ScanStep
    angle = params.Scan.Angle
    dx = params.dx
    x, y = size(dps)[3:4]
    trans = Array{Float64, 3}(undef,x, y, 2)
    for j = 1:y
        for i = 1:x
            trans[i, j, 1] = (i * scanstep / dx * cos(angle * pi / 180) - j * scanstep / dx * sin(angle * pi / 180))
            trans[i, j, 2] = (i * scanstep / dx * sin(angle * pi / 180) + j * scanstep / dx * cos(angle * pi / 180))
        end
    end
    return trans
end

function initialise_recon(params::Parameters,dps::DiffractionPatterns)
    trans = get_rotational_offset(params,dps)
    trans_exec = similar(trans)
    trans_exec[:, :, 1] = trans[:, :, 1] .- floor(minimum(trans[:, :, 1])) .+ params.ObjPadding
    trans_exec[:, :, 2] = trans[:, :, 2] .- floor(minimum(trans[:, :, 2])) .+ params.ObjPadding
    length_x = maximum(trans_exec[:, :, 1]) - minimum(trans_exec[:, :, 1])
    length_y = maximum(trans_exec[:, :, 2]) - minimum(trans_exec[:, :, 2])
    recon_size = [ceil(length_x), ceil(length_y)] .+ size(dps)[1:2] .+ params.ObjPadding * 2
    obj = init_obj(recon_size)
    probe = init_probe(params, dps)
    return Reconstruction(obj, probe, 0), round.(Int, trans_exec)
end

@kernel function ePIE_iteration_kernel!()

end

function ePIE_iteration(recon::Reconstruction, params::Parameters, dps::DiffractionPatterns, trans_exec::AbstractArray, backend::B=CPU(), precision::DataType = Float32) where {B<:Backend}
    obj = recon.Object.ObjectMatrix
    probe = to_backend(backend, Complex{precision}, recon.Probe.ProbeMatrix)
    dps = dps.DPs
    current_rmse = 0.0
    dims = size(dps)
    x, y = dims[3:4]
    for j = 1:y
        for i = 1:x
            dp = to_backend(backend, precision, dps[:,:,i,j])
            trans = @view trans_exec[j,i,:]
            sx = trans[1]:(trans[1]+dims[1]-1)
            sy = trans[2]:(trans[2]+dims[2]-1)
            uobj = to_backend(backend, Complex{precision}, obj[sx,sy])
            ew = uobj .* probe
            ewf = Ptycho_fft2(ew)
            ewfn = dp .* exp.(im.*angle.(ewf))
            uew = Ptycho_ifft2(ewfn)
            probe .+= params.ProbeUpdate.*(uew.-ew) .* conj.(uobj) ./maximum(abs.(uobj).^2)
            uobj .+= params.ObjUpdate .*(uew.-ew) .* conj.(probe) ./maximum(abs.(probe).^2)
            obj[sx,sy] .= to_backend(CPU(), Complex{precision}, uobj)
            current_rmse+=rmse(abs.(ewf),dp)
        end
    end
    current_rmse /= x*y
    return Reconstruction(Object(obj), Probe(probe), recon.Iteration+1), current_rmse
end
