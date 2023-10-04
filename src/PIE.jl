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

"""
Given a 4d data array of the dimension (x, y, z, w), where x and y are the dimensions of the diffraction patterns, z and w are the dimensions of the scan trajectory. Corresponding index notations are (i, j, k, l).
"""
#@kernel function ePIE_iteration_kernel!(dps, probe, object, trans, probeup, objup)
#    k,l = @index(Global, Ntuple)
#    for j in axes(dps, 2)
#        for i in axes(dps, 1)
            # main logic
#        end
#    end
#end

@kernel function get_ew_kernel!(ew, obj, probe, sx, sy)
    i, j = @index(Global, NTuple)
    ew[i,j] = probe[i,j] * obj[sx+i-1,sy+j-1]
end

function get_ew!(ew, obj, probe, sx, sy, backend)
    kernel! = get_ew_kernel!(backend)
    kernel!(ew, obj, probe, sx, sy, ndrange = size(probe))
end

@kernel function probe_update_kernel!(obj, probe, ew, uew, sx, sy, pup, obj_max)
    i, j = @index(Global, NTuple)
    probe[i,j] +=pup * (uew[i,j] - ew[i,j]) * conj(obj[i+sx-1,j+sy-1]) / obj_max
end

function probe_update!(obj, probe, ew, uew, sx, sy, pup, obj_max, backend)
    kernel! = probe_update_kernel!(backend)
    kernel!(obj, probe, ew, uew, sx, sy, pup, obj_max, ndrange=size(probe))
end

@kernel function object_update_kernel!(obj, probe, ew, uew, sx, sy, oup, probe_max)
    i, j = @index(Global, NTuple)
    obj[i+sx-1,j+sy-1] +=oup * (uew[i,j] - ew[i,j]) * conj(probe[i,j]) / probe_max
end

function object_update!(obj, probe, ew, uew, sx, sy, oup, probe_max, backend)
    kernel! = object_update_kernel!(backend)
    kernel!(obj, probe, ew, uew, sx, sy, oup, probe_max, ndrange=size(probe))
end

function ePIE_iteration(recon::Reconstruction, params::Parameters, dps::DiffractionPatterns, trans_exec::AbstractArray, backend::B=CPU(), precision::DataType = Float32) where {B<:Backend}
    complexprecision = Complex{precision}
    obj = to_backend(backend, complexprecision, recon.Object.ObjectMatrix)
    probe = to_backend(backend, complexprecision, recon.Probe.ProbeMatrix)
    dps = dps.DPs
    
    ew = similar(probe)
    fp = plan_fft(ew)
    ifp = plan_ifft(ew)

    current_rmse = 0.0
    dims = size(dps)
    x, y = dims[3:4]
    for j = 1:y
        for i = 1:x
            dp = to_backend(backend, precision, dps[:,:,i,j])
            sx, sy = @view trans_exec[j,i,:]
            get_ew!(ew, obj, probe, sx, sy, backend)

            ewf = ifftshift(fp*(fftshift(ew))) ./ sqrt(length(ew))

            current_rmse+=rmse(abs.(ewf), dp)
            ewfn = dp .* exp.(im.*angle.(ewf))

            uew = ifftshift(ifp*(fftshift(ewfn))) .* sqrt(length(ewfn))

            obj_max = maximum(abs.(obj).^2)
            probe_update!(obj, probe, ew, uew, sx, sy, params.ProbeUpdate, obj_max, backend)
            probe_max = maximum(abs.(probe).^2)    
            object_update!(obj, probe, ew, uew, sx, sy, params.ObjUpdate, probe_max, backend)
        end
    end
    return Reconstruction(Object(obj), Probe(probe), recon.Iteration+1), current_rmse/(x*y)
end
