export initialise_recon, ePIE_iterations, get_rotational_offset

function get_rotational_offset(params::Parameters, dps::DiffractionPatterns)
    scanstep = params.Scan.ScanStep
    angle = params.Scan.Angle
    dx = params.dx
    x, y = size(dps)[3:4]
    trans = Array{Float64,3}(undef, x, y, 2)
    for j = 1:y
        for i = 1:x
            trans[i, j, 1] = (
                i * scanstep / dx * cos(angle * pi / 180) -
                j * scanstep / dx * sin(angle * pi / 180)
            )
            trans[i, j, 2] = (
                i * scanstep / dx * sin(angle * pi / 180) +
                j * scanstep / dx * cos(angle * pi / 180)
            )
        end
    end
    return trans
end

function initialise_recon(params::Parameters, dps::DiffractionPatterns)
    trans = get_rotational_offset(params, dps)
    trans_exec = similar(trans)
    trans_exec[:, :, 1] =
        trans[:, :, 1] .- floor(minimum(trans[:, :, 1])) .+ params.ObjPadding
    trans_exec[:, :, 2] =
        trans[:, :, 2] .- floor(minimum(trans[:, :, 2])) .+ params.ObjPadding
    length_x = maximum(trans_exec[:, :, 1]) - minimum(trans_exec[:, :, 1])
    length_y = maximum(trans_exec[:, :, 2]) - minimum(trans_exec[:, :, 2])
    recon_size = [ceil(length_x), ceil(length_y)] .+ size(dps)[1:2] .+ params.ObjPadding * 2
    obj = init_obj(recon_size)
    probe = init_probe(params, dps)
    return Reconstruction(obj, probe, 0), round.(Int, trans_exec)
end

@kernel function abs_max_kernel!(abs_max, obj, probe)
    i, j = @index(Global, NTuple)
    tmp = abs(obj[i, j])^2
    if abs_max[1] < tmp
        abs_max[1] = tmp
    end
    tmp = abs(probe[i, j])^2
    if abs_max[2] < tmp
        abs_max[2] = tmp
    end
end

@kernel function update_kernel!(obj, probe, ew, uew, oup, pup, abs_max)
    i, j = @index(Global, NTuple)
    obj[i, j] += oup * (uew[i, j] - ew[i, j]) * conj(probe[i, j]) / abs_max[2]
    probe[i, j] += pup * (uew[i, j] - ew[i, j]) * conj(obj[i, j]) / abs_max[1]
end

function ePIE_iterations(
    params::Parameters,
    dps::DiffractionPatterns,
    iterations::Integer,
    backend::B = CPU(),
    precision::DataType = Float32,
) where {B<:Backend}
    println("Initialsing...")
    @time "Initialisation finished in" begin
    recon, trans_exec = initialise_recon(params, dps)
    a, b, x, y = size(dps)
    sqlen = sqrt(a * b)
    complexprecision = Complex{precision}
    object = to_backend(backend, complexprecision, recon.Object.ObjectMatrix)
    probe = to_backend(backend, complexprecision, recon.Probe.ProbeMatrix)
    dps = to_backend(backend, eltype(dps.DPs), dps.DPs)

    ew = similar(probe)
    ewf = similar(probe)
    uewf = similar(probe)
    uew = similar(probe)
    tmp = similar(probe)

    #uew = similar(probe)
    fp = plan_fft(ew)
    ifp = plan_ifft(ew)

    abs_max = allocate(backend, precision, 2)

    objup = params.ObjUpdate
    probeup = params.ProbeUpdate

    abs_max! = abs_max_kernel!(backend)
    update! = update_kernel!(backend)

    rmse_ls = Vector{Float64}(undef, iterations)
    end

    println("Start reconstruction...")
    for iter = 1:iterations
        println("Current iteration : ", iter)
        @time "Iteration finished in" begin
            current_rmse = 0.0
            abs_max .= 0.0
            for j = 1:y
                for i = 1:x
                    sx, sy = @view trans_exec[j, i, :]
                    dp = view(dps, :, :, i, j)
                    obj = view(object, sx:sx+a-1, sy:sy+b-1)
                    
                    dot!(ew, obj, probe)

                    fftshift!(ewf, ew)
                    mul!(tmp, fp, ewf)
                    ifftshift!(ewf, tmp)
                    ewf ./= sqlen

                    current_rmse += rmse(abs.(ewf), dp)
                    dot!(uewf, dp, exp.(im .* angle.(ewf)))

                    fftshift!(uew, uewf)
                    mul!(tmp, ifp, uew)
                    ifftshift!(uew, tmp)
                    uew .*= sqlen

                    abs_max!(abs_max, obj, probe, ndrange = size(probe))
                    update!(obj, probe, ew, uew, objup, probeup, abs_max, ndrange = size(probe))
                end
            end
            current_rmse /= x * y
            rmse_ls[iter] = current_rmse
            println("Current RMSE: $current_rmse")
        end
    end
    println("Reconstruction finished")
    return Reconstruction(Object(object), Probe(probe), iterations), rmse_ls
end
