export DM_iterations

function fourier_projection!(Ψ_fp, Ψ, dps, fplan, ifplan, backend)
    fft_mul! = fft_mul_kernel!(backend)
    fftshift!(Ψ_fp, Ψ)
    fft_mul!(Ψ, fplan, Ψ_fp, ndrange=size(Ψ))
    ifftshift!(Ψ_fp, Ψ)
    dot!(Ψ, dps, exp.(im.*angle.(Ψ_fp)))
    fftshift!(Ψ_fp, Ψ)
    fft_mul!(Ψ, ifplan, Ψ_fp, ndrange=size(Ψ))
    ifftshift!(Ψ_fp, Ψ)
end

@kernel function fft_mul_kernel!(result, plan, arr)
    i,j,k,l = @index(Global, NTuple)
    tmp = zero(eltype(arr))
    for m in 1:size(arr, 1)
        tmp += plan[i,m] * arr[m,j,k,l]
    end
    result[i,j,k,l] = tmp
end

@kernel function overlap_projection_kernel!(Ψ, probe, object, trans)
    i,j,k,l = @index(Global, NTuple)
    sx, sy = trans[l, k, :]
    Ψ[i,j,k,l] = probe[i,j] * object[i+sx-1,j+sy-1]
end

@kernel function update_kernel!(Ψ, probe, object, object_tmp, trans, x, y)
    i,j = @index(Global, NTuple)
    probe_tmp = 0.0
    for l in 1:y
        for k in 1:x
            sx, sy = trans[l, k, :]
            object_tmp[i+sx-1,j+sy-1] += conj(probe[i,j]) * Ψ[i,j,k,l] / abs(probe[i,j])^2
            probe_tmp += conj(object[i+sx-1,j+sy-1]) * Ψ[i,j,k,l] / abs(object[i+sx-1,j+sy-1])^2
        end
    end
    object[i+sx-1,j+sy-1] = object_tmp[i+sx-1,j+sy-1]
    probe[i,j] = probe_tmp
end

function DM_iterations(params::Parameters,
    dps::DiffractionPatterns,
    iterations::Integer,
    backend::B = CPU(),
    precision::DataType = Float32,
) where {B<:Backend}
    println("Initialsing...")
    @time "Initialisation finished in" begin
    recon, trans = initialise_recon(params, dps)
    a, b, x, y = size(dps)

    complexprecision = Complex{precision}
    object = to_backend(backend, complexprecision, recon.Object.ObjectMatrix)
    probe = to_backend(backend, complexprecision, recon.Probe.ProbeMatrix)
    dps = to_backend(backend, eltype(dps.DPs), dps.DPs)
    obj_tmp = similar(object)

    overlap_projection! = overlap_projection_kernel!(backend)
    update! = update_kernel!(backend)

    Ψ = allocate(backend, complexprecision, (a, b, x, y))
    Ψ_diff = similar(Ψ)
    Ψ_op = similar(Ψ)
    Ψ_fp = similar(Ψ)
    overlap_projection!(Ψ, probe, object, trans, ndrange = (a,b,x,y))

    fplan = plan_fft(Ψ[:,:,1,1])
    ifplan = plan_ifft(Ψ[:,:,1,1])

    dm_error = Vector{Float64}(undef, iterations)
    end

    println("Start reconstruction...")
    for iter = 1:iterations
        println("Current iteration : ", iter)
        @time "Iteration finished in" begin
            overlap_projection!(Ψ_op, probe, object, trans, ndrange = (a,b,x,y))
            Ψ_diff .= fourier_projection!(Ψ_fp, (2 .* Ψ_op .- Ψ), dps, fplan, ifplan, backend) .- Ψ_op
            dm_error[iter] = sum(abs.(Ψ_diff))
            Ψ .+= Ψ_diff
            object_tmp .= 0.0
            update!(Ψ, probe, object, obj_tmp, trans, x, y, ndrange = (a,b))
            copyto!(object, obj_tmp)
        end
    end
    println("Reconstruction finished")
    return Reconstruction(Object(object), Probe(probe), iterations), dm_error
end
