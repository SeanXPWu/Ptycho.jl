export Probe, Object, Reconstruction
export init_probe, init_trans, generate_dpList, prestart, iterate

struct Probe{T<:Complex}
    ProbeMatrix::AbstractArray{T,2}
end

function Ptycho_fft2(x)
    xif = ifftshift(fft(fftshift(x))) ./ sqrt(length(x))
    return (xif)
end

function Ptycho_ifft2(x)
    xif = ifftshift(ifft(fftshift(x))) .* sqrt(length(x))
    return (xif)
end

function init_probe(params::Parameters, dps::DiffractionPatterns)
    x, y = size(dps)[1:2]
    lamda = 1.23984244 / sqrt(params.Voltage * (2 * 510.99906 + params.Voltage)) * 1e-9
    temp1 = (-1/2*x):(1/2*x-1)
    temp1 = temp1 * 1e10 * lamda / params.dx / x
    temp2 = (-1/2*y):(1/2*y-1)
    temp2 = temp2 * 1e10 * lamda / params.dx / y
    k1 = temp1' .* ones(length(temp2))
    k2 = ones(length(temp1))' .* temp2
    w = k1 + im * k2
    wi = k1 - im * k2
    C1 = params.Aberrations.Defocus * exp(0) * 1e-9
    kai = real(1 / 2 * w .* wi * C1)
    aberr = -2 * pi / lamda * kai
    aperture = zeros(x, y)
    aperture[(k1 .^ 2+k2 .^ 2).<=(params.Semiangle*1e-3)^2] .= 1
    temp = sqrt.((k1 .^ 2 + k2 .^ 2))
    Nedge = 2
    dEdge = Nedge * lamda * 1e10 / (params.Semiangle * 1e-3) / x / params.dx
    ind = findall(
        (i / (params.Semiangle * 1e-3) > 1 - dEdge) &&
        (i / (params.Semiangle * 1e-3) < 1 + dEdge) for i in temp
    )
    aperture[ind] .=
        0.5 * (1 .- sin.(pi / (2 * dEdge) * (temp[ind] / (params.Semiangle * 1e-3) .- 1)))
    probe = exp.(im .* aberr) .* aperture
    probe = Ptycho_fft2(probe)
    probe =
        probe ./ sqrt.(sum(sum(abs.(probe .* conj.(probe))))) *
        sqrt.(sum(sum(abs.(dps.DPs[:, :, 1, 1] .* conj.(dps.DPs[:, :, 1, 1])))))
    return Probe(probe)
end

function init_probe(ds::DataSet)
    return init_probe(ds.Params, ds.DPs)
end

struct Object{T<:Complex}
    ObjectMatrix::AbstractArray{T,2}
end

struct Reconstruction
    Object::Object
    Probe::Probe
    Iteration::Integer
end

function init_trans(params::Parameters, dps::DiffractionPatterns)
    count = 0
    x, y = size(dps)[3:4]
    trans_related = Array{Float64, 2}(undef,x * y, 2)
    for i = 1:x
        for j = 1:y
            count = count + 1
            trans_related[count, 1] =
                i * params.Scan.ScanStep / params.dx * cos(params.Scan.Angle * pi / 180) -
                j * params.Scan.ScanStep / params.dx * sin(params.Scan.Angle * pi / 180)
            trans_related[count, 2] =
                i * params.Scan.ScanStep / params.dx * sin(params.Scan.Angle * pi / 180) +
                j * params.Scan.ScanStep / params.dx * cos(params.Scan.Angle * pi / 180)
        end
    end
    return trans_related
end

function generate_dpList(dps::DiffractionPatterns)
    dpList = collect(1:length(dps))
    return dpList
end

function init_obj(recon_size)
    return Object(
        ones(ComplexF32, ceil(Integer, recon_size[1]), ceil(Integer, recon_size[2])),
    )
end

function prestart(params::Parameters, dps::DiffractionPatterns)
    trans_related = init_trans(params, dps)
    dpList = generate_dpList(dps)
    trans_exec = similar(trans_related)
    trans_exec[:, 1] =
        trans_related[:, 1] .- floor(minimum(trans_related[:, 1])) .+ params.Adjustment
    trans_exec[:, 2] =
        trans_related[:, 2] .- floor(minimum(trans_related[:, 2])) .+ params.Adjustment
    length_x = maximum(trans_exec[:, 1]) - minimum(trans_exec[:, 1])
    length_y = maximum(trans_exec[:, 2]) - minimum(trans_exec[:, 2])
    recon_size = [ceil(length_x), ceil(length_y)] .+ size(dps)[1:2] .+ params.Adjustment * 2
    obj = init_obj(recon_size)
    probe = init_probe(params, dps)
    recon = Reconstruction(obj, probe, 0)
    return recon, trans_exec, dpList
end

function rmse(arr1::T, arr2::T) where {T<:AbstractArray}
    return sqrt(sum((arr1 .- arr2) .^ 2))
end

function iterate(recon::Reconstruction,trans_exec,dps::DiffractionPatterns)
        current_rmse=0
        obj = recon.Object.ObjectMatrix
        probe = recon.Probe.ProbeMatrix
        for count_dp=1:length(dps)
                sx = round(Int, trans_exec[count_dp,1]):round(Int,trans_exec[count_dp,1]+size(dps)[1]-1)
                sy = round(Int,trans_exec[count_dp,2]):round(Int,trans_exec[count_dp,2]+size(dps)[2]-1)
                y =(count_dp-1) ÷ size(dps)[3]+1
                x = count_dp-(y-1)*size(dps)[3]
                #dp_current=Float32.(dps.DPs[:,:,count_dp])
                dp_current=Float32.(dps.DPs[:,:,x,y])
                ew = obj[sx,sy] .* probe
                ewf = Ptycho_fft2(ew)
                ewfn = dp_current .* exp.(im.*angle.(ewf))
                ew1 = Ptycho_ifft2(ewfn)
                probe0 = probe
                probe = probe .+ params.ProbeUpdate.*(ew1.-ew) .* conj.(obj[sx,sy]) ./maximum(abs.(obj[sx,sy]).^2)
                obj[sx,sy] = obj[sx,sy] .+ params.ObjUpdate .*(ew1.-ew) .* conj.(probe0) ./maximum(abs.(probe0).^2)
                current_rmse=current_rmse+rmse(abs.(ewf),Float64.(dp_current))
        end
        current_rmse=current_rmse/length(dps.DPs[1,1,:,:])
        return Reconstruction(Object(obj), Probe(probe), recon.Iteration+1), current_rmse
end
