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

function generate_dpList(dps::DiffractionPatterns,trans_related,mode,OtherPara)
    if mode=="all"
    dpList = collect(1:length(dps))
    elseif mode=="frc"
        if mod(length(dps),2) == 0
            dpList = zeros(length(dps)/2,1);
            count = 1;
            for ii = 1:size(dps)[3]
                dpList[count:count-1+size(dps)[4]/2] = ((ii-1)*size(dps)[4].+1:2:ii*size(dps)[4]) .+ mod(ii+OtherPara[1],2)
                count = count + size(dps)[4]/2
            end
        else
            dpList = (1:2:(length(dps))) .+ mod(OtherPara[1]+1,2)
        end
    elseif mode=="rows"
        dpList = zeros(length(OtherPara)*size(dps)[4],1)
        count = 1;
        OtherPara=Int.(OtherPara)
        for ii = OtherPara
            dpList[count:size(dps)[4]+count-1] = ((ii-1)*size(dps)[4]+1:1:ii*size(dps)[4])
            count = count + size(dps)[4]
        end
    elseif mode=="columns"
        dpList = zeros(size(dps)[3]*length(OtherPara),1)
        count = 1
        OtherPara=Int.(OtherPara)
        for ii = OtherPara
            dpList[count:size(dps)[3]+count-1] = (ii:size(dps)[4]:ii+(size(dps)[3]-1)*size(dps)[4])'
            count = count + size(dps)[3]
        end
    elseif mode=="singledp"
        dpList = OtherPara
    elseif mode=="subarray"
        BigMtx = 1:(length(dps))
        BigMtx = reshape(BigMtx,size(dps)[4],size(dps)[3])'
        SmlMtx = BigMtx[OtherPara[1]:OtherPara[2],OtherPara[3]:OtherPara[4]]
        dpList = zeros((OtherPara[4]-OtherPara[3] + 1)*(OtherPara[2]-OtherPara[1]+1),1)
        count = 0;
        for ii = 1:OtherPara[2]-OtherPara[1]+1
            for jj = 1: OtherPara[4]-OtherPara[3] + 1
                count = count + 1;
                dpList[count] = SmlMtx[ii,jj]
            end
        end
    elseif mode=="double-skip"
        NewPxlNum = Int.(ceil.(size(dps)[3:4]./(OtherPara+1)))
        dpList = zeros(NewPxlNum[1]*NewPxlNum[2],1)
        count = 0
        for ii = 1:NewPxlNum[1]
            for jj = 1:NewPxlNum[2]
                count = count + 1
                dpList[count] = (ii-1)*(OtherPara+1)*size(dps)[4] + (jj-1)*(OtherPara+1)+1
            end
        end
    end
    trans_new = Array{Float64, 2}(undef,length(dpList), 2)
    dpList=Int.(dpList)
    if mode=="singledp"
        trans_new[:,1]=trans_related[1,1]
        trans_new[:,2]=trans_related[1,2]
    else
        trans_new[:,1]=trans_related[dpList,1]
        trans_new[:,2]=trans_related[dpList,2]
    end
    return dpList, trans_new
end

function init_obj(recon_size)
    return Object(
        ones(ComplexF32, ceil(Integer, recon_size[1]), ceil(Integer, recon_size[2])),
    )
end

function prestart(params::Parameters, dps::DiffractionPatterns)
    trans_related = init_trans(params, dps)
    dpList, trans_related = generate_dpList(dps,trans_related,"double-skip",1)
    trans_exec = similar(trans_related)
    trans_exec[:, 1] =
        trans_related[:, 1] .- floor(minimum(trans_related[:, 1])) .+ params.ObjPadding
    trans_exec[:, 2] =
        trans_related[:, 2] .- floor(minimum(trans_related[:, 2])) .+ params.ObjPadding
    length_x = maximum(trans_exec[:, 1]) - minimum(trans_exec[:, 1])
    length_y = maximum(trans_exec[:, 2]) - minimum(trans_exec[:, 2])
    recon_size = [ceil(length_x), ceil(length_y)] .+ size(dps)[1:2] .+ params.ObjPadding * 2
    obj = init_obj(recon_size)
    probe = init_probe(params, dps)
    recon = Reconstruction(obj, probe, 0)
    return recon, trans_exec, dpList
end

function rmse(arr1::T, arr2::T) where {T<:AbstractArray}
    return sqrt(sum((arr1 .- arr2) .^ 2))
end

function iterate(recon::Reconstruction,trans_exec,params::Parameters, dps::DiffractionPatterns, dpList)
        current_rmse=0
        obj = recon.Object.ObjectMatrix
        probe = recon.Probe.ProbeMatrix
        for count_dp=1:length(dpList)
                sx = round(Int, trans_exec[count_dp,1]):round(Int,trans_exec[count_dp,1]+size(dps)[1]-1)
                sy = round(Int,trans_exec[count_dp,2]):round(Int,trans_exec[count_dp,2]+size(dps)[2]-1)
                y =(dpList[count_dp]-1) ÷ size(dps)[3]+1
                x = dpList[count_dp]-(y-1)*size(dps)[3]
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
