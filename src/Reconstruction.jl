export Probe, Object, UpdateParameters, Reconstruction
export init_probe, save_recon

struct Probe{T<:Real}
    ProbeMatrix::Array{T,3}
    RecordStep::Integer
end

function Ptycho_fft2(x)
    xif = ifftshift(ifft(fftshift(x)))*sqrt.(length(x));
    return(xif)
end

function init_probe(params::Parameters, dps::DiffractionPatterns)
    x, y = size(dps)[1:2]
    lamda= 1.23984244/sqrt(params.Voltage*(2*510.99906+params.Voltage))*1e-9
    temp1= (-1/2*x):(1/2*x-1)
    temp1=temp1*1e10*lamda/p.dp.dxy/p.dp.size(1)
    temp2= (-1/2*y):(1/2*y-1)
    temp2=temp2*1e10*lamda/params.dx/y
    k1=temp1' .* ones(length(temp2))
    k2 = ones(length(temp1))' .* temp2
    w=k1+im*k2
    wi=k1-im*k2
    C1   = params.Aberrations.Defocus*exp(0)*1e-9
    kai=real(1/2*w.*wi*C1)
    aberr = -2*pi/lamda*kai
    apeture=(k1.^2+k2.^2)<=(params.Semiangle*1e-3).^2
    temp=sqrt((k1.^2+k2.^2))
    Nedge=2
    dEdge = Nedge*lamda*1e10/(params.Semiangle*1e-3)/x/params.dx;
    ind = findall((i/(params.Semiangle*1e-3) > 1-dEdge) && (i/(params.Semiangle*1e-3) < 1+dEdge) for i in temp)
    apeture[ind] = 0.5*(1-sin(pi/(2*dEdge)*(temp(ind)/(params.Semiangle*1e-3)-1)))
    probe = exp(im.*aberr).*apeture
    probe = Ptycho_fft2(probe);
    probe=probe/sqrt(sum(sum(abs(probe.*conj(probe)))))*sqrt(sum(sum(abs(dp_ref.*conj(dp_ref)))))
    return Probe(Array{ComplexF64,2}(undef, x, y), 1)
end

function init_probe(ds::DataSet)
    return init_probe(ds.Params, ds.DPs)
end

struct Object{T<:Real}
    ObjectMatrix::Array{T,3}
    RecordStep::Integer
end

struct UpdateParameters

end

struct Reconstruction
    Probe::Probe
    Object::Object
    Iteration::Integer
    UpdateParameters::UpdateParameters
end

function save_recon(filepath::String, recon::Reconstruction, ext::String)
    if ext == "mat"
        file = matopen(filepath, "w")
    elseif ext == "jld2"
        file = jldopen(filepath, "w")
    end
    write(file, "Probe", recon.Probe)
    write(file, "Object", recon.Object)
    write(file, "Iteration", recon.Iteration)
    write(file, "UpdateParameters", recon.UpdateParameters)
    close(file)
end
