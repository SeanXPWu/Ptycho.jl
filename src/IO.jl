export load_dp, load_dps, load_dps_3d, save_dataset, save_recon

function load_dp(filename::String, varname::String = "dp")
    ext = split(filename, ".")[end]
    if ext == "mat"
        file = matopen(filename)
    else
        error("Unsupported file extension: $ext")
    end
    dp = read(file, varname)
    close(file)
    return dp
end

"""
Load diffraction patterns from a collection of data files in a directory.

The scansize is a tuple of the form (x, y) where x and y are the number of scan points in the x and y directions respectively.
"""
function load_dps(
    dirpath::String,
    scansize::T,
    varname::String = "dp",
) where {T<:Tuple{Integer,Integer}}
    files = readdir(dirpath, join = true, sort = false)
    if length(files) != prod(scansize)
        error("Number of files in $dirpath does not match scan dimensions.")
    end
    tmp = load_dp(files[1], varname)
    x, y = Base.size(tmp)
    dps = Array{UInt8,4}(undef, x, y, scansize[1], scansize[2])
    for i =1:length(files)
        filename="$dirpath/dp_$i.mat"
        tmp = load_dp(filename, varname)
        y = (i-1) ÷ scansize[1] + 1
        x = i - (y-1)*scansize[1]
        dps[:, :, x, y] = tmp
    end
    return DiffractionPatterns(dps)
end

function load_dps_3d(
    dirpath::String,
    scansize::T,
    varname::String = "dp",
) where {T<:Tuple{Integer,Integer}}
    files = readdir(dirpath, join = true, sort = false)
    if length(files) != prod(scansize)
        error("Number of files in $dirpath does not match scan dimensions.")
    end
    tmp = load_dp(files[1], varname)
    x, y = Base.size(tmp)
    dps = Array{UInt8,3}(undef, x, y, scansize[1] * scansize[2])
    for i = 1:length(files)
        filename = "$dirpath/dp_$i.mat"
        tmp = load_dp(filename, varname)
        dps[:, :, i] = tmp
    end
    return DiffractionPatterns(dps)
end


function save_dataset(
    filepath::String,
    data,
    ext::String = "jld2",
    mode::String = "Ptycho",
)
    if mode == "Ptycho" && ext != "jld2"
        error("Ptycho struct mode only supports JLD2 format.")
    end
    if ext == "mat"
        file = matopen(filepath, "w")
    elseif ext == "jld2"
        file = jldopen(filepath, "w")
    end
    write(file, "Params", data.Params)
    write(file, "DPs", data.DPs)
    close(file)
end


"""
Save the reconstruction data to a file with the given extension.
"""
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
