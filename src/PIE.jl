function trans(params::Parameters,dps::DiffractionPatterns)
    scanstep = params.Scan.ScanStep
    angle = params.Scan.Angle
    dx = params.dx
    x, y = size(dps)[3,4]
    trans = Array{Float64, 3}(undef,x, y, 2)
    for j = 1:y
        for i = 1:x
            trans[i, j, 1] = (i * scanstep / dx * cos(angle * pi / 180) - j * scanstep / dx * sin(angle * pi / 180))
            trans[i, j, 2] = (i * scanstep / dx * sin(angle * pi / 180) + j * scanstep / dx * cos(angle * pi / 180))
        end
    end
    return trans
end

@kernel function ePIE_iteration_kernel!()

end

function ePIE_iteration!()

end
