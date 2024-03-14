using DelimitedFiles
using ChannelMixingDiscretization
includet("discretization.jl")

function broadenkk(specin, specoutre, specoutim, paramdict; dir="spectral", prefix="spec_FDM_dens_")
    Nz = paramdict["Nz"]
    min = paramdict["broaden_min"]
    max = paramdict["broaden_max"]
    ratio = paramdict["broaden_ratio"]
    T = paramdict["T"]
    alpha = paramdict["broaden_alpha"]
    gamma = paramdict["broaden_gamma"]

    name_bin = prefix * specin * ".bin"

    specoutre_dat = specoutre * ".dat" 
    specoutim_dat = specoutim * ".dat" 
    mkpath(dir)
    if !isfile("$dir/$specoutim_dat")
        run(`broaden -x $gamma -m $min -M $max -r $ratio $name_bin $Nz $alpha $T 1e-100`)
        # real part
        run(`kk spec.dat spec-re.dat`)
        run(`mv spec.dat $dir/$specoutim_dat`)
        run(`mv spec-re.dat $dir/$specoutre_dat`)
    else
        println(specoutim_dat * " already exists!")
    end
end

function readspectral(realfile, imagfile; dir="spectral")
    real = readdlm(dir * "/$realfile")
    imag = readdlm(dir * "/$imagfile")

    return -pi*(real[:, 2] .+ im*imag[:, 2])
end

function discretize(Δint, zs, J, mesh_min, mesh_ratio)
    ω0 = big"1e-100"
    D0 = big"1."

    posfreqs = BigFloat[]
    let ω = D0
        while ω > mesh_min
            pushfirst!(posfreqs, ω)
            ω /= mesh_ratio
        end
    end

    # println(posfreqs)
    posfreqs = vcat([ω0], posfreqs)
    ωs = vcat(-reverse(posfreqs), posfreqs)

    values = Δint(ωs)
    # this works for one channel!!!
    weights = values
    @show length(ωs)
    open("omegas.dat", "w") do io
        writedlm(io, ωs)
    end
    # @show weights

    posdiscs = calculatediscretizers(ωs, values, weights, 1, :fixed)
    negdiscs = calculatediscretizers(ωs, values, weights, -1, :fixed)
    @info "Evaluating coefficients!"
    Ts, Es = evaluatecoefficients(posdiscs, negdiscs, Δint, zs, J)
    starH = StarHamiltonian(Ts, Es, zs)

    chains = maptochainsONEBAND(starH; m=2J)
    savechains(chains)
    ChannelMixingDiscretization.nrgfilesONEBAND(chains)
end

function nrg!(solver, path)
    # 1. Perform the discretization of the current hybridization Δ using ChannelMixingDiscretization: output are the Wilson chain coeffients for each z-number.
    # 2. Use NRG Ljubljana to solve the impurity problem.
    # 3. Post-process the results of the impurity solver.
    paramdict = solver.paramdict
    
    Nz = parse(Int, paramdict["Nz"])
    J = parse(Int, paramdict["J"])
    mesh_max = solver.mesh_max
    mesh_min = solver.mesh_min
    mesh_ratio = solver.mesh_ratio

    paramfile = solver.paramfile

    mkpath(path)
    cp(paramfile, joinpath(path, paramfile), force=true)
    cp("model.m", joinpath(path, "model.m"), force=true)
    cp("modeloperators.m", joinpath(path, "modeloperators.m"), force=true)
    cd(path)
    
    # 1. Discretize the hybridization function
    # save input Delta.dat
    open("Delta.dat", "w") do io
        writedlm(io, [solver.ωs real.(solver.Δ) imag.(solver.Δ)])
    end
    
    # scalar case!
    Δ = -imag.(solver.Δ)
    # Scale and fix Delta min!
    scale = convert(Float64, mesh_max)
    xs = solver.ωs ./ scale
    ys = Δ .* scale
    ys[ys .< solver.Delta_min] .= solver.Delta_min
    # interpolate hybridization function rescaled to [-1., 1.]
    Δint = LinearInterpolation(ys, xs)
    zs = range(1/Nz, 1, Nz)
    discretize(Δint, zs, J, mesh_min, mesh_ratio)

    # 2. Solve the impurity problem
    # Write the current impurity energy to the param file (eps = -μ)
    ε = -solver.μ
    for (i, z) in enumerate(zs)
        @info "Performing NRG run for i = $(i), z = $(z)"
        cmd_m4 = `m4 -DEPS=$ε -DZZZ=$z param.m4`
        write("param", read(cmd_m4))
        mv("param", "$i/param")
        cd("$i")
        run(`nrginit`)
        run(`nrgrun`)
        cd("..")
    end

    # 3. Process the results (broadening of spectral functions and Σ-trick)
    # Broadening
    @info "Broadening spectral functions!"
    broadenkk("A_d-A_d", "c-reG", "c-imG", paramdict)
    broadenkk("self_d-A_d", "c-reFl", "c-imFl", paramdict)
    broadenkk("A_d-self_d", "c-reFr", "c-imFr", paramdict)
    broadenkk("self_d-self_d", "c-reI", "c-imI", paramdict)
    
    # Gather expectation values
    @info "Gathering expectation values!"
    run(`gatherlastlines custom`)
    write("custom.avg", read(`columnavg_comment custom`))
    run(`gatherlastlines customfdm`)
    write("customfdm.avg", read(`columnavg_comment customfdm`))

    expectdir = "expectation"
    mkpath(expectdir)

    # move them to a separate folder

    mv("custom", expectdir * "/custom", force=true)
    mv("custom.avg", expectdir * "/custom.avg", force=true)
    mv("customfdm", expectdir * "/customfdm", force=true)
    mv("customfdm.avg", expectdir * "/customfdm.avg", force=true)
    
    @info "Performing the Σ-trick!"

    G = readspectral("c-reG.dat", "c-imG.dat")
    Fl = readspectral("c-reFl.dat", "c-imFl.dat")
    Fr = readspectral("c-reFr.dat", "c-imFr.dat")
    I = readspectral("c-reI.dat", "c-imI.dat")

    # read the constant Hartree part
    sigmaH = parse(Float64, readchomp(`extractcolumn $expectdir/custom.avg SigmaHd`))

    write("spectral/sigmaH.dat", readchomp(`extractcolumn $expectdir/custom.avg SigmaHd`))

    # Set new self-energy!
    solver.Σ = sigmaH .+ I .- (Fl .^ 2) ./ G .+ solver.μ

    open("spectral/imsigma.dat", "w") do io
        writedlm(io, [real.(solver.ωs) imag.(solver.Σ)])
    end
    
    open("spectral/resigma.dat", "w") do io
        writedlm(io, [real.(solver.ωs) real.(solver.Σ)])
    end

    @info "Finished post-processing of NRG run!"
    # NRG finished move one directory up!
    cd("..")
    nothing
end