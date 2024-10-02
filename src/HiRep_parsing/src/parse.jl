function gaugegroup(file)
    if HDF5.ishdf5(file)
        hdf5 = h5open(file, "r")
        return read(hdf5, "gauge group")
    else
        return gaugegroup_log(file)
    end
end
function quarkmasses(file)
    if HDF5.ishdf5(file)
        hdf5 = h5open(file, "r")
        return read(hdf5, "quarkmasses")
    else
        return quarkmasses_log(file)
    end
end
function quarkmasses_chimera(file)
    if HDF5.ishdf5(file)
        hdf5 = h5open(file, "r")
        mf = read(hdf5, "quarkmasses_fundamental")
        mas = read(hdf5, "quarkmasses_antisymmetric")
        return mf, mas
    else
        return quarkmasses_chimera_log(file)
    end
end
function APE_smearing(file)
    if HDF5.ishdf5(file)
        hdf5 = h5open(file, "r")
        APE_eps = read(hdf5, "APE_eps")
        APE_level = read(hdf5, "APE_levels")
        return APE_eps, APE_level
    else
        return APE_smearing_logfile(file)
    end
end
function Wuppertal_smearing_mixed(file)
    if HDF5.ishdf5(file)
        hdf5 = h5open(file, "r")
        antisymmetric_eps = read(hdf5, "Wuppertal_eps_anti")
        fundamental_eps = read(hdf5, "Wuppertal_eps_fund")
        return antisymmetric_eps, fundamental_eps
    else
        return Wuppertal_smearing_mixed_logfile(file)
    end
end
function latticesize(file)
    if HDF5.ishdf5(file)
        hdf5 = h5open(file, "r")
        return read(hdf5, "lattice")
    else
        return latticesize_log(file)
    end
end
function plaquettes(file)
    if HDF5.ishdf5(file)
        hdf5 = h5open(file, "r")
        return read(hdf5, "plaquette")
    else
        return plaquettes_log(file)
    end
end
function inverse_coupling(file)
    if HDF5.ishdf5(file)
        hdf5 = h5open(file, "r")
        return read(hdf5, "beta")
    else
        return inverse_coupling_log(file)
    end
end
function correlators(file, type, key; kws...)
    if HDF5.ishdf5(file)
        hdf5 = h5open(file, "r")
        return read(hdf5, key)
    else
        return correlators_logfile(file, type, key; kws...)
    end
end
function correlators(file, type, key, nhits; kws...)
    if HDF5.ishdf5(file)
        hdf5 = h5open(file, "r")
        corr = read(hdf5, key)
        return corr
    else
        d = correlators_logfile(file, type, key; masses, mass, kws...)
        return reshape_disconnected(d, nhits; rescale = 1)
    end
end
function reshape_disconnected(d::Matrix{S}, nhits; rescale = 1) where {S}
    # the rounding automatically removes incomplete calcluations
    T, N = size(d)
    n = div(N, nhits)
    m = zeros(eltype(d[1]), (n, nhits, T))
    for i = 1:n, j = 1:nhits, t = 1:T
        m[i, j, t] = d[t, (i-1)*nhits+j]
    end
    @. m = rescale * m
    return m
end
function gaugegroup_log(file)
    for line in eachline(file)
        if occursin("Gauge group", line)
            pos = findlast(' ', line)
            return strip(line[pos:end])
        end
    end
end
function quarkmasses_log(file; pattern = "[MAIN][0]Mass[0]")
    masses = Float64[]
    for line in eachline(file)
        if occursin(pattern, line)
            s = split(line, (","))
            for i in eachindex(s)
                m = parse(Float64, split(s[i], "=")[2])
                append!(masses, m)
            end
            return masses
        end
    end
end
function quarkmasses_chimera_log(file)
    mf = quarkmasses_log(file; pattern = "[MAIN][0]mf[0]")
    mas = quarkmasses_log(file; pattern = "[MAIN][0]mas[0]")
    return mf, mas
end
function latticesize_log(file)
    for line in eachline(file)
        if occursin("Global size is", line)
            pos = last(findfirst("Global size is", line)) + 1
            sizestring = lstrip(line[pos:end])
            latticesize = parse.(Int, split(sizestring, "x"))
            return latticesize
        end
    end
end
function plaquettes_log(file)
    plaquettes = Float64[]
    for line in eachline(file)
        if occursin("Plaquette", line)
            line = replace(line, "=" => " ")
            line = replace(line, ":" => " ")
            p = parse(Float64, split(line)[end])
            append!(plaquettes, p)
        end
    end
    return plaquettes
end
function _match_config_name(filename)
    regex =
        r".*/(?<run>[^/]*)_(?<T>[0-9]+)x(?<L>[0-9]+)x[0-9]+x[0-9]+nc[0-9]+(?:r[A-Z]+)?(?:nf[0-9]+)?b(?<beta>[0-9]+\.[0-9]+)?(?:m-?[0-9]+\.[0-9]+)?n(?<conf>[0-9]+)"
    return match(regex, filename)
end
function inverse_coupling_log(file)
    try
        l = split(file, "beta")[end]
        β = parse(Float64, split(l, "m")[1])
        return β
    catch
        for line in eachline(file)
            if occursin("Configuration from", line)
                match = _match_config_name(line)
                β = parse(Float64, match[:beta])
                return β
            end
        end
    end
end
function APE_smearing_logfile(file)
    APE_level = Int64[]
    APE_eps = Float64[]
    # start with reference level for no smearing
    N = -1
    for line in eachline(file)
        if startswith(line, "[APE][0]APE smearing with val")
            eps = parse(Float64, last(split(line, "=")))
            append!(APE_eps, eps)
        end
        if startswith(line, "[APE][0]N=")
            pos1 = length("[APE][0]N=") + 1
            pos2 = findfirst('<', line) - 1
            N0 = parse(Int, line[pos1:pos2])
            N = max(N, N0)
        end
        if startswith(line, "[APE][0]APE smearing END")
            append!(APE_level, N)
            N = -1 # reset to reference value for no smearing
        end
    end
    return APE_eps, APE_level
end
function Wuppertal_smearing_mixed_logfile(file)
    fundamental_eps = Float64[]
    antisymmetric_eps = Float64[]

    #default reference values for no smearing
    epsA = 0.0
    epsF = 0.0

    for line in eachline(file)
        if startswith(line, "[SMEAR][0]source smearing epsilon =")
            pos1 = length("[SMEAR][0]source smearing epsilon =") + 1
            pos2 = first(findfirst("iterations:", line)) - 1
            epsA = parse(Float64, line[pos1:pos2])
        end
        if startswith(line, "[SMEAR][0]Fundamental source smearing epsilon =")
            pos1 = length("[SMEAR][0]Fundamental source smearing epsilon =") + 1
            pos2 = first(findfirst("iterations:", line)) - 1
            epsF = parse(Float64, line[pos1:pos2])
        end
        if occursin("analysed", line)
            append!(antisymmetric_eps, epsA)
            append!(fundamental_eps, epsF)
        end
    end
    return antisymmetric_eps, fundamental_eps
end
function correlators_logfile(file, type, key; kws...)
    corrs = parse_spectrum(file, type; filterkey = true, key_pattern = key, kws...)
    return reduce(hcat, getindex.(corrs, key))
end
function parse_spectrum(
    file,
    type;
    disconnected = false,
    masses = false,
    mass = "",
    filterkey = false,
    key_pattern = "",
    nhits = 1,
    with_progress = false,
)
    T = latticesize(file)[1]
    corr = zeros(T) # preallocate array for parsing of correlator
    dict = Dict{String,Vector{Float64}}()
    dictarray = Dict{String,Vector{Float64}}[]
    conf0 = 0
    src0 = 0
    # keep track of position in file for progress meter
    with_progress &&
        (p = Progress(countlines(file); dt = 1, desc = "Match $type: Progress:"))
    for line in eachline(file)
        with_progress && next!(p)
        if occursin(type, line)
            if masses
                occursin("mass=$mass", line) || continue
            end
            if filterkey
                occursin(key_pattern, line) || continue
            end
            # get configuration number
            pos_num = findfirst('#', line)
            end_num = findnext(' ', line, pos_num)
            conf = parse(Int, line[pos_num+1:end_num-1])
            # find number of the source if available
            if disconnected
                pos_src = last(findfirst("src", line)) + 1
                end_src = findnext(' ', line, pos_src + 1)
                src = parse(Int, line[pos_src:end_src])
            else
                src = 0
            end
            # find last '=' sign which separates values from Γ structure
            # TODO this does not work for momenta
            pos_eq = findlast('=', line)
            #key_st = findprev(' ',line,pos_eq)
            key_st = last(findfirst(type, line)) + 1
            key = line[key_st+1:pos_eq-1]
            if disconnected
                # create new entry if configuration or source number changes
                # if we need to parse more than one source at a time per configuration
                if conf0 != conf || src0 != src
                    if !isempty(dict)
                        push!(dictarray, dict)
                        dict = Dict{String,Vector{Float64}}()
                    end
                end
            end
            # parse corrrelator values
            pos_0 = findnext(' ', line, pos_eq)
            for t = 1:T
                pos_1 = findnext(' ', line, pos_0 + 1)
                corr[t] = Parsers.parse(Float64, line[pos_0:pos_1])
                pos_0 = pos_1
            end
            dict[key] = copy(corr)
            conf0 = conf
            src0 = src
        end
        if !disconnected
            # If we only have one source at a time and possibly one configuration
            # at a time: the method used to separate distinct measurements fails.
            # In this case the end of measurement on a given confiuration is
            # signalled by a line that reads:
            # [MAIN][0]Configuration #N: analysed in [a sec b usec]
            if occursin("analysed", line)
                if !isempty(dict)
                    push!(dictarray, dict)
                    dict = Dict{String,Vector{Float64}}()
                end
            end
        end
    end
    if !isempty(dict)
        push!(dictarray, dict)
    end
    return _reshape_connected(dictarray; disconnected, nhits)
end
function _reshape_connected(dict; disconnected = false, nhits = 1)
    corrs = Dict{String,Array{Float64}}()
    for k in keys(dict[1])
        if disconnected
            corrs[k] = reshape_disconnected(reduce(hcat, getindex.(dict, k)), nhits)
            # rename keys for disconnected measurements
            corrs = _rename_disconnected_keys(corrs)
        else
            corrs[k] = permutedims(reduce(hcat, getindex.(dict, k)))
        end
    end
    return corrs
end
function _rename_disconnected_keys(corrs_discon)
    new_corrs_discon = Dict{String,Array{Float64}}()
    for k in keys(corrs_discon)
        val = corrs_discon[k]
        k_new = replace(k, "_re" => "", "_disc" => "")
        new_corrs_discon[k_new] = val
    end
    return new_corrs_discon
end
#####################################################
# Parsing using regular expressions (for smearing)  #
#####################################################
function parse_spectrum_with_regexp(
    file,
    type;
    disconnected = false,
    masses = false,
    mass = "",
    filterkey = false,
    key_pattern = "",
    nhits = 1,
    with_progress = false,
)
    T = latticesize(file)[1]
    corr = zeros(T) # preallocate array for parsing of correlator
    dict = Dict{String,Vector{Float64}}()
    dictarray = Dict{String,Vector{Float64}}[]
    conf0 = 0
    src0 = 0
    # keep track of position in file for progress meter
    if with_progress
        p = Progress(countlines(file); dt = 1, desc = "Match $type: Progress:")
        linecount = 0
    end
    for line in eachline(file)
        if with_progress
            linecount += 1
            # don't update the progress meter after every line
            linecount % 100 == 0 && update!(p, linecount)
        end
        if occursin(type, line)
            m = match(type, line)
            matched_type = m.match
            if masses
                occursin("mass=$mass", line) || continue
            end
            if filterkey
                occursin(key_pattern, line) || continue
            end
            # get configuration number
            pos_num = findfirst('#', line)
            end_num = findnext(' ', line, pos_num)
            conf = parse(Int, line[pos_num+1:end_num-1])
            # find number of the source if available
            if disconnected
                pos_src = last(findfirst("src", line)) + 1
                end_src = findnext(' ', line, pos_src + 1)
                src = parse(Int, line[pos_src:end_src])
            else
                src = 0
            end
            # find last '=' sign which separates values from Γ structure
            # TODO this does not work for momenta
            pos_eq = findlast('=', line)
            key_st = last(findfirst(matched_type, line)) + 1
            key = joinpath(matched_type, line[key_st+1:pos_eq-1])
            if disconnected
                # create new entry if configuration or source number changes
                # if we need to parse more than one source at a time per configuration
                if conf0 != conf || src0 != src
                    if !isempty(dict)
                        push!(dictarray, dict)
                        dict = Dict{String,Vector{Float64}}()
                    end
                end
            end
            # parse corrrelator values
            pos_0 = findnext(' ', line, pos_eq)
            for t = 1:T
                pos_1 = findnext(' ', line, pos_0 + 1)
                corr[t] = Parsers.parse(Float64, line[pos_0:pos_1])
                pos_0 = pos_1
            end
            dict[key] = copy(corr)
            conf0 = conf
            src0 = src
        end
        if !disconnected
            # If we only have one source at a time and possibly one configuration
            # at a time: the method used to separate distinct measurements fails.
            # In this case the end of measurement on a given confiuration is
            # signaled by a line that reads:
            # [MAIN][0]Configuration #N: analysed in [a sec b usec]
            if occursin("analysed", line)
                if !isempty(dict)
                    push!(dictarray, dict)
                    dict = Dict{String,Vector{Float64}}()
                end
            end
        end
    end
    with_progress && finish!(p)
    if !isempty(dict)
        push!(dictarray, dict)
    end
    return _reshape_connected(dictarray; disconnected, nhits)
end
#################################################
# Disconnected Measurements from /Disocnnected  #
#################################################
function dilution(file)
    for line in eachline(file)
        if occursin("will be used", line)
            eo = occursin("eo", lowercase(line))
            time = occursin("time", lowercase(line))
            spin = occursin("spin", lowercase(line))
            color = occursin("color", lowercase(line))
            return eo, time, spin, color
        end
    end
end
function ncolors(file)
    for line in eachline(file)
        if occursin("Gauge group", line)
            pos1 = findfirst('(', line) + 1
            pos2 = findfirst(')', line) - 1
            colors = parse(Int, line[pos1:pos2])
            return colors
        end
    end
end
function hits(file)
    for line in eachline(file)
        if occursin("Number of noise vector : nhits", line)
            pos = findfirst('=', line)
            nhits = parse(Int, line[pos+1:end])
            return nhits
        end
    end
end
function nconfigs(file)
    nconf = 0
    for line in eachline(file)
        if occursin("read", line)
            nconf += 1
        end
    end
    return nconf
end
function parse_disconnected(file)
    nGamma = 16
    T = latticesize(file)[1]
    cut = length("[CORR][0]")
    nconf = nconfigs(file)
    nhits = hits(file)
    # Query whether even-odd dilution, time dilution, spin dilution
    # or color dilution has been used in the measurement
    eo, time, spin, color = dilution(file)
    # True of no dilution has been used (pure volume source)
    pure = (!eo * !time * !spin * !color)
    # Two indices are needed if even-odd dilution has been used
    neo = eo ? 2 : 1
    # Nc indices are used in color dilution has been used
    ncolor = color ? ncolors(file) : 1
    data = zeros(nGamma, nconf, nhits, ncolor, neo, T)
    conf = 0
    col = 1
    eoi = 1
    p = Progress(countlines(file), 1)
    for line in eachline(file)
        if startswith(line, "[CORR][0]")
            # remove lines containing the statements
            # "Contraction done" and "Start to perform contractions"
            # "Number of noise vector"
            if isdigit(line[cut+1])
                pos1 = findnext(' ', line, cut)
                pos2 = findnext(' ', line, pos1 + 1)
                Γ = Parsers.parse(Int, line[pos1:pos2]) + 1
                t = Parsers.parse(Int, line[cut+1:pos1]) + 1
                pos3 = findnext(' ', line, pos2 + 1)
                hit = Parsers.parse(Int, line[pos2:pos3]) + 1
                if color
                    pos2 = findnext(' ', line, pos3 + 1)
                    pos2, pos3 = pos3, pos2
                    col = Parsers.parse(Int, line[pos2:pos3]) + 1
                end
                if eo
                    pos2 = findnext(' ', line, pos3 + 1)
                    pos2, pos3 = pos3, pos2
                    eoi = Parsers.parse(Int, line[pos2:pos3]) + 1
                end
                if pure
                    for _ = 1:3
                        pos3 = findnext(' ', line, pos3 + 1)
                    end
                end
                pos4 = findnext(' ', line, pos3 + 1)
                Re = Parsers.parse(Float64, line[pos3:pos4])
                if !(hit > nhits)
                    data[Γ, conf, hit, col, eoi, t] = Re
                end
            end
        end
        if startswith(line, "[MAIN][0]Configuration from")
            conf += 1
        end
        next!(p)
    end
    return average_dilution(data)
end
function average_dilution(d::Array{S,6}) where {S}
    nGamma, n, nhits, ncol, neo, T = size(d)
    data = zeros(nGamma, n, nhits, T)
    # averages over different eo pattern and color dilution
    for Γ = 1:nGamma, t = 1:T, j = 1:nhits, eo = 1:neo, col = 1:ncol, i = 1:n
        data[Γ, i, j, t] += d[Γ, i, j, col, eo, t]
    end
    # put everything into a dictionary format to match the output of the other routines
    disc = Dict{String,Array{Float64}}()
    disc["g5"] = data[1, :, :, :]
    disc["g1"] = data[2, :, :, :]
    disc["g2"] = data[3, :, :, :]
    disc["g3"] = data[4, :, :, :]
    disc["-ig0g5"] = data[5, :, :, :]
    disc["-ig0g1"] = data[6, :, :, :]
    disc["-ig0g2"] = data[7, :, :, :]
    disc["-ig0g4"] = data[8, :, :, :]
    disc["id"] = data[9, :, :, :]
    disc["-ig5g1"] = data[10, :, :, :]
    disc["-ig5g2"] = data[11, :, :, :]
    disc["-ig5g3"] = data[12, :, :, :]
    disc["g0"] = data[13, :, :, :]
    disc["-ig5g0g1"] = data[14, :, :, :]
    disc["-ig5g0g2"] = data[15, :, :, :]
    disc["-ig5g0g3"] = data[16, :, :, :]
    return disc
end

function confignames(file)
    fns = AbstractString[]
    for line in eachline(file)
        if occursin("read", line)
            if occursin("Configuration", line)
                pos1 = findlast('/', line)
                pos2 = findnext(']', line, pos1)
                push!(fns, line[pos1+1:pos2-1])
            end
        end
    end
    return fns
end
