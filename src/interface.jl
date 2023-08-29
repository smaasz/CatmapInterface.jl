@kwdef struct TState
    symbol::String
    beta::Float64
end

@kwdef struct ParsedReaction
    educts::Vector{Pair{String, Int}}
    products::Vector{Pair{String, Int}}
    tstate::Union{Nothing, TState}
end

abstract type AbstractSpecies end

struct FictiousSpecies <:AbstractSpecies
    species_name::String
    formation_energy::Float64
    pressure::Float64
    function FictiousSpecies(; species_name, formation_energy, pressure)
        if pressure < 0.0
            throw(DomainError("pressure must be nonnegative"))
        end
        new(species_name, formation_energy, pressure)
    end
end

struct GasSpecies <: AbstractSpecies
    species_name::String
    formation_energy::Float64
    pressure::Float64
    frequencies::Vector{Float64}
    function GasSpecies(; species_name, formation_energy, pressure, frequencies)
        if pressure < 0.0
            throw(DomainError("pressure must be nonnegative"))
        end
        if any(frequencies .<= 0.0)
            throw(DomainError("all frequencies must be positive"))
        end
        new(species_name, formation_energy, pressure, frequencies)
    end
end

struct AdsorbateSpecies <: AbstractSpecies
    species_name::String
    formation_energy::Float64
    coverage::Float64
    site::String
    surface_name::String
    frequencies::Vector{Float64}
    sigma_params::@NamedTuple{a::Float64, b::Float64}
    function AdsorbateSpecies(; species_name, formation_energy, coverage, site, surface_name, frequencies, sigma_params::@NamedTuple{a::Float64, b::Float64})
        if coverage < 0.0 || coverage > 1.0
            throw(DomainError("coverage must be between 0 and 1"))
        end
        if any(frequencies .<= 0.0)
            throw(DomainError("all frequencies must be positive"))
        end
        new(species_name, formation_energy, coverage, site, surface_name, frequencies, sigma_params)
    end
end
#AdsorbateSpecies(; formation_energy, coverage, site, surface_name, frequencies, sigma_params::Vector{Float64}) = AdsorbateSpecies(; formation_energy, coverage, site, surface_name, frequencies, sigma_params=(; a=sigma_params[1], b=sigma_params[2]))

struct TStateSpecies <: AbstractSpecies
    species_name::String
    formation_energy::Float64
    coverage::Float64
    site::String
    surface_name::String
    frequencies::Vector{Float64}
    sigma_params::@NamedTuple{a::Float64, b::Float64}
    β::Float64
    function TStateSpecies(; species_name, formation_energy, coverage, site, surface_name, frequencies, sigma_params::@NamedTuple{a::Float64, b::Float64}, β)
        if coverage < 0.0 || coverage > 1.0
            throw(DomainError("coverage must be between 0 and 1"))
        end
        if any(frequencies .<= 0.0)
            throw(DomainError("all frequencies must be positive"))
        end
        new(species_name, formation_energy, coverage, site, surface_name, frequencies, sigma_params, β)
    end
end
#TStateSpecies(; formation_energy, coverage, site, surface_name, frequencies, sigma_params::Vector{Float64}) = TStateSpecies(; formation_energy, coverage, site, surface_name, frequencies, sigma_params=(; a=sigma_params[1], b=sigma_params[2]))

struct SiteSpecies <: AbstractSpecies
    energy::Float64
    site_name::String
end

"""
Main data structure for interfacing CatMAP
"""
struct CatmapParams
    reactions::Vector{ParsedReaction}
    prefactors::Vector{Float64}
    species_list::Dict{String, AbstractSpecies}
    gas_thermo_mode::Symbol
    adsorbate_thermo_mode::Symbol
    electrochemical_thermo_mode::Symbol
    bulk_pH::Float64
    Upzc::Float64
    potential_reference_scale::String
    T::Float64
    function CatmapParams(; reactions, prefactors, species_list, gas_thermo_mode, adsorbate_thermo_mode, electrochemical_thermo_mode, bulk_pH, Upzc, potential_reference_scale, T)       
        if !(length(prefactors) == length(reactions))
            throw(ArgumentError("The number of prefactors must match the number of reactions"))
        end
        if any(prefactors .<= 0.0)
            throw(DomainError("all prefectors must be positive"))
        end
        # check that for each species included in the reactions it is contained in the species list
        species_set = Set{String}()
        for (; educts, products) in reactions
            union!(species_set, first.(educts))
            union!(species_set, first.(products))
        end
        for sp in species_set
            if !(haskey(species_list, sp))
                throw(ArgumentError("$sp not specified in the species_list"))
            end
        end
        # check that the provided modes are defined
        for mode in [gas_thermo_mode, adsorbate_thermo_mode, electrochemical_thermo_mode]
            if !isdefined(@__MODULE__, mode)
                throw(ArgumentError("$(String(:mode))=$mode is not implemented"))
            end
        end
        # check that reference scale is either RHE or SHE
        if !(potential_reference_scale == "RHE" || potential_reference_scale == "SHE")
            throw(ArgumentError("$potential_reference_scale must be either SHE or RHE"))
        end
        if T < 0.0
            throw(ArgumentError("temperature T=$T must be positive"))
        end
        new(reactions, prefactors, species_list, gas_thermo_mode, adsorbate_thermo_mode, electrochemical_thermo_mode, bulk_pH, Upzc, potential_reference_scale, T)
    end
end

function remove_whitespaces(s)
    s_vec = collect(s)
    join(s_vec[map(!isspace, s_vec)])
end

function parse_reactant_sum(rs)
    reactant_sum = r"^(?:[^+]+)(?:\+[^+]+)*$"
    if !isnothing(match(reactant_sum, rs))
        rts = split(rs, "+")
    else
        throw(ArgumentError("$rs is not a valid sum of reactants"))
    end
    reactant = r"^(?<factor>[1-9][0-9]*)?(?<species>[A-Za-z0-9]*)\*?_(?<site>[a-z])$"
    reactants = Pair{String, Int}[]
    for rt in rts
        match_reactant = match(reactant, rt)
        if !isnothing(match_reactant)
            factor = isnothing(match_reactant[:factor]) ? 1 : parse(Int, match_reactant[:factor])
            symbol = "$(match_reactant[:species])_$(match_reactant[:site])"
            push!(reactants, symbol => factor)
        else
            throw(ArgumentError("$rt is not a valid reactant"))
        end
    end
    reactants
end

function parse_reaction(r)
    r = remove_whitespaces(r)
    rxn = r"^(?<educts>[^<>]+)<->(?<products>[^<>]+)$"
    match_rxn = match(rxn, r)
    rxn_with_TS = r"^(?<educts>[^<>]+)<->(?<tstate>[^<>]+)<->(?<products>[^<>]+);beta=(?<beta>[0-9.]+)$"
    match_rxn_with_TS = match(rxn_with_TS, r)
    
    educts      = []
    products    = []
    tstate      = nothing
    if !isnothing(match_rxn)
        educts      = parse_reactant_sum(match_rxn[:educts])
        products    = parse_reactant_sum(match_rxn[:products]) 
    elseif !isnothing(match_rxn_with_TS)
        educts      = parse_reactant_sum(match_rxn_with_TS[:educts])
        products    = parse_reactant_sum(match_rxn_with_TS[:products])
        beta = 0.0
        try
            beta = parse(Float64, match_rxn_with_TS[:beta])    
        catch e
            throw(ArgumentError("$(match_rxn_with_TS[:beta]) is not a valid float"))
        end
        
        tstate = TState(symbol=match_rxn_with_TS[:tstate], beta=beta)
    else
        throw(ArgumentError("$r is not a valid reaction equation"))
    end
    return ParsedReaction(educts, products, tstate)
end

function parse_vector(s, Tval)::Vector{Tval}
    if isnothing(match(r"^\[.*?\]$", s))
        throw(ArgumentError("$s is not a valid vector"))        
    end
    entries = strip.(split(s[2:end-1], ","))
    entries = entries[(!isempty).(entries)]
    if Tval <: Number
        map(x->parse(Tval, x), entries)
    else
        entries
    end
end

function parse_energy_table(input_file_path)
    @local_unitfactors eV cm
    @local_phconstants h c_0
    (dc, dh) = readdlm(input_file_path, '\t', String; header=true)
    entry_types = [
        :surface_name      => String, 
        :site_name         => String,
        :species_name      => String,
        :formation_energy  => Float64, #in eV
        :bulk_structure    => String,
        :frequencies       => Vector{Float64}, # in (h / cm / eV * c_0)^-1
        :other_parameters  => Vector{String},
        :reference         => String,
    ]
    if !(dh[1,:] == string.(first.(entry_types)))
        throw(ArgumentError("headers are not valid"))
    end


    TRow = NamedTuple{(first.(entry_types)...,), Tuple{last.(entry_types)...}}
    table = TRow[]
    for (nrow, row_strings) in enumerate(eachrow(dc))
        row = []
        for ((entry_header, entry_type), entry) in zip(entry_types, row_strings)
            try
                if entry_type == String
                    push!(row, entry_header => entry)
                elseif entry_header == :formation_energy
                    push!(row, entry_header => parse(entry_type, entry) * eV)
                elseif entry_header == :frequencies
                    push!(row, entry_header => parse_vector(entry, eltype(entry_type)) .* h / cm / eV * c_0)
                elseif entry_type <: Number
                    push!(row, entry_header => parse(entry_type, entry))
                elseif entry_type <: Vector
                    push!(row, entry_header => parse_vector(entry, eltype(entry_type)))
                end
            catch e
                if isa(e, ArgumentError)
                    throw(ArgumentError("in row $nrow: $(e.msg)"))
                else
                    throw(e)
                end
            end
        end
        push!(table, (; row...))  
    end
    table
end

function parse_catmap_input(input_file_path)
    @pyinclude(input_file_path)
    
    reactions = ParsedReaction[]
    for r in py"rxn_expressions"
        push!(reactions, parse_reaction(r))
    end

    prefactors = py"prefactor_list"

    species_definitions = py"species_definitions"

    energies_file_path  = joinpath(dirname(input_file_path), py"input_file")
    energy_table        = parse_energy_table(energies_file_path)

    surface_name = let 
        surface_names = py"surface_names"
        if length(surface_names) > 1
            throw(ArgumentError("implementation only works for one surface for now"))
        end
        surface_names[1]
    end
    
    species_list = specieslist(reactions, species_definitions, energy_table, surface_name)


    bulk_pH = py"bulk_ph"
    
    Upzc                        = py"Upzc"
    potential_reference_scale   = py"potential_reference_scale"

    
    gas_thermo_mode             = Symbol(py"gas_thermo_mode")
    adsorbate_thermo_mode       = Symbol(py"adsorbate_thermo_mode")
    electrochemical_thermo_mode = Symbol(py"electrochemical_thermo_mode")

    T = 298

    CatmapParams(; reactions, prefactors, species_list, gas_thermo_mode, adsorbate_thermo_mode, electrochemical_thermo_mode, bulk_pH, Upzc, potential_reference_scale, T)
end


# for each species included in the reactions
#   1. check that there exists an entry in species_defs
#   2. determine the type of the species
#   3a. if type == gas, complete species_def by using energy table
#   3b. if type == adsorbate, then determine the site, for each surface complete species_def by using energy table
#   3c. if type == tstate, then do same as for adsorbate
#   3d. if type == fictious, then ...
function specieslist(reactions::Vector{ParsedReaction}, species_defs, energy_table, surface_name)
    # collect all species in a set
    species = Set{String}()
    for (; educts, products, tstate) in reactions
        union!(species, first.(educts))
        union!(species, first.(products))
    end
    specieslist = Dict{String, AbstractSpecies}()
    for s in species
        match_fictious  = match(r"^(?<species_name>ele|OH)_g$", s)
        match_gas       = match(r"^(?<species_name>[A-Za-z0-9]+)_g$", s)
        match_adsorbate = match(r"^(?<species_name>[A-Za-z0-9]+)_(?<site>[^g])$", s)
        match_site      = match(r"^_(?<site>[^g])$", s)
        if !isnothing(match_fictious)
            species_name            = match_fictious[:species_name]
            (; pressure)            = findspecies(species_name, "g", species_defs)
            (; formation_energy)    = findspecies(species_name, energy_table)
            pressure                = species_defs[s]["pressure"] 
            specieslist[s]          = FictiousSpecies(; species_name, formation_energy, pressure)
        elseif !isnothing(match_gas)
            species_name                        = match_gas[:species_name]
            (; pressure)                        = findspecies(species_name, "g", species_defs)
            (; formation_energy, frequencies)   = findspecies(species_name, energy_table)
            specieslist[s]                      = GasSpecies(; species_name, formation_energy, pressure, frequencies)
        elseif !isnothing(match_adsorbate)
            species_name                        = match_adsorbate[:species_name]
            site                                = match_adsorbate[:site]
            (; sigma_params)                    = findspecies(species_name, site, species_defs)
            sigma_params                        = (; a = sigma_params[1], b = sigma_params[2])
            coverage                            = 0.0
            (; site_names)                      = findspecies("", site, species_defs)
            site_name                           = site_names[1]
            (; formation_energy, frequencies)   = findspecies(species_name, energy_table; surface_name, site_name)
            specieslist[s]                      = AdsorbateSpecies(; species_name, formation_energy, coverage, site, surface_name, frequencies, sigma_params)
        elseif !isnothing(match_site)
            site            = match_site[:site]
            (; site_names)  = findspecies("", site, species_defs)
            specieslist[s]  = SiteSpecies(0.0, site_names[1])
        else
            throw(ArgumentError("species $s is not a valid ficitious gas, gas, adsorbate, or site"))
        end
    end
    for (; tstate) in reactions
        if isnothing(tstate)
            continue
        end
        match_tstate = match(r"^(?<species_name>[A-Za-z0-9\-]+)_(?<site>[^g])$", tstate.symbol)
        if !isnothing(match_tstate)
            species_name                        = match_tstate[:species_name]
            site                                = match_tstate[:site]
            (; sigma_params)                    = findspecies(species_name, site, species_defs)
            sigma_params                        = (; a = sigma_params[1], b = sigma_params[2])
            coverage                            = 0.0
            (; site_names)                      = findspecies("", site, species_defs)
            site_name                           = site_names[1]
            (; formation_energy, frequencies)   = findspecies(species_name, energy_table; surface_name, site_name)
            specieslist[tstate.symbol]          = TStateSpecies(; species_name, formation_energy, coverage, site, surface_name, frequencies, sigma_params, β=tstate.beta)
        else
            throw(ArgumentError("$(tstate.symbol) is not a valid transition state"))
        end
    end
    # special case to make sure that H2_g and H2O_g are contained in the list
    if !haskey(specieslist, "H2_g")
        pressure = 1.0
        local formation_energy, frequencies
        try
            (; formation_energy, frequencies) = findspecies("H2", energy_table)
        catch e
            if !isa(e, ArgumentError)
                throw(e)
            end
        else
            specieslist["H2_g"] = GasSpecies(; species_name="H2", formation_energy, pressure, frequencies) 
        end
    end
    if !haskey(specieslist, "H2O_g")
        pressure = 1.0
        local formation_energy, frequencies
        try
            (; formation_energy, frequencies) = findspecies("H2O", energy_table)
        catch e
            if !isa(e, ArgumentError)
                throw(e)
            end
        else
            specieslist["H2O_g"] = GasSpecies(; species_name="H2O", formation_energy, pressure, frequencies) 
        end
    end
    @assert haskey(specieslist, "H2_g") && haskey(specieslist, "H2O_g")
    specieslist
end

function findspecies(species_name::AbstractString, site::AbstractString, species_defs)
    species = isempty(species_name) ? site : "$(species_name)_$(site)"
    try
        (; [Symbol(k) => v for (k,v) in species_defs[species]]...)
    catch e
        if isa(e, KeyError)
            throw(ArgumentError("species $(species_name)_$(site) is not specified in species_definitions"))
        else
            throw(e)
        end
    end
end

function findspecies(species_name::AbstractString, energy_table; surface_name="None", site_name="gas")
    found = nothing
    for row in energy_table
        if (row.surface_name == surface_name) && (row.site_name == site_name) && (row.species_name == species_name)
            found = row
        end
    end
    if isnothing(found)
        throw(ArgumentError("species $(species_name) is not specified in the energy table"))
    end
    found
end