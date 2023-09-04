@kwdef struct TState
    symbol::String
    beta::Float64
end

@kwdef struct ParsedReaction
    educts::Vector{Pair{String, Int}}
    products::Vector{Pair{String, Int}}
    tstate::Union{Nothing, TState}
end


"""
$(TYPEDEF)

Main data structure for interfacing CatMAP

$(TYPEDFIELDS)
"""
struct CatmapParams
    """
    List of parsed reactions
    """
    reactions::Vector{ParsedReaction}
    """
    Prefactors in the Arrhenius relation specifying the activation energy of the corresponding reactions
    """
    prefactors::Vector{Float64}
    """
    List of species including all reactants
    """
    species_list::Dict{String, AbstractSpecies}
    """
    Mode of the thermodynamical correction to the formation energy of the gases
    """
    gas_thermo_mode::Symbol
    """
    Mode of the thermodynamical correction to the formation energy of the adsorbates and transition states
    """
    adsorbate_thermo_mode::Symbol
    """
    Mode of the electrochemical correction to the formation energy of the species
    """
    electrochemical_thermo_mode::Symbol
    """
    pH value in the bulk of the electrolyte
    """
    bulk_pH::Float64
    """
    Potential of zero charge with respect to the specified potential reference scale
    """
    Upzc::Float64
    """
    Either RHE or SHE
    """
    potential_reference_scale::String
    """
    Temperature in the bulk
    """
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


const re_reactant_sum  = r"^(?:[^+]+)(?:\+[^+]+)*$"
const re_reactant      = r"^(?<factor>[1-9][0-9]*)?(?<species>[A-Za-z0-9]*)\*?_(?<site>[a-z])$"
"""
$(SIGNATURES) 

Parse a string as a sum of reactants into a list of pairs reactant/stoichiometric factor.
"""
function parse_reactant_sum(rs::AbstractString)
    if !isnothing(match(re_reactant_sum, rs))
        rts = split(rs, "+")
    else
        throw(ArgumentError("$rs is not a valid sum of reactants"))
    end
    reactants = Pair{String, Int}[]
    for rt in rts
        match_reactant = match(re_reactant, rt)
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


const re_rxn            = r"^(?<educts>[^<>]+)<->(?<products>[^<>]+)$"
const re_rxn_with_TS    = r"^(?<educts>[^<>]+)<->(?<tstate>[^<>]+)<->(?<products>[^<>]+);beta=(?<beta>[0-9.]+)$"
"""
$(SIGNATURES) 

Parse a specification of a chemical reaction into a [`ParsedReaction`](@ref).

# Example

```jldoctest
julia> CatmapInterface.parse_reaction("CO*_t <-> CO_g + *_t")
CatmapInterface.ParsedReaction(["CO_t" => 1], ["CO_g" => 1, "_t" => 1], nothing)

julia> CatmapInterface.parse_reaction("COOH*_t + H2O_g + ele_g <-> COOH-H2O-ele_t <-> CO*_t + H2O_g + OH_g + *_t; beta=0.5")
CatmapInterface.ParsedReaction(["COOH_t" => 1, "H2O_g" => 1, "ele_g" => 1], ["CO_t" => 1, "H2O_g" => 1, "OH_g" => 1, "_t" => 1], CatmapInterface.TState("COOH-H2O-ele_t", 0.5))
```
"""
function parse_reaction(r::AbstractString)
    r = remove_whitespaces(r)
    
    match_rxn           = match(re_rxn, r)
    match_rxn_with_TS   = match(re_rxn_with_TS, r)
    
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


"""
$(SIGNATURES) 

Parse an energy table into a list of named tuples.

The entries of an energy table must be separated by tabs.
The first row of an energy table includes the headers of the columns.
The following are recognized (the order of the columns does not matter): 
- surface_name
- site_name
- species_name
- formation_energy: in eV
- bulk_structure
- frequencies: list of enclosed with brackets and comma-separated, in cm⁻¹
- other_parameters
- reference
"""
function parse_energy_table(input_file_path)
    @local_unitfactors eV cm
    (dc, dh) = readdlm(input_file_path, '\t', String; header=true)
    entry_types = [
        :surface_name      => String, 
        :site_name         => String,
        :species_name      => String,
        :formation_energy  => Float64, #in J/mole
        :bulk_structure    => String,
        :frequencies       => Vector{Float64}, # in m⁻¹
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
                    push!(row, entry_header => parse_vector(entry, eltype(entry_type)) ./ cm)
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


"""
$(SIGNATURES) 

Parse the content of a CatMAP input file needed for a Poisson-Nernst-Planck model.

The input file is first interpreted by a Python interpreter and then parsed.
The following information is used and must be specified:
- rxn_expressions
- prefactor_list
- species_definitions (including ...)
- input_file (including the energy table)
- surface_names
- bulk_ph
- Upzc
- potential_reference_scale
- gas_thermo_mode
- adsorbate_thermo_mode
- electrochemical_thermo_mode
See [CatMAP documentation](https://catmap.readthedocs.io/en/latest/index.html) for details.
"""
function parse_catmap_input(input_file_path::AbstractString)
    @assert isfile(input_file_path)
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


const re_fictious_gas   = r"^(?<species_name>ele|OH)_g$"
const re_gas            = r"^(?<species_name>[A-Za-z0-9]+)_g$"
const re_adsorbate      = r"^(?<species_name>[A-Za-z0-9]+)_(?<site>[^g])$"
const re_site           = r"^_(?<site>[^g])$"
const re_tstate         = r"^(?<species_name>[A-Za-z0-9\-]+)_(?<site>[^g])$"

"""
$(SIGNATURES) 

Collect the specifications of all reactants in a list.
"""
function specieslist(reactions::Vector{ParsedReaction}, species_defs, energy_table, surface_name)
    @local_unitfactors μA cm
    # collect all species in a set
    species = Set{String}()
    for (; educts, products, tstate) in reactions
        union!(species, first.(educts))
        union!(species, first.(products))
    end
    species_list = Dict{String, AbstractSpecies}()
    for s in species
        match_fictious  = match(re_fictious_gas,s)
        match_gas       = match(re_gas,         s)
        match_adsorbate = match(re_adsorbate,   s)
        match_site      = match(re_site,        s)
        if !isnothing(match_fictious)
            species_name            = match_fictious[:species_name]
            (; pressure)            = findspecies(species_name, "g", species_defs)
            (; formation_energy)    = findspecies(species_name, energy_table)
            pressure                = species_defs[s]["pressure"] 
            species_list[s]         = FictiousSpecies(; species_name, formation_energy, pressure)
        elseif !isnothing(match_gas)
            species_name                        = match_gas[:species_name]
            (; pressure)                        = findspecies(species_name, "g", species_defs)
            (; formation_energy, frequencies)   = findspecies(species_name, energy_table)
            species_list[s]                     = GasSpecies(; species_name, formation_energy, pressure, frequencies)
        elseif !isnothing(match_adsorbate)
            species_name                        = match_adsorbate[:species_name]
            site                                = match_adsorbate[:site]
            (; sigma_params)                    = findspecies(species_name, site, species_defs)
            sigma_params                        = (; a = sigma_params[2] / (μA/cm^2), b = sigma_params[1] / (μA/cm^2)^2)
            coverage                            = 0.0
            (; site_names)                      = findspecies("", site, species_defs)
            site_name                           = site_names[1]
            (; formation_energy, frequencies)   = findspecies(species_name, energy_table; surface_name, site_name)
            species_list[s]                     = AdsorbateSpecies(; species_name, formation_energy, coverage, site, surface_name, frequencies, sigma_params)
        elseif !isnothing(match_site)
            site            = match_site[:site]
            (; site_names)  = findspecies("", site, species_defs)
            species_list[s]  = SiteSpecies(0.0, site_names[1])
        else
            throw(ArgumentError("species $s is not a valid ficitious gas, gas, adsorbate, or site"))
        end
    end
    for (; educts, products, tstate) in reactions
        if isnothing(tstate)
            continue
        end
        match_tstate = match(re_tstate, tstate.symbol)
        if !isnothing(match_tstate)
            species_name                        = match_tstate[:species_name]
            site                                = match_tstate[:site]
            (; sigma_params)                    = findspecies(species_name, site, species_defs)
            sigma_params                        = (; a = sigma_params[2] / (μA/cm^2), b = sigma_params[1] / (μA/cm^2)^2)
            coverage                            = 0.0
            (; site_names)                      = findspecies("", site, species_defs)
            site_name                           = site_names[1]
            (; formation_energy, frequencies)   = findspecies(species_name, energy_table; surface_name, site_name)
            species_list[tstate.symbol]         = TStateSpecies(; species_name, formation_energy, coverage, site, surface_name, frequencies, sigma_params, β=tstate.beta, between_species=[first.(educts); first.(products)])
        else
            throw(ArgumentError("$(tstate.symbol) is not a valid transition state"))
        end
    end
    # special case to make sure that H2_g and H2O_g are contained in the list
    if !haskey(species_list, "H2_g")
        pressure = 1.0
        local formation_energy, frequencies
        try
            (; formation_energy, frequencies) = findspecies("H2", energy_table)
        catch e
            if !isa(e, ArgumentError)
                throw(e)
            end
        else
            species_list["H2_g"] = GasSpecies(; species_name="H2", formation_energy, pressure, frequencies) 
        end
    end
    if !haskey(species_list, "H2O_g")
        pressure = 1.0
        local formation_energy, frequencies
        try
            (; formation_energy, frequencies) = findspecies("H2O", energy_table)
        catch e
            if !isa(e, ArgumentError)
                throw(e)
            end
        else
            species_list["H2O_g"] = GasSpecies(; species_name="H2O", formation_energy, pressure, frequencies) 
        end
    end
    @assert haskey(species_list, "H2_g") && haskey(species_list, "H2O_g")
    species_list
end

"""
$(SIGNATURES) 

Extract specification of a species at a given site from the species_defs of a CatMAP input file.
"""
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

"""
$(SIGNATURES) 

Extract specification of a species from an energy table.
"""
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