abstract type AbstractSpecies end

struct AdsorbateSpecies <: AbstractSpecies
    energy::Float64
    frequencies::Vector{Float64}
    pressure::Float64
    sigma_params::@NamedTuple{a::Float64, b::Float64}
    function Species(; pressure, sigma_params)
        if pressure < 0.0
            throw(DomainError("pressure must be nonnegative"))
        end
        if any(frequencies .<= 0.0)
            throw(DomainError("all frequencies must be positive"))
        end
        new(energy, frequencies, pressure, sigma_params)
    end
end

struct GasSpecies <: AbstractSpecies
    energy::Float64
    frequencies::Vector{Float64}
    pressure::Float64

end


"""
Main data structure for interfacing CatMAP
"""
struct CatmapParams
    rn::ReactionSystem
    prefactors::Vector{Float64}
    species_definitions::Dict{String, AbstractSpecies}
    gas_thermo_mode::Symbol
    adsorbate_thermo_mode::Symbol
    electrochemical_thermo_mode::Symbol
    function CatmapParams(; rn, prefactors, species, gas_thermo_mode, adsorbate_thermo_mode, electrochemical_thermo_mode)       
        if !(length(prefactors == numreactions(rn)))
            throw(ArgumentError("The number of prefactors must match the number of reactions"))
        end
        if any(prefactors .<= 0.0)
            throw(DomainError("all prefectors must be positive"))
        end
        # check that for each species in the reaction network parameters are provided
        function isspecified(sp)
            any(map(isequal(sp), keys(species)))
        end
        for sp in species(rn)
            if !isspecified(sp) 
                throw(ArgumentError("$sp not specified in the argument species_definitions"))
            end
        end
        # check that the provided modes are defined
        for mode in [gas_thermo_mode, adsorbate_thermo_mode, electrochemical_thermo_mode]
            if !isdefined(@__MODULE__, mode)
                throw(ArgumentError("$(String(:mode))=$mode is not implemented"))
            end
        end
        new(rn, prefactors, species_definitions, gas_thermo_mode, adsorbate_thermo_mode, electrochemical_thermo_mode)
    end
end

function isspecified(sp::SymbolicUtils.BasicSymbolic{Real}, params::CatmapParams, )
    sp_name = string(operation(sp))
    any(map(isequal(spname), keys(params.species_definitions)))
end

function compile_rn(rn::ReactionSystem, params::CatmapParams)
    @assert all(map(isspecified(params), species(rn))) "all species used in the reaction network must be specified"
    

end




@kwdef struct TState
    symbol::String
    beta::Float64
end

@kwdef struct ParsedReaction
    educts::Vector{Pair{String, Int}}
    products::Vector{Pair{String, Int}}
    tstate::Union{Nothing, TState}
end


function parse_reactant_sum(rs)
    reactant_sum = r"^(?:[^+]+)(?:\+[^+]+)*$"
    if !isnothing(match(reactant_sum, rs))
        rts = split(rs, "+")
    else
        throw(ArgumentError("$rs is not a valid sum of reactants"))
    end
    reactant = r"^(?<factor>[1-9][0-9]*)?(?<symbol>[A-Za-z0-9]*\*?_[a-z])$"
    reactants = Pair{String, Int}[]
    for rt in rts
        match_reactant = match(reactant, rt)
        if !isnothing(match_reactant)
            factor = isnothing(match_reactant[:factor]) ? 1 : parse(Int, match_reactant[:factor])
            push!(reactants, match_reactant[:symbol] => factor)
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
    (dc, dh) = readdlm(input_file_path, '\t', String; header=true)
    entry_types = [
        "surface_name"      => String, 
        "site_name"         => String,
        "species_name"      => String,
        "formation_energy"  => Float64,
        "bulk_structure"    => String,
        "frequencies"       => Vector{Float64},
        "other_parameters"  => Vector{String},
        "reference"         => String,
    ]
    if !(dh[1,:] == first.(entry_types))
        throw(ArgumentError("headers are not valid"))
    end

    table = []
    for (nrow, row_strings) in enumerate(eachrow(dc))
        row = []
        for ((entry_header, entry_type), entry) in zip(entry_types, row_strings)
            try
                if entry_type == String
                    push!(row, entry)
                elseif entry_type <: Number
                    push!(row, parse(entry_type, entry))
                elseif entry_type <: Vector
                    push!(row, parse_vector(entry, eltype(entry_type)))
                end
            catch e
                if isa(e, ArgumentError)
                    throw(ArgumentError("in row $nrow: $(e.msg)"))
                else
                    throw(e)
                end
            end
        end
        push!(table, row)  
    end
    table
end

function parse_catmap_input(input_file_path)
    @pyinclude(input_file_path)
    
    reactions = ParsedReaction[]
    for r in py"rxn_expressions"
        push!(reactions, parse_reaction(r))
    end

    species_definitions = py"species_definitions"

    energies_file_path  = py"input_file"
    energy_table        = parse_energy_table(energies_file_path)
    
    species = specieslist(reactions, species_definitions, energy_table)


    bulk_pH = py"bulk_ph"
    
    Upzc                        = py"Upzc"
    potential_reference_scale   = py"potential_reference_scale"

    
    gas_thermo_mode             = py"gas_thermo_mode"
    adsorbate_thermo_mode       = py"adsorbate_thermo_mode"
    electrochemical_thermo_mode = py"electrochemical_thermo_mode"



end


# for each species included in the reactions
#   1. check that there exists an entry in species_defs
#   2. determine the type of the species
#   3a. if type == gas, complete species_def by using energy table
#   3b. if type == adsorbate, then determine the site, for each surface complete species_def by using energy table
#   3c. if type == tstate, then do same as for adsorbate
#   3d. if type == fictious, then ...
function specieslist(reactions::Vector{ParsedReaction}, species_defs, energy_table)
    rfictious = r"ele_g|OH_g"
    rgas = r"(?<species>[A-Za-z0-9]*)_g"
    radsorbate = r"(?<species>[A-Za-z0-9]*)\*_[^g]"
    rtstate = r"(?<species>[A-Za-z0-9\-]*)\*_[^g]"

    # collect all species in a set
    species = Set{String}()
    for (; educts, products) in reactions
        union!(species, first.(educts))
        union!(species, first.(products))
    end
    
    specieslist = Dict{String, AbstractSpecies}()
    for s in species
        if isfictious(s)
            specieslist[s] = FictiousSpecies()
        elseif isgas(s)
            specieslist[s] = GasSpecies()
        elseif isadsorbate(s)
            specieslist[s] = AdsorbateSpecies()
        elseif iststate(s)
            specieslist[s] = TstateSpecies()
        end
    end
        
    end
end