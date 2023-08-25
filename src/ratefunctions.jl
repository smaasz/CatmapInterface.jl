function isgas(sp::BasicSymbolic{Real})
    repr(operation(sp))[-2:end] == "_g"
end

function gases(sps::Vector{BasicSymbolic{Real}})
    @views sps[isgas.(sps)]
end

function isadsorbate(sp::BasicSymbolic{Real})
    repr(operation(sp))[-2:end] == "_t"
end

function adsorbates(sps::Vector{BasicSymbolic{Real}})
    @views sps[isadsorbate.(sps)]
end

function istransition_state(sp::BasicSymbolic{Real})
    repr(operation(sp))[-3:end] == "_TS"
end

function transition_states(sps::Vector{BasicSymbolic{Real}})
    @views sps[istransition_state.(sps)]
end

function ratefunctions(params::CatmapParams)
    sps = species(params.rn)
    gases = gases(sps)
    adsorbates = adsorbates(sps)
    transition_states = transition_states(sps)


    return rfs
end