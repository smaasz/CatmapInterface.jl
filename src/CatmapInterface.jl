__precompile__()
module CatmapInterface

using PyCall

function __init__()  
    pyimport_conda("ase.thermochemistry", "ase")
    pyimport_conda("ase.build", "ase")
    py"""
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from ase.build import molecule

def get_thermal_correction_adsorbate(T, frequencies):
    thermo = HarmonicThermo(frequencies)
    return thermo.get_helmholtz_energy(T, verbose=False)

def get_thermal_correction_ideal_gas(T, frequencies, symmetrynumber, geometry, spin, name):
    thermo = IdealGasThermo(frequencies, geometry, atoms=molecule(name), symmetrynumber=symmetrynumber, spin=spin) 
    H = thermo.get_enthalpy(T, verbose=False)
    S = thermo.get_entropy(T, 1.0e5, verbose=False)

    free_energy = H-T*S
    return free_energy
	"""
    @pyinclude(joinpath(@__DIR__, "../data/parameter_data.py"))
end

using Catalyst
using ModelingToolkit
using DelimitedFiles
using LessUnitful
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)
include("species.jl")
export AbstractSpecies, GasSpecies, AdsorbateSpecies, SiteSpecies, TStateSpecies, FictiousSpecies
include("interface.jl")
export parse_catmap_input
include("corrections.jl")
include("reaction_network.jl")
export create_reaction_network, generate_function

end