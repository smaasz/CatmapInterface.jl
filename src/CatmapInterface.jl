module CatmapInterface

ENV["PYCALL_JL_RUNTIME_PYTHON"] = Sys.which("python3")
using PyCall

using DifferentialEquations
using ProgressBars
using NonlinearSolve

function __init__()
    py"""
    from catmap import ReactionModel

    def catmap_kinetic_model(setup_file):

        model = ReactionModel(setup_file=setup_file)

        # some solver parameters have to be set manually (?!)
        import mpmath as mp
        model.solver._mpfloat = mp.mpf
        model.solver._math = mp
        model.solver._matrix = mp.matrix

        model.solver.compile() # compiles all templates, here (rate_constants) are needed

        # define descriptor grid
        if isinstance(model.resolution, int):
            model.resolution = [model.resolution] * len(model.descriptor_names)
        
        def linspace(start, stop, nitems):
            
            if start >= stop:
                return [start]
        
            delta = (stop - start)/(nitems-1)
            r = []
            for i in range(nitems):
                r.append(start + i * delta)
            return r

        def get_idx_tuple(idx, lengths):
            idx_tuple = []
            rest = idx
            for length in lengths:
                idx_tuple.append(rest % length)
                rest = rest // length
            return idx_tuple

        ndescriptor_tuples = 1
        for res in model.resolution:
            ndescriptor_tuples *= res
        
        descriptor_ranges = [linspace(_range[0], _range[1], res) for _range, res in zip(model.descriptor_ranges, model.resolution)]
        
        descriptor_grid = [0] * ndescriptor_tuples
        rxn_params_grid = [0] * ndescriptor_tuples
        for i in range(ndescriptor_tuples):
            descriptor_values = [descriptor_range[idx] for idx, descriptor_range in zip(get_idx_tuple(i, model.resolution), descriptor_ranges)]
            descriptor_grid[i] = descriptor_values
            rxn_params_grid[i] = model.scaler.get_rxn_parameters(descriptor_values)


        # the rate constants should not depend on the coverage 
        # as long as the no adsorbate interactions are assumed (model.adsorbate_interaction_model == "ideal")
        coverages = [0.0] * len(model.adsorbate_names)

        rate_constants_grid = [0] * ndescriptor_tuples
        for i in range(ndescriptor_tuples):
            params = rxn_params_grid[i]
            kfs, krs, dkfs, dkrs = model.solver.rate_constants(
                params, coverages, model.solver._gas_energies, model.solver._site_energies, 
                model.solver.temperature, model.solver.interaction_function, 
                model.solver._mpfloat, model.solver._matrix, 
                model.solver._math.exp, model.solver._math.sqrt
            )
            rate_constants_grid[i] = [kfs, krs, dkfs, dkrs]

        ss_eqs = model.solver.rate_equations()

        # def indent_string(string,levels):
        #     lines = string.split('\n')
        #     indention = '\n'+'    '*levels
        #     return indention.join(lines)

        # steady_state_expressions = indent_string('\n    '.join(ss_eqs), 0)
        steady_state_expressions = '\n    '.join(ss_eqs)

        return descriptor_grid, rate_constants_grid, steady_state_expressions, model.elementary_rxns, model.gas_names
    """
end


include("kineticmodel.jl")
export get_catmap_output
export convert_to_julia

include("solver.jl")
export compute_coverage_map
export compute_coverage

end
