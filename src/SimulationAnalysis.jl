
module SimulationAnalysis
import Base.show, Quickhull
using HDF5, Polyester, Tullio, LoopVectorization, Base.Threads, Parsers, DelimitedFiles, Random, IfElse, Dierckx, ProgressMeter, OffsetArrays, ChunkSplitters, CellListMap, StaticArrays


abstract type Simulation end
import Base.step

include("Simulation.jl")
include("LoadData.jl")
include("Kspace.jl")
include("CorrelationFunction.jl")
include("DensityModes.jl")
include("IntermediateScatteringFunction.jl")
include("StructureFactors.jl")
include("Forces.jl")
include("RadialDistributionFunction.jl")
include("MeanSquaredDisplacement.jl")
include("F4_diagonal.jl")
include("OverlapFunction.jl")
include("Neighborlists.jl")
include("BondBreakingParameter.jl")



function show(io::IO,  ::MIME"text/plain", s::Union{KSpace, Simulation})
    println(io, "This is a $(typeof(s)).")
    println(io, "It contains ")
    for fieldname in fieldnames(typeof(s))
        if getfield(s, fieldname) isa Union{Int, Float64, String}
            println(io, "$(fieldname): $(getfield(s, fieldname))")
        else
            println(io, "$(fieldname): $(typeof(getfield(s, fieldname)))")
        end
    end
end

end # module

