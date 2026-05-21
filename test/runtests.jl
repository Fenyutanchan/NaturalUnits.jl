using NaturalUnits
using Test

const test_files = [
    "test_constructors",
    "test_conversion",
    "test_arithmetic",
    "test_math_functions",
    "test_comparison",
    "test_unit_conversion",
    "test_constants",
    "test_edge_cases",
    "test_broadcast",
]

@testset "NaturalUnits.jl" begin
    for f in test_files
        @testset "$f" begin
            include("$f.jl")
        end
    end
end
