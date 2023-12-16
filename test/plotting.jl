using Test

using PyPlot
import SymmetryReduceBZ
import SymmetryReduceBZ.Plotting: plot_convexhulls


@testset "plotting" begin
    @testset "plot_convexhulls" begin

        cub_latvecs = SymmetryReduceBZ.Lattices.genlat_CUB(1)
        real_latvecs = cub_latvecs
        atom_types = [0]
        coords = "Cartesian"
        atom_pos = Array([0 0 0]')
        convention = "ordinary"
        ax = plot_convexhulls(real_latvecs,atom_types,atom_pos,coords,
            true,convention)

        ax = plot_convexhulls(real_latvecs,atom_types,atom_pos,coords,
            false,convention)
        @test true

        sqr_latvecs = SymmetryReduceBZ.Lattices.genlat_SQR(1)
        real_latvecs = sqr_latvecs
        atom_types = [0]
        coords = "Cartesian"
        atom_pos = Array([0 0]')
        convention = "ordinary"
        ax = plot_convexhulls(real_latvecs,atom_types,atom_pos,coords,
            true,convention)
        ax = plot_convexhulls(real_latvecs,atom_types,atom_pos,coords,
            false,convention)
        @test true
    end
end
