using Test

import IBZ
import IBZ.Plotting: plot_convexhulls


@testset "plotting" begin
    @testset "plot_convexhulls" begin

        cub_latvecs = IBZ.Lattices.genlat_CUB(1)
        real_latvecs = cub_latvecs
        atom_types = [0]
        coords = "Cartesian"
        atom_pos = Array([0 0 0]')
        convention = "ordinary"
        (fig,ax) = plot_convexhulls(real_latvecs,atom_types,atom_pos,coords,
            convention,true)
        (fig,ax) = plot_convexhulls(real_latvecs,atom_types,atom_pos,coords,
            convention,false)
        @test true

        sqr_latvecs = IBZ.Lattices.genlat_SQR(1)
        real_latvecs = sqr_latvecs
        atom_types = [0]
        coords = "Cartesian"
        atom_pos = Array([0 0]')
        convention = "ordinary"
        (fig,ax) = plot_convexhulls(real_latvecs,atom_types,atom_pos,coords,
            convention,true)
        (fig,ax) = plot_convexhulls(real_latvecs,atom_types,atom_pos,coords,
            convention,false)
        @test true
    end
end
