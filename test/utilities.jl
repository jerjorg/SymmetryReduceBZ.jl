using Test

import SymmetryReduceBZ.Utilities: affine_trans, contains, edgelengths,
	get_uniquefacets, mapto_xyplane, remove_duplicates, sample_circle,
	sample_sphere, shoelace, sortpts_perm, unique_points

import SymmetryReduceBZ.Lattices: genlat_SQR, genlat_REC, genlat_RECI, genlat_HXG,
    genlat_OBL, genlat_CUB

import SymmetryReduceBZ.Symmetry: calc_bz

# Lattice vectors
# 2D
a=5/7.
b=11/6
θ=π/3.2
sqr_latvecs=genlat_SQR(a)
rec_latvecs=genlat_REC(a,b)
reci_latvecs=genlat_RECI(a,b)
hxg_latvecs=genlat_HXG(a)
obl_latvecs=genlat_OBL(a,b,θ)

listreal_latvecs = [sqr_latvecs, rec_latvecs, reci_latvecs, hxg_latvecs,
    obl_latvecs]
conventions=["ordinary", "angular"]
bzformat = "convex hull"

@testset "Utilities" begin

    @testset "affine_trans" begin

        pts = Array([0.5 0.5 0.5; 0.5 -0.5 0.5; -0.5 0.5 0.5; -0.5 -0.5 0.5]')
        A = affine_trans(pts)
        transpts=A*[pts; [1 1 1 1]]

        @test transpts[3,:] ≈ zeros(4)
        @test transpts[4,:] ≈ ones(4)

        pts = Array([0.5 0.5 -0.5; 0.5 -0.5 -0.5; -0.5 0.5 -0.5;
            -0.5 -0.5 -0.5]')
        A = affine_trans(pts)
        transpts=A*[pts; [1 1 1 1]]

        @test transpts[3,:] ≈ zeros(4)
        @test transpts[4,:] ≈ ones(4)
    end

    @testset "contains" begin
        pt = [1,2,3]
        pts = Array([1 2 3; 1.1 2.2 3.3; 2 3 4]')
        @test contains(pt,pts) == true

        pt = [1+1e-10,2,3]
        pts = Array([1 2 3; 1.1 2.2 3.3; 2 3 4]')
        @test contains(pt,pts) == true

        pt = [1+1e-6,2,3]
        pts = Array([1 2 3; 1.1 2.2 3.3; 2 3 4]')
        @test contains(pt,pts) == false

        rtol = 1e-6
        pt = [1+1e-6,2,3]
        pts = Array([1 2 3; 1.1 2.2 3.3; 2 3 4]')
        @test contains(pt,pts,rtol=rtol) == true
    end

    @testset "edgelengths" begin
        radii = [0.2, 1, π, 1e10, 1e-10]
        basis = [1 0; 0 1]
        for radius=radii
            lengths = edgelengths(basis,radius)
            @test lengths ≈ 2*radius.*[1,1]
        end

        radius = -π
        @test_throws ArgumentError edgelengths(basis,radius)

        radius = 1
        basis = [1 2 3; 0 0 1]
        @test_throws ArgumentError edgelengths(basis,radius)
    end

    @testset "get_uniquefacets" begin
        a=1
        real_latvecs=genlat_CUB(a)
		atom_types = [0]
		atom_pos = Array([0 0 0]')
        convention = "ordinary"
        bzformat = "convex hull"
		coords = "Cartesian"
		makeprim = false

		bz = calc_bz(real_latvecs,atom_types,atom_pos,coords,bzformat,makeprim,
			convention)

        facetsᵢ = get_uniquefacets(bz)
        facets = [bz.points[facetsᵢ[i],:] for i=1:length(facetsᵢ)]
		truefacets = [[0.5 -0.5 0.5; 0.5 -0.5 -0.5; 0.5 0.5 -0.5; 0.5 0.5 0.5],
			[-0.5 -0.5 -0.5; 0.5 -0.5 -0.5; 0.5 0.5 -0.5; -0.5 0.5 -0.5],
			[-0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 -0.5; -0.5 0.5 -0.5],
			[-0.5 -0.5 -0.5; 0.5 -0.5 -0.5; 0.5 -0.5 0.5; -0.5 -0.5 0.5],
			[-0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 -0.5 0.5; -0.5 -0.5 0.5],
			[-0.5 -0.5 0.5; -0.5 -0.5 -0.5; -0.5 0.5 -0.5; -0.5 0.5 0.5]]

        @test truefacets ≈ facets
    end

    @testset "mapto_xyplane" begin
        facets = [[0.5 -0.5 -0.5; 0.5 0.5 -0.5; 0.5 0.5 0.5; 0.5 -0.5 0.5],
            [0.5 0.5 -0.5; -0.5 0.5 -0.5; -0.5 -0.5 -0.5; 0.5 -0.5 -0.5],
            [0.5 0.5 -0.5; -0.5 0.5 -0.5; -0.5 0.5 0.5; 0.5 0.5 0.5],
            [0.5 -0.5 -0.5; 0.5 -0.5 0.5; -0.5 -0.5 0.5; -0.5 -0.5 -0.5],
            [0.5 0.5 0.5; 0.5 -0.5 0.5; -0.5 -0.5 0.5; -0.5 0.5 0.5],
            [-0.5 0.5 -0.5; -0.5 0.5 0.5; -0.5 -0.5 0.5; -0.5 -0.5 -0.5]]

        @test all([
            mapto_xyplane(facets[6]') ≈ [0.0 1.0 1.0 0.0; 0.0 0.0 1.0 1.0] for
            i=1:length(facets)])
        end

    @testset "unique_points" begin
        pts = Array([1.1 1.2 1.3; 1.1+1e-10 1.2+1e-11 1.3+1e-9; 0.1 0.2 0.3]')
        upts = unique_points(pts)
        @test size(upts) == (3,2)

        pts = Array([1.1 1.2 1.3; 1.1+1e-7 1.2+1e-11 1.3+1e-10; 0.1 0.2 0.3]')
        upts = unique_points(pts)
        @test size(upts) == (3,3)

        rtol=1e-7
        pts = Array([1.1 1.2 1.3; 1.1+1e-7 1.2+1e-11 1.3+1e-10; 0.1 0.2 0.3]')
        upts = unique_points(pts,rtol=rtol)
        @test size(upts) == (3,2)
    end

    @testset "remove_duplicates" begin
        pts = [1,1,2,2,3,3,4,4]
        @test isapprox(remove_duplicates(pts),[1,2,3,4])
    end

    @testset "sample_circle" begin
        basis = [1 0; 0 1]
        radius = 1
        offset = [0,0]
        pts=sample_circle(basis,radius,offset)
        cpts = Array([0 0; 1 0; -1 0; 0 1; 0 -1]')
        @test all([contains(cpts[:,i],pts) for i=1:size(cpts,2)])
        @test size(pts,2) == size(cpts,2)

        basis = [1 0; 0 1]
        radius = √(2)
        offset = [0,0]
        pts=sample_circle(basis,radius,offset)
        cpts = Array([0 0; 1 0; -1 0; 0 1; 0 -1; 1 1; 1 -1; -1 1; -1 -1]')
        @test all([contains(cpts[:,i],pts) for i=1:size(cpts,2)])
        @test size(pts,2) == size(cpts,2)
    end

    @testset "sample_sphere" begin
        basis = [1 0 0; 0 1 0; 0 0 1]
        radius = 1
        offset = [0,0,0]
        pts=sample_sphere(basis,radius,offset)
        cpts = Array([0 0 0; 1 0 0 ; -1 0 0 ; 0 1 0 ; 0 -1 0; 0 0 1; 0 0 -1]')
        @test all([contains(cpts[:,i],pts) for i=1:size(cpts,2)])
        @test size(pts,2) == size(cpts,2)

        basis = [1 0 0; 0 1 0; 0 0 1]
        radius = √(2)
        offset = [0,0,0]
        pts=sample_sphere(basis,radius,offset)
        cpts = Array([0 0 0; 1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1;
                1 1 0; 1 -1 0; -1 1 0; -1 -1 0; 1 0 1; -1 0 1; 1 0 -1; -1 0 -1;
                0 1 1; 0 1 -1; 0 -1 1; 0 -1 -1]')

        @test all([contains(cpts[:,i],pts) for i=1:size(cpts,2)])
        @test size(pts,2) == size(cpts,2)
    end

    @testset "shoelace" begin
        for (real_latvecs,convention)=Iterators.product(listreal_latvecs,
                conventions)
			atom_types = [0]
			atom_pos = Array([0 0 0]')
			coords = "Cartesian"
			makeprim = false
			bz = calc_bz(real_latvecs,atom_types,atom_pos,coords,bzformat,
				makeprim,convention)

            svol = shoelace(bz.points[bz.vertices,:]')
            @test svol ≈ bz.volume
        end
    end

    @testset "sortpts_perm" begin
        pts = [0.5 -0.5 0.5; 0.5 -0.5 -0.5; 0.5 0.5 -0.5; 0.5 0.5 0.5]'
        @test sortpts_perm(pts) == [1,2,3,4]

        pts = [-0.5 -0.5 0.5; -0.5 -0.5 -0.5; -0.5 0.5 -0.5; -0.5 0.5 0.5]'
        @test sortpts_perm(pts) == [1,2,3,4]
    end

end
