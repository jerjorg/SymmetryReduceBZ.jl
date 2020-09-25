using Test

import IBZ.Lattices
const lt = Lattices

import IBZ.Symmetry: calc_spacegroup, calc_pointgroup, calc_bz, calc_ibz,
    mapto_unitcell
import IBZ.Utilities: remove_duplicates
import QHull: chull
import PyCall: pyimport
sympy=pyimport("sympy")
import SymPy: symbols
x,y,z=symbols("x,y,z")


# Lattice vectors
# 2D
a=ℯ
b=π
θ=π/3.2
sqr_latvecs=lt.genlat_SQR(a)
rec_latvecs=lt.genlat_REC(a,b)
reci_latvecs=lt.genlat_RECI(a,b)
hxg_latvecs=lt.genlat_HXG(a)
obl_latvecs=lt.genlat_OBL(a,b,θ)

# 3D
a=π
cub_latvecs=lt.genlat_CUB(a)
fcc_latvecs=lt.genlat_FCC(a)
bcc_latvecs=lt.genlat_BCC(a)

c=ℯ
#1 c < a
#2 c > a
tet_latvecs=lt.genlat_TET(a,c)
bct_latvecs1=lt.genlat_BCT(a,c)
bct_latvecs2=lt.genlat_BCT(c,a)

b=5/7+(π+ℯ)/3
orc_latvecs=lt.genlat_ORC(a,b,c)

# Three face-centered orthorhombic cases
# 1/a² <,=,> 1/b²+1/c²
a=π
b=5/7+(π+ℯ)/3
c=ℯ
orcf_latvecs1=lt.genlat_ORCF(a,b,c)
a=π/2
b=5/7+(π+ℯ)/3
c=ℯ
orcf_latvecs2=lt.genlat_ORCF(a,b,c)
b=5/7
c=ℯ
a=b*c/√(b^2 + c^2)
orcf_latvecs3=lt.genlat_ORCF(a,b,c)

# One body-centered orthorhombic case
a=ℯ
b=π
c=4.0+1/3
orci_latvecs=lt.genlat_ORCI(a,b,c)
# One base-centered orthorhombic case
orcc_latvecs=lt.genlat_ORCC(a,b,c)
# One hexagonal case
hex_latvecs=lt.genlat_HEX(a,c)

# Two rhombohedral case α < π/2 and α > π/2
α=π/7+π/2
rhl_latvecs1=lt.genlat_RHL(a,α)
α=π/7
rhl_latvecs2=lt.genlat_RHL(a,α)

# One monoclinic case
mcl_latvecs=lt.genlat_MCL(a,b,c,α)

# Five base-centered monoclinic cases
# kᵧ > π/2
a=1+1/3
b=a+1/3
c=b+1/3
α=7*π/20
mclc_latvecs1=lt.genlat_MCLC(a,b,c,α)
# rlatvecs=get_recip_latvecs(mclc_latvecs1)
# ((a,b,c),(α,β,γ))=lt.get_latparams(rlatvecs)
# @show γ>π/2

# kᵧ == π/2
a=ℯ
b=π
c=4+1/3
α=1.0456607472366295
mclc_latvecs2=lt.genlat_MCLC(a,b,c,α)
# rlatvecs=get_recip_latvecs(mclc_latvecs2)
# ((a,b,c),(α,β,γ))=lt.get_latparams(rlatvecs)
# @show γ==π/2

# kᵧ < π/2, b*cos(α)/c + b^2*sin(α)^2/a^2 <1
a=ℯ
b=π
c=4+1/3
α=.3
mclc_latvecs3=lt.genlat_MCLC(a,b,c,α)
#@show b*cos(α)/c + b^2*sin(α)^2/a^2 <1
# rlatvecs=get_recip_latvecs(mclc_latvecs3)
# ((a,b,c),(α,β,γ))=lt.get_latparams(rlatvecs)
# @show γ<π/2

# kᵧ < π/2, b*cos(α)/c + b^2*sin(α)^2/a^2 > 1
a=1.2
b=1.3
c=1.4
α=0.6
mclc_latvecs4=lt.genlat_MCLC(a,b,c,α)
# @show b*cos(α)/c + b^2*sin(α)^2/a^2 > 1
# rlatvecs=get_recip_latvecs(mclc_latvecs4)
# ((a,b,c),(α,β,γ))=lt.get_latparams(rlatvecs)
# @show γ<π/2

# kᵧ < π/2, b*cos(α)/c + b^2*sin(α)^2/a^2 > 1
b=1.1
c=1.2
α=1.1
a=(b*√(c)*sin(α))/√(c - b*cos(α));
mclc_latvecs5=lt.genlat_MCLC(a,b,c,α)
# @show b*cos(α)/c + b^2*sin(α)^2/a^2 == 1
# rlatvecs=get_recip_latvecs(mclc_latvecs5)
# ((a,b,c),(α,β,γ))=lt.get_latparams(rlatvecs)
# @show γ==π/2

# Four triclinic cases
# kₐ>π/2,kᵦ>π/2,kᵧ>π/2 and kₐ = min(kₐ,kᵦ,kᵧ)
a=1.1
b=1.2
c=1.3
α=3*π/20
β=3.5*π/20
γ=4*π/20
tri_latvecs1=lt.genlat_TRI(a,b,c,α,β,γ)
rlatvecs=get_recip_latvecs(tri_latvecs1)
((a,b,c),(α,β,γ))=lt.get_latparams(rlatvecs)
# @show (α,β,γ)
# @show α>π/2
# @show β>π/2
# @show γ>π/2
# @show γ==minimum([α,β,γ])

# kₐ>π/2,kᵦ>π/2,kᵧ>π/2 and kₐ = max(kₐ,kᵦ,kᵧ)
a=1.1
b=1.2
c=1.3
α=3.5*π/20 + π/3
β=5*π/20 + π/3
γ=3.3*π/20 + π/3
tri_latvecs2=lt.genlat_TRI(a,b,c,α,β,γ)
rlatvecs=get_recip_latvecs(tri_latvecs2)
((a,b,c),(α,β,γ))=lt.get_latparams(rlatvecs)
# @show (α,β,γ)
# @show α<π/2
# @show β<π/2
# @show γ<π/2
# @show γ==maximum([α,β,γ])

# kₐ>π/2,kᵦ>π/2,kᵧ==π/2
a=1.1
b=1.2
c=1.3
α=3.5*π/20
β=5*π/20
γ=6*π/20-0.0188221060748125
tri_latvecs3=lt.genlat_TRI(a,b,c,α,β,γ)
rlatvecs=get_recip_latvecs(tri_latvecs3)
((a,b,c),(α,β,γ))=lt.get_latparams(rlatvecs)
# @show γ-π/2
# @show (α,β,γ)
# @show α>π/2
# @show β>π/
# @show γ≈π/2

# kₐ<π/2,kᵦ<π/2,kᵧ==π/2
a=1.1
b=1.2
c=1.3
α=2.5*π/20
β=2*π/20
γ=3.1*π/20 + 0.01079768761229742
tri_latvecs4=lt.genlat_TRI(a,b,c,α,β,γ)
rlatvecs=get_recip_latvecs(tri_latvecs4)
((a,b,c),(α,β,γ))=lt.get_latparams(rlatvecs)
# @show γ-π/2
# @show (α,β,γ)
# @show α>π/2
# @show β>π/2
# @show γ≈π/2

listreal_latvecs = [sqr_latvecs, rec_latvecs, reci_latvecs, hxg_latvecs,
    obl_latvecs, cub_latvecs, fcc_latvecs, bcc_latvecs, tet_latvecs,
    bct_latvecs1, bct_latvecs2, orc_latvecs, orcf_latvecs1, orcf_latvecs2,
    orcf_latvecs3, orci_latvecs, orcc_latvecs,hex_latvecs, rhl_latvecs1,
    rhl_latvecs2, mcl_latvecs, mclc_latvecs1, mclc_latvecs2, mclc_latvecs3,
    mclc_latvecs4, mclc_latvecs5, tri_latvecs1, tri_latvecs2, tri_latvecs3,
    tri_latvecs4]

pgsizes = [8, 4, 4, 12, 2, 48, 48, 48, 16, 16, 16, 8, 8, 8, 8, 8, 8, 24, 12, 12,
    4, 4, 4, 4, 4, 4, 2, 2, 2, 2]

atom_types=[0]
coords="Cartesian"
convention="ordinary"
bzformat="convex hull"
ibzformat="convex hull"

@testset "symmetry" begin
    @testset "calc_pointgroup" begin
        for (i,real_latvecs)=enumerate(listreal_latvecs)
            pointgroup=calc_pointgroup(real_latvecs)
            @test length(pointgroup) == pgsizes[i]
        end

        @test_throws ArgumentError calc_pointgroup([1 0])

    end

    @testset "calc_bz" begin
        convention="ordinary"
        bzformat="convex hull"
        for real_latvecs=listreal_latvecs
            bz=calc_bz(real_latvecs,convention,bzformat)
            pointgroup=calc_pointgroup(real_latvecs)
            # Rotate BZ vertices
            bzverts=reduce(hcat,[op*(bz.points[i,:]) for op=pointgroup
                for i=1:size(bz.points,1)])
            # BZ vertices should map to other BZ vertices under point group
            # operations
            @test size(remove_duplicates(bzverts),2) == size(bz.points,1)
        end

        bzformat="bad format"
        @test_throws ArgumentError calc_bz(listreal_latvecs[1],
            convention, bzformat)
    end

    @testset "calc_ibz" begin
        ibzformat="convex hull"
        for real_latvecs=listreal_latvecs
            if size(real_latvecs) == (2,2)
                atom_pos=Array([0 0]')
            else
                atom_pos=Array([0 0 0]')
            end
            pointgroup=calc_pointgroup(real_latvecs)
            bz=calc_bz(real_latvecs,convention,bzformat)
            ibz=calc_ibz(real_latvecs,atom_types,atom_pos,coords,ibzformat,
                convention);

            # Unfold IBZ
            unfoldpts=reduce(hcat,[op*(ibz.points[i,:]) for op=pointgroup
                        for i=1:size(ibz.points,1)])
            unfold_chull = chull(Array(unfoldpts'))
            unfoldpts=unfold_chull.points[unfold_chull.vertices,:]
            @test size(unfoldpts,1) == size(bz.points[bz.vertices,:],1)
            @test all([any([isapprox(unfoldpts[i,:],bz.points[j,:])
                for i=1:size(unfoldpts,1)]) for j=1:size(bz.points,1)])
        end

        ibzformat="bad format"
        @test_throws ArgumentError calc_ibz(listreal_latvecs[1], [0],
            Array([0 0]'), "lattice", ibzformat, "ordinary")

        passed = false
        ibz = calc_ibz(listreal_latvecs[1], [0,0], Array([0 0; 0.5 0.5]'),
            "Cartesian", "half-space", "ordinary")
        passed = true
        @test passed
    end

    @testset "mapto_unitcell" begin
        pt = [2,3]
        basis = [1 0; 0 1]
        inv_basis = [1 0; 0 1]
        coords = "Cartesian"
        @test mapto_unitcell(pt,basis,inv_basis,coords) ≈ [0,0]

        pt = [2,3]
        basis = [1/2 1/2; 1/2 -1/2]
        inv_basis = inv(basis)
        coords = "Cartesian"
        @test mapto_unitcell(pt,basis,inv_basis,coords) ≈ [0,0]

        pt = [1,1]
        basis = [1/2 1/2; 1/2 -1/2]
        inv_basis = inv(basis)
        coords = "lattice"
        @test mapto_unitcell(pt,basis,inv_basis,coords) ≈ [0,0]

        coords = "spherical"
        @test_throws ArgumentError mapto_unitcell(pt,basis,inv_basis,coords)
    end

    @testset "symmetry" begin
        @testset "calc_spacegroup" begin

            function compareSG(real_latvecs, atom_types, atom_pos, coords,
                findsymSG,atol=1e-8,rtol=1e-8)

                # Make a matrix and vector of each space group operator that can
                # easily be compared.
                findsymP=[]
                findsymT=[]
                for opᵢ=1:length(findsymSG)
                    tmp=reduce(hcat,[Array{Float64,1}(
                        [sympy.poly(findsymSG[opᵢ][i],x,y,z).coeff_monomial(var)
                        for var=[x,y,z,1]]) for i=1:3])
                    append!(findsymP,[tmp[1:3,1:3]])
                    append!(findsymT,[tmp[end,:]])
                end

                # We first calculate the space group in Cartesian coordinates
                # (the point operators operate on points in Cartesian
                # coordinates and the translations are in Cartesian
                # coordinates).
                (ibzT,ibzP)=calc_spacegroup(real_latvecs,atom_types,atom_pos,
                    coords)
                # We put the space group in lattice coordinates (the point
                # operators operate on points in lattice coordinates and the
                # translations are in lattice coordinates).
                ibzSG=[round.(inv(real_latvecs)*ibzP[i]*real_latvecs,
                    digits=8)*[x,y,z] + round.(inv(real_latvecs)*ibzT[i],
                    digits=8) for i=1:length(ibzP)]
                # Make a matrix and vector of each space group operator that can
                # easily be compared.
                ibzP=[]
                ibzT=[]
                for opᵢ=1:length(ibzSG)
                    tmp=reduce(hcat,[Array{Float64,1}([sympy.poly(ibzSG[opᵢ][i],
                        x,y,z).coeff_monomial(var) for var=[x,y,z,1]])
                        for i=1:3])
                    append!(ibzP,[tmp[1:3,1:3]])
                    append!(ibzT,[tmp[end,:]])
                end

                same = []
                for i=1:length(ibzP), j=1:length(findsymP)
                    if isapprox(ibzP[i],findsymP[j],atol=atol,rtol=rtol) &&
            		   isapprox(ibzT[i],findsymT[j],atol=atol,rtol=rtol)
                        append!(same,i)
                    end
                end
                sort(same) == range(1,stop=length(ibzT)) &&
                    length(ibzP) == length(findsymP)
            end

            # Number of atom types and positions are different.
            a=1
            real_latvecs=lt.genlat_CUB(a)
            atom_pos = Array([0 0 0; 0.25 0.25 0]')
            # findsym FCC convention
            creal_latvecs=reduce(hcat,[real_latvecs*[-1,0,0],
                real_latvecs*[0,1,-1], real_latvecs*[0,-1,-1]])
            catom_pos = Array([0 0 0; 0 0 0.75]')
            atom_types=[0]
            coords = "lattice"

            @test_throws ArgumentError calc_spacegroup(creal_latvecs,atom_types,
                atom_pos,coords)

            # SC
            a=1
            real_latvecs=lt.genlat_CUB(a)
            atom_pos = Array([0 0 0; 0.25 0.25 0]')
            # findsym FCC convention
            creal_latvecs=reduce(hcat,[real_latvecs*[-1,0,0],
                real_latvecs*[0,1,-1], real_latvecs*[0,-1,-1]])
            catom_pos = Array([0 0 0; 0 0 0.75]')
            atom_types=[0,1]
            coords = "lattice"

            findsymSG =
            [[x,y,z],
            [-x,-y,z],
            [-x,y,z],
            [x,-y,z]]

            sctest = compareSG(creal_latvecs,atom_types,catom_pos,coords,
                findsymSG)
            @test sctest

            # BCC
            a=1
            real_latvecs=lt.genlat_BCC(a)
            # findsym BCC convention
            creal_latvecs=reduce(hcat,[real_latvecs*[-1,0,1],
                real_latvecs*[1,-1,0],real_latvecs*[1,1,1]])
            catom_pos = Array([0 0 0; 0 0 0.5]')
            atom_types=[0,1]
            coords = "lattice";

            findsymSG =
            [[x,y,z],
            [x-y,x,z],
            [-y,x-y,z],
            [-x,-y,z],
            [-x+y,-x,z],
            [y,-x+y,z],
            [x-y,-y,-z],
            [x,x-y,-z],
            [y,x,-z],
            [ -x+y,y,-z],
            [ -x,-x+y,-z],
            [ -y,-x,-z],
            [ -x,-y,-z],
            [ -x+y,-x,-z],
            [ y,-x+y,-z],
            [ x,y,-z],
            [ x-y,x,-z],
            [ -y,x-y,-z],
            [ -x+y,y,z],
            [ -x,-x+y,z],
            [ -y,-x,z],
            [ x-y,-y,z],
            [ x,x-y,z],
            [ y,x,z]]
            bcctest = compareSG(creal_latvecs, atom_types, catom_pos, coords,
                findsymSG)
            @test bcctest

            # FCC
            a=1
            real_latvecs=lt.genlat_FCC(a)
            # findsym FCC convention
            creal_latvecs=reduce(hcat,[real_latvecs*[0,0,1],
                real_latvecs*[1,1,-1], real_latvecs*[-1,1,0]])
            catom_pos = Array([0 0 0; 0.5 0 0]')
            atom_types=[0,1]
            coords = "lattice";

            findsymSG =
            [[x,y,z],
            [x,-y,-z],
            [-x,y,-z],
            [-x,-y,z],
            [-x,-y,-z],
            [-x,y,z],
            [x,-y,z],
            [x,y,-z]]

            fcctest = compareSG(creal_latvecs, atom_types, catom_pos, coords, findsymSG)
            @test fcctest

            # BCT
            a=2.2
            c=1
            real_latvecs=lt.genlat_BCT(a,c)
            atom_pos=Array([0 0 0; 0.5 0 0]')
            # findsym FCC convention
            creal_latvecs=reduce(hcat,[real_latvecs*[1,1,0],
                real_latvecs*[-1,1,0], real_latvecs*[0,0,1]])
            catom_pos = Array([0 0 0; 0 0.5 0.5]')
            atom_types=[0,1]
            coords = "lattice";

            findsymSG =
            [[x,y,z],
            [-x,y,-z],
            [-x,-y,-z],
            [x,-y,z]]

            bcttest = compareSG(creal_latvecs, atom_types, catom_pos, coords,
                findsymSG)
            @test bcttest

            # HEX
            a=1.3
            c=1.5
            real_latvecs=lt.genlat_HEX(a,c)'
            atom_pos=Array([0 0 0; 0 0.5 0.5]')
            # findsym conventional lattice vectors
            creal_latvecs=reduce(hcat,[real_latvecs*[1,0,0],
                real_latvecs*[0,1,0],real_latvecs*[0,0,1]])
            catom_pos = Array([0 0 0; 0.0 0.5 0.5]')
            atom_types=[0,1]
            coords = "lattice";


            findsymSG =
            [[x,y,z],
            [x,-y,-z],
            [-x,y,-z],
            [-x,-y,z],
            [-x,-y,-z],
            [-x,y,z],
            [x,-y,z],
            [x,y,-z]]

            hextest = compareSG(creal_latvecs, atom_types, catom_pos, coords,
                findsymSG)
            @test hextest

            # MCL
            a=1.3
            b=2.4
            c=3.5
            α=0.3
            real_latvecs=lt.genlat_MCL(a,b,c,α)
            atom_pos=Array([0 0 0; 0.5 0 0.5]')
            # findsym conventional lattice vectors
            creal_latvecs=reduce(hcat,[real_latvecs*[0,-1,1],
                real_latvecs*[-1,0,0],real_latvecs*[0,-2,1]])
            catom_pos = Array([0 0 0; 0 0.5 0.5]')
            atom_types=[0,1]
            coords = "lattice";

            findsymSG =
            [[x,y,z],
            [-x,y,-z],
            [-x,-y,-z],
            [x,-y,z]]

            mcltest = compareSG(creal_latvecs,atom_types,catom_pos,coords,
                findsymSG)
            @test mcltest

            # MCLC
            a=1.2
            b=1.4
            c=1.6
            α=0.7
            real_latvecs=lt.genlat_MCLC(a,b,c,α)
            atom_pos=Array([0 0 0; 0 0.5 0]')
            # findsym conventional lattice vectors
            creal_latvecs=reduce(hcat,[real_latvecs*[-1,0,0],
                real_latvecs*[0,-1,0],real_latvecs*[-1,-1,1]])
            catom_pos = Array([0 0 0; 0 0.5 0]')
            atom_types=[0,1]
            coords = "lattice";

            findsymSG =
            [[x,y,z],
            [-x,-y,-z]]

            mclctest = compareSG(creal_latvecs, atom_types, catom_pos, coords,
                findsymSG)
            @test mclctest

            # ORC
            a=1.2
            b=1.3
            c=1.5
            real_latvecs=lt.genlat_ORC(a,b,c)
            atom_pos=Array([0 0 0; 0.5 0.5 0.5]')
            creal_latvecs=reduce(hcat,[real_latvecs*[1,0,0],
                real_latvecs*[0,1,0],real_latvecs*[0,0,1]])
            catom_pos = Array([0 0 0; 0.5 0.5 0.5]')

            findsymSG =
            [[x,y,z],
            [x,-y,-z],
            [-x,y,-z],
            [-x,-y,z],
            [-x,-y,-z],
            [-x,y,z],
            [x,-y,z],
            [x,y,-z]]

            orctest = compareSG(creal_latvecs, atom_types, catom_pos, coords,
                findsymSG)
            @test orctest

            # ORCC
            a=1.2
            b=1.3
            c=1.5
            real_latvecs=lt.genlat_ORCC(a,b,c)
            atom_pos=Array([0 0 0; 0.5 0.5 0.5]')
            creal_latvecs=reduce(hcat,[real_latvecs*[-1,1,0],
                real_latvecs*[-1,-1,0],real_latvecs*[0,0,1]])
            catom_pos = Array([0 0 0; 0.5 0 0.5]')

            findsymSG =
            [[x,y,z],
            [x,-y,-z],
            [-x,y,-z],
            [-x,-y,z],
            [-x,-y,-z],
            [-x,y,z],
            [x,-y,z],
            [x,y,-z]]

            orcctest = compareSG(creal_latvecs, atom_types, catom_pos, coords,
                findsymSG)
            @test orcctest

            # ORCF
            a=1.2
            b=1.3
            c=1.5
            real_latvecs=lt.genlat_ORCF(a,b,c)
            atom_pos=Array([0 0 0; 0.5 0.5 0.5]')
            creal_latvecs=reduce(hcat,[real_latvecs*[1,-1,1],
                real_latvecs*[-1,-1,1],real_latvecs*[1,-1,-1]])
            catom_pos = Array([0 0 0; 0.5 0.5 0.5]')

            findsymSG =
            [[x,y,z],
            [x,-y,-z],
            [-x,y,-z],
            [-x,-y,z],
            [-x,-y,-z],
            [-x,y,z],
            [x,-y,z],
            [x,y,-z]]

            orcftest = compareSG(creal_latvecs, atom_types, catom_pos, coords,
                findsymSG)
            @test orcftest

            # ORCI
            a=1.2
            b=1.3
            c=1.5
            real_latvecs=lt.genlat_ORCI(a,b,c)
            atom_pos=Array([0 0 0; 0.5 0.5 0.5]')
            creal_latvecs=reduce(hcat,[real_latvecs*[1,0,0],
                real_latvecs*[0,1,0],real_latvecs*[0,0,1]])
            catom_pos = Array([0 0 0; 0.5 0.5 0.5]')

            findsymSG =
            [[x,y,z],
            [-x,-y,-z]]

            orcitest = compareSG(creal_latvecs, atom_types, catom_pos, coords,
                findsymSG)
            @test orcitest

            # RHL
            a=1.2
            α=0.8
            real_latvecs=lt.genlat_RHL(a,α)
            atom_pos=Array([0 0 0; 0.5 0.5 0.5]')
            creal_latvecs=reduce(hcat,[real_latvecs*[0,1,-1],
                real_latvecs*[-1,0,1],real_latvecs*[1,1,1]])
            catom_pos = Array([0 0 0; 0 0 0.5]')

            findsymSG =
            [[x,y,z],
            [x-y,x,z],
            [-y,x-y,z],
            [-x,-y,z],
            [-x+y,-x,z],
            [y,-x+y,z],
            [x-y,-y,-z],
            [x,x-y,-z],
            [y,x,-z],
            [ -x+y,y,-z],
            [ -x,-x+y,-z],
            [ -y,-x,-z],
            [ -x,-y,-z],
            [ -x+y,-x,-z],
            [ y,-x+y,-z],
            [ x,y,-z],
            [ x-y,x,-z],
            [ -y,x-y,-z],
            [ -x+y,y,z],
            [ -x,-x+y,z],
            [ -y,-x,z],
            [ x-y,-y,z],
            [ x,x-y,z],
            [ y,x,z]]

            rhltest = compareSG(creal_latvecs, atom_types, catom_pos, coords,
                findsymSG)
            @test rhltest

            # TET
            a=1.2
            c=1.5
            real_latvecs=lt.genlat_TET(a,c)
            atom_pos=Array([0 0 0; 0.5 0.5 0.5]')
            creal_latvecs=reduce(hcat,[real_latvecs*[1,0,0],
                real_latvecs*[0,1,0],real_latvecs*[0,0,1]])
            catom_pos = Array([0 0 0; 0.5 0.5 0.5]')

            findsymSG =
            [[x,y,z],
            [x,-y,-z],
            [-x,y,-z],
            [-x,-y,z],
            [-y,-x,-z],
            [-y,x,z],
            [y,-x,z],
            [y,x,-z],
            [-x,-y,-z],
            [ -x,y,z],
            [ x,-y,z],
            [ x,y,-z],
            [ y,x,z],
            [ y,-x,-z],
            [ -y,x,-z],
            [ -y,-x,z]]

            tettest = compareSG(creal_latvecs, atom_types, catom_pos, coords,
                findsymSG)
            @test tettest

            # TRI
            a=1.2
            b=1.4
            c=1.5
            α=0.5
            β=0.7
            γ=0.85
            real_latvecs=lt.genlat_TRI(a,b,c,α,β,γ)
            atom_pos=Array([0 0 0; 0.5 0.5 0.5]')
            creal_latvecs=reduce(hcat,[real_latvecs*[0,-1,1],
                real_latvecs*[-1,0,1],real_latvecs*[-1,-1,1]])
            catom_pos = Array([0 0 0; 0 0 0.5]')

            findsymSG =
            [[x,y,z],
            [-x,-y,-z]]

            tritest = compareSG(creal_latvecs, atom_types, catom_pos, coords,
                findsymSG)
            @test tritest
        end
    end
end
