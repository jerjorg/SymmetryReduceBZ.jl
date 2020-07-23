using Test,IBZ

# Lattice vectors
# 2D
a=ℯ
b=π
θ=π/3.2
sqr_latvecs=genlat_SQR(a)
rec_latvecs=genlat_REC(a,b)
reci_latvecs=genlat_RECI(a,b)
hxg_latvecs=genlat_HXG(a)
obl_latvecs=genlat_OBL(a,b,θ)

# 3D
a=π
cub_latvecs=genlat_CUB(a)
fcc_latvecs=genlat_FCC(a)
bcc_latvecs=genlat_BCC(a)

c=ℯ
#1 c < a
#2 c > a
tet_latvecs=genlat_TET(a,c)
bct_latvecs1=genlat_BCT(a,c)
bct_latvecs2=genlat_BCT(c,a)

b=5/7+(π+ℯ)/3
orc_latvecs=genlat_ORC(a,b,c)

# Three face-centered orthorhombic cases
# 1/a² <,=,> 1/b²+1/c²
a=π
b=5/7+(π+ℯ)/3
c=ℯ
orcf_latvecs1=genlat_ORCF(a,b,c)
a=π/2
b=5/7+(π+ℯ)/3
c=ℯ
orcf_latvecs2=genlat_ORCF(a,b,c)
b=5/7
c=ℯ
a=b*c/√(b^2 + c^2)
orcf_latvecs3=genlat_ORCF(a,b,c)

# One body-centered orthorhombic case
a=ℯ
b=π
c=4.0+1/3
orci_latvecs=genlat_ORCI(a,b,c)
# One base-centered orthorhombic case
orcc_latvecs=genlat_ORCC(a,b,c)
# One hexagonal case
hex_latvecs=genlat_HEX(a,c)

# Two rhombohedral case α < π/2 and α > π/2
α=π/7+π/2
rhl_latvecs1=genlat_RHL(a,α)
α=π/7
rhl_latvecs2=genlat_RHL(a,α)

# One monoclinic case
mcl_latvecs=genlat_MCL(a,b,c,α)

# Five base-centered monoclinic cases
# kᵧ > π/2
a=1+1/3
b=a+1/3
c=b+1/3
α=7*π/20
mclc_latvecs1=genlat_MCLC(a,b,c,α)
# rlatvecs=get_recip_latvecs(mclc_latvecs1)
# ((a,b,c),(α,β,γ))=IBZ.Lattices.get_latparams(rlatvecs)
# @show γ>π/2

# kᵧ == π/2
a=ℯ
b=π
c=4+1/3
α=1.0456607472366295
mclc_latvecs2=genlat_MCLC(a,b,c,α)
# rlatvecs=get_recip_latvecs(mclc_latvecs2)
# ((a,b,c),(α,β,γ))=IBZ.Lattices.get_latparams(rlatvecs)
# @show γ==π/2

# kᵧ < π/2, b*cos(α)/c + b^2*sin(α)^2/a^2 <1
a=ℯ
b=π
c=4+1/3
α=.3
mclc_latvecs3=genlat_MCLC(a,b,c,α)
#@show b*cos(α)/c + b^2*sin(α)^2/a^2 <1
# rlatvecs=get_recip_latvecs(mclc_latvecs3)
# ((a,b,c),(α,β,γ))=IBZ.Lattices.get_latparams(rlatvecs)
# @show γ<π/2

# kᵧ < π/2, b*cos(α)/c + b^2*sin(α)^2/a^2 > 1
a=1.2
b=1.3
c=1.4
α=0.6
mclc_latvecs4=genlat_MCLC(a,b,c,α)
# @show b*cos(α)/c + b^2*sin(α)^2/a^2 > 1
# rlatvecs=get_recip_latvecs(mclc_latvecs4)
# ((a,b,c),(α,β,γ))=IBZ.Lattices.get_latparams(rlatvecs)
# @show γ<π/2

# kᵧ < π/2, b*cos(α)/c + b^2*sin(α)^2/a^2 > 1
b=1.1
c=1.2
α=1.1
a=(b*√(c)*sin(α))/√(c - b*cos(α));
mclc_latvecs5=genlat_MCLC(a,b,c,α)
# @show b*cos(α)/c + b^2*sin(α)^2/a^2 == 1
# rlatvecs=get_recip_latvecs(mclc_latvecs5)
# ((a,b,c),(α,β,γ))=IBZ.Lattices.get_latparams(rlatvecs)
# @show γ==π/2

# Four triclinic cases
# kₐ>π/2,kᵦ>π/2,kᵧ>π/2 and kₐ = min(kₐ,kᵦ,kᵧ)
a=1.1
b=1.2
c=1.3
α=3*π/20
β=3.5*π/20
γ=4*π/20
tri_latvecs1=genlat_TRI(a,b,c,α,β,γ)
rlatvecs=get_recip_latvecs(tri_latvecs1)
((a,b,c),(α,β,γ))=IBZ.Lattices.get_latparams(rlatvecs)
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
tri_latvecs2=genlat_TRI(a,b,c,α,β,γ)
rlatvecs=get_recip_latvecs(tri_latvecs2)
((a,b,c),(α,β,γ))=IBZ.Lattices.get_latparams(rlatvecs)
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
tri_latvecs3=genlat_TRI(a,b,c,α,β,γ)
rlatvecs=get_recip_latvecs(tri_latvecs3)
((a,b,c),(α,β,γ))=IBZ.Lattices.get_latparams(rlatvecs)
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
tri_latvecs4=genlat_TRI(a,b,c,α,β,γ)
rlatvecs=get_recip_latvecs(tri_latvecs4)
((a,b,c),(α,β,γ))=IBZ.Lattices.get_latparams(rlatvecs)
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
            @test size(IBZ.Utilities.unique(bzverts),2) == size(bz.points,1)
        end
    end

    @testset "calc_ibz" begin
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
    end
end
