using IBZ, Test

@testset "lattices" begin
  @testset "minkowski_reduce" begin
    # Test 1
    basis = [0.997559590093 0.327083693383 0.933708574257;
             0.246898693832 0.853865606330 0.534685117471;
             0.470554239317 0.503012665094 0.122191325919]
    rbasis=minkowski_reduce(basis)
    @test check_reduced(rbasis)

    # Test 2
    basis = [0.566824707611  0.230580598095  0.644795207463;
             0.062164441073  0.069592500822  0.444446612167;
             0.400346380548  0.334526702097  0.184730262940]
    rbasis=minkowski_reduce(basis)
    @test check_reduced(rbasis)

    # Test 3
    basis = [0.965915375496  0.913559885031  0.715882418258;
             0.907651021009  0.074522909849  0.841837872220;
             0.273128578139  0.826489413457  0.040093668849]
    rbasis=minkowski_reduce(basis)
    @test check_reduced(rbasis)

    # Test 4
    basis = [0.747927685471  0.288681880791  0.388321675735;
             0.519739433266  0.818689570053  0.616453247920;
             0.886262688190  0.399796817611  0.970347244023]
    rbasis=minkowski_reduce(basis)
    @test check_reduced(rbasis)

    # Test 5
    basis = [0.367390897375  0.660775660474  0.211656825515;
             0.084605591875  0.481734488916  0.429855314823;
             0.585439077728  0.797785266051  0.636905287099]
    rbasis=minkowski_reduce(basis)
    @test check_reduced(rbasis)

    # Test 6
    basis = [0.480636898755  0.799208677212  0.981910531646;
             0.395947716585  0.894180214969  0.229684979406;
             0.466195536761  0.321402105618  0.558040043342]
    rbasis=minkowski_reduce(basis)
    @test check_reduced(rbasis)

    # Test 7
    basis = [0.073754263634  0.919171496982  0.084716875229;
             0.101486421912  0.827727284628  0.362098268238;
             0.199341503544  0.463651019283  0.902481585522]
    rbasis=minkowski_reduce(basis)
    @test check_reduced(rbasis)

    # Test 8
    basis = [0.005355229278  0.995622271352  0.278599636291;
             0.047892267806  0.312276051592  0.364509680531;
             0.921399596836  0.679640515605  0.122704568964]
    rbasis=minkowski_reduce(basis)
    @test check_reduced(rbasis)

    # Test 9
    basis = [0.347081288518  0.164663213375  0.877540564643;
             0.239622809292  0.463104682188  0.998022806449;
             0.154273976230  0.627982335854  0.966889709417]
    rbasis=minkowski_reduce(basis)
    @test check_reduced(rbasis)

    # Test 10
    basis = [0.342243816073  0.811478344011  0.567158427035;
             0.643460104906  0.993952140553  0.710611242854;
             0.941907678956  0.015738060731  0.598940997350]
    rbasis=minkowski_reduce(basis)
    @test check_reduced(rbasis)

  end
end
