using Kavr
using Test

@testset "roots" begin
    # test Schwarzschild
    @test all(isapprox.(get_roots(27, 0., 0.), (-6, 0, 3, 3), rtol=1e-5))
    @test all(isapprox.(get_roots(0., √27., 0.), (-6, 0, 3, 3), rtol=1e-5))
    # test case 1 or 2
    @test all(isapprox.(get_roots(-1.20503, 11.4, -0.82), (-12.332812618682777, -0.002731857389424548, 2.4271811157108694, 9.908363360361331), rtol=1e-5))
    # test case 3
    @test begin 
        roots = get_roots(2.47908, 2.4, -0.915)
        roots_ans = (-3.820331170598019, 0.07875611823150777, 1.8707875261832558 - 1.8435163296601766im, 1.8707875261832558 + 1.8435163296601766im )
        print(roots)
        (roots[3] ≈ conj(roots[4])) && all(isapprox.(map(x->real(x)+abs(imag(x))im, roots), map(x->real(x)+abs(imag(x))im, roots_ans), rtol=1e-5))
        
    end
    # test for case 4
    @test begin
        roots = get_roots(-0.087975, -1, -1.)
        roots_ans = (-0.3710972685097377 - 0.5479363974641974im, -0.3710972685097377 + 0.5479363974641974im, 0.3710972685097377 - 0.2513305984649188im, 0.3710972685097377 + 0.2513305984649188im)
        (roots[1] ≈ conj(roots[2])) && (roots[3] ≈ conj(roots[4]))all(isapprox.(map(x->real(x)+abs(imag(x))im, roots), map(x->real(x)+abs(imag(x))im, roots_ans), rtol=1e-5))
    end
end
