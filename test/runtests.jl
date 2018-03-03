using PkgName
using Base.Test

# write your own tests here
@test "Test for A*V[:,1:m] is almost equal to V[:,1:m+1]*H" begin include("approximate.jl") end
@test "Test for orthogonality" begin include("ortho.jl") end
