using PkgName, Base.Test

function equal(A::AbstractMatrix{T},V::AbstractMatrix{T},H::AbstractMatrix{T},m::Int) where{T}
    x = norm(A*V[:,1:m]-V*H)
    y = 10*eps(T)
    if x < y
        return true
    else
        return false
    end
end

A = sprand(100,100,0.3)
H, V = hessenberg(A,30)
@test equal(A,V,H,30) == true
