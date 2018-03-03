using PkgName, Base.Test

function ortho(V::AbstractMatrix{T})where{T}
    n = size(V,2)
    t = true
    for i in 1:n-1
        for j in (i+1):n
            if(vecdot(view(V, :, i),view(V, :, j))>(10*eps(T)))
                t = false
                break
            end
        end
    end
    return t
end

@test ortho(hessenberg(sprand(100,100,0.3),30)[2]) == true
