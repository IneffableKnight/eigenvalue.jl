function eigen(A::AbstractMatrix{T},m::Int,tol::Float64)where{T}

    H, V = hessenberg(A,m)

    Vm = view(V, :, 1:(m))
    Hm = view(H, 1:(m), :)

    vec = eig(Hm)
    values = vec[1]
    vectors = vec[2] #coloumn eigen vector

    eigenvectors = zeros(eltype(vectors),size(A,1),m)
    identity = falses(m)

    for i in 1:m
        A_mul_B!(view(eigenvectors, :, i),Vm,view(vectors, :, i))
        Res = H[m+1,m]*eigenvectors[m,i]
        if abs(Res) < tol
            identity[i] = true
        end
    end
    return values[identity] , eigenvectors[:,identity]
end
