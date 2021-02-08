quad(A::AbstractMatrix, b::AbstractVector) = dot(b, A * b)
invquad(A::AbstractMatrix, b::AbstractVector) = dot(b, A \ b)