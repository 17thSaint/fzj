using LinearAlgebra

function lanczos(a,iters)
    alpha = zeros(iters)
    beta = zeros(iters-1)
    v = zeros(size(a,1),iters)
    v[:,1] = rand(size(a,1))
    v[:,1] = v[:,1]./norm(v[:,1])
    w = zeros(size(a,1),iters)
    for i in 1:iters
        if i == 1
            wp = a*v[:,i]
            alpha[i] = transpose(conj(wp)) * v[:,i]
            w[:,i] = wp - alpha[i].*v[:,i]
        else
            beta[i-1] = norm(w[:,i-1])
            v[:,i] = w[:,i-1]./beta[i-1]
            v[:,i] = v[:,i]./norm(v[:,i])
            wp = a*v[:,i]
            alpha[i] = transpose(conj(wp)) * v[:,i]
            w[:,i] = wp - alpha[i].*v[:,i] - beta[i-1].*v[:,i-1]
        end
        #=println("Finished Iteration ",i)
        display(v[:,i])
        println("Alpha: ",alpha[i]," Beta: ", i > 1 ? beta[i-1] : "1st Iteration")=#
    end
    t = Matrix(Tridiagonal(beta,alpha,beta))
    return t,v
end


iters = 3
a = rand(10,10)
t,v = lanczos(a,iters)
display(round.(t - (transpose(v) * a * v),digits=4))

for i in 1:iters
    for j in 1:i
        println("$i,$j: ",dot(v[:,i],v[:,j]))
    end
end































"fin"