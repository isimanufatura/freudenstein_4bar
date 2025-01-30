using LinearAlgebra
using Plots

function solve_alpha(l0,l1,l2,l3,θ,ϕ) 

    x = l0 - l1*cos(θ) + l3*cos(ϕ)
    y = - l1*sin(θ) + l3*sin(ϕ)

    α = atan(y , x)
    #α = θ - α

    return α
end

function angs(l0,l1,l2,l3,θ)
    k1 = -2*l1*l3*sin(θ)
    k2 = 2*l3*(l0-(l1*cos(θ)))
    k3 = l0^2 + l1^2 -l2^2 + l3^2 - 2*l0*l1*cos(θ)

    Δ = k1^2 + k2^2 - k3^2
    a = +k1 + sqrt(Δ)
    
    ϕ = (2 * atan(k3 - k2 , a))+π
    α = solve_alpha(l0,l1,l2,l3,θ,ϕ) 
    λ = θ-α
    return ϕ,α
end

function posicao(l0, l1, l2, l3)
    θs = range(0, stop=2π, length=40)
    ϕs = []
    αs = []

    for θ in θs
        phi,alpha = angs(l0, l1, l2, l3, θ)
        push!( ϕs, phi )
        push!( αs, alpha )
    end
    
    return ϕs, αs, θs
end

function A_i(l1,l2,l3,u1,u2,v1)
    E = [0 -1 ; 1 0] 

    A = [l1*E*u1 zeros(2) l2*E*u2;
     1 1 -1]

    v = [(l3*E*v1)-(l2*E*u2);1]
    
    return A,v
end

function unitario(v,θ)
    vi = v*cos(θ)
    vj = v*sin(θ)
    nm = [vi,vj]
    
    n = nm/norm(nm)
    
    return n
end

function velocidades()
    l0 = 80
    l3 = 20
    l2 = 66
    l1 = 56

    ϕs, αs, θs = posicao(l0,l3,l2,l1)

    u1s = []
    u2s = []
    v1s = []
    
    vs = []
    
    tamanho = length(θs)
    
    for (ϕ, α, θ) in zip(ϕs, αs, θs)
        push!(u1s, unitario(l3, ϕ))
        push!(u2s, -unitario(l2, α))
        push!(v1s, unitario(l1, θ))
    end
    
    φv = 25.0 #rad / s
    
    for i in 1:tamanho
        Ai,vi = A_i(l1,l2,l3,u1s[i],u2s[i],v1s[i])

        vel_i = (Ai \ vi) * φv
        push!(vs,vel_i[1])
    end
    
    return vs,θs
end

function main()
  ϕs, αs, θs = posicao(l0,l1,l2,l3)
  vs,θs = velocidades()
  display(plot(rad2deg.(θs), rad2deg.(ϕs), label="φ(θ)", xlabel="θ", ylabel="φ(θ)", title="Variação de φ com θ"))
  display(plot(rad2deg.(θs), rad2deg.(αs), label="α(θ)", xlabel="θ", ylabel="α(θ)", title="Variação de α com θ"))
  display(plot(rad2deg.(θs), vs, label=false, xlabel="θ°", ylabel="dϕ/dt (rad/s)", size=(1000/1.6,550/1.6),grid=false,minorticks=5))
  
end
