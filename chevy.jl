using LinearAlgebra
using Plots

# ---------------- Chebyshev 기본 ----------------
# T_m(x) in monomial basis (coeffs of x^0..x^m)
function chebT_coeffs(m::Int; T=Float64)
    m == 0 && return T[1]
    m == 1 && return T[0, 1]              # x
    Tm2 = T[1]                             # T0
    Tm1 = T[0, 1]                          # T1
    for n in 2:m
        two_x_Tm1 = vcat(zero(T), 2*T(1) .* Tm1)  # 2x*T_{n-1}
        L = max(length(two_x_Tm1), length(Tm2))
        Tn = zeros(T, L)
        Tn[1:length(two_x_Tm1)] .+= two_x_Tm1
        Tn[1:length(Tm2)]      .-= Tm2
        Tm2, Tm1 = Tm1, Tn
    end
    return Tm1
end

# p(x) → p(a*x + b)  (both in monomial coeffs)
function poly_affine_subst(p::AbstractVector{T}, a::T, b::T) where {T<:Number}
    deg = length(p) - 1
    out = zeros(T, deg + 1)
    for k in 0:deg
        ck = p[k+1]; iszero(ck) && continue
        for i in 0:k   # (a x + b)^k
            out[i+1] += ck * T(binomial(k,i)) * (a^i) * (b^(k-i))
        end
    end
    return out
end

# Solve A:  T_n(a x + b) = Σ_{k=0}^n A_k T_k(x)
function cheb_affine_expand(n::Int, a, b; T=Float64)
    rhs   = poly_affine_subst(chebT_coeffs(n; T=T), T(a), T(b))
    basis = [chebT_coeffs(k; T=T) for k in 0:n]
    maxdeg = n + 1
    M = zeros(T, maxdeg, n+1)
    for k in 0:n
        tk = basis[k+1]
        M[1:length(tk), k+1] .= tk
    end
    rhs = vcat(rhs, zeros(T, maxdeg - length(rhs)))
    return M \ rhs                      # A[1]=A0, ..., A[n+1]=An
end

# ---------------- a,b 추출 (슬라이드의 엔드포인트 매칭) ----------------
# x ∈ [x_min, 1]  ↔  x' ∈ [-1, x0'],   x_min = cos(2π d/λ),  x0' = cosh( acosh(R)/n )
function linear_shift_params(n::Int; d_over_lambda::Real, R::Real, T=Float64)
    x_min = cos(2*T(pi) * T(d_over_lambda))
    x0p   = cosh(acosh(T(R)) / T(n))
    a = (x0p + 1) / (1 - x_min)
    b = x0p - a
    return a, b
end

# ---------------- 편의: Dolph-Chebyshev 가중치 A_k ----------------
# N = 2n+1 (odd, symmetric)
function dolph_cheb_A(N::Int; d_over_lambda::Real, R::Real, T=Float64)
    @assert isodd(N) "Use odd N (symmetric), N=2n+1."
    n = (N-1) ÷ 2
    a,b = linear_shift_params(n; d_over_lambda=d_over_lambda, R=R, T=T)
    A = cheb_affine_expand(n, a, b; T=T)           # A0..An
    return A, a, b, n
end

# ---------------- 전류(중심~끝)로 변환 ----------------
# AF(ψ) = A0 + Σ_{k=1}^n A_k cos(kψ) = I0 + 2 Σ_{m=1}^n I_m cos(mψ)
# ⇒ I0 = A0,  I_m = A_m/2 (m≥1)
function currents_from_A(A::AbstractVector)
    I = similar(A)
    I[1] = A[1]
    for m in 2:length(A)
        I[m] = A[m]/2
    end
    return I
end

# ---------------- 패턴 계산 ----------------
# ψ = 2π (d/λ) cos(φ),  φ ∈ [0, π]
function array_factor_from_currents(I::AbstractVector; d_over_lambda::Real,
                                    phi = range(0, stop=pi, length=1801))
    n = length(I)-1
    ψ = 2pi*d_over_lambda .* cos.(phi)
    AF = fill(I[1], length(phi))
    for m in 1:n
        AF .+= 2I[m+1] .* cos.(m .* ψ)
    end
    AFamp = abs.(AF)
    AFamp ./= maximum(AFamp)         # 0 dB 정규화
    AFdB = 20 .* log10.(AFamp .+ 1e-12)
    return phi, AFamp, AFdB
end

# ---------------- 플로팅 ----------------
function plot_patterns(I; d_over_lambda::Real, title_str="")
    phi, AFamp, AFdB = array_factor_from_currents(I; d_over_lambda=d_over_lambda)

    # Cartesian (dB)
    p1 = plot(rad2deg.(phi), AFdB,
              xlabel="ϕ (deg)", ylabel="Normalized AF (dB)",
              ylim=(-60, 0), xlim=(0,180), legend=false, grid=true,
              title=isempty(title_str) ? "Array Pattern (Cartesian, dB)" : title_str)
    savefig("cheby_c.png")

    # Polar (linear) — 대칭 반사로 0..2π 채움
    θ = vcat(phi, 2pi .- reverse(phi))
    r = vcat(AFamp, reverse(AFamp))
    p2 = plot(θ, r, proj=:polar, legend=false, title="Array Pattern (Polar, linear)")
    savefig("cheby_p.png")
    return p1, p2
end

# ===================== 예제 =====================
# N=9(=2*4+1), d/λ=0.5, SLL=-30 dB → R=10^(30/20)
N = 9
d_over_lambda = 0.5
R = 10^(10/20)           # 31.6227766017...

A, a, b, n = dolph_cheb_A(N; d_over_lambda=d_over_lambda, R=R)
I = currents_from_A(A)   # 중심~끝 전류 (대칭이며 실제 소자에선 ±m로 복제)

@show a b A I

p1, p2 = plot_patterns(I; d_over_lambda=d_over_lambda,
                       title_str="Dolph-Chebyshev (N=$N, d/λ=$(d_over_lambda), SLL=-30 dB)")
display(p1); display(p2)
