using Plots
using FFTW

# ================= Core types =================
struct ArrayAntenna
    n::Int
    d::Float64      # spacing [m]
    λ::Float64      # wavelength [m]
    I::Vector{Float64}
    β::Vector{Float64}
end

ArrayAntenna(n::Int, d::Float64, λ::Float64) = ArrayAntenna(n, d, λ, zeros(n), zeros(n))

set_currents!(arr::ArrayAntenna, I::Vector{Float64}) = (arr.I .= I; arr)
set_phases!(arr::ArrayAntenna,  β::Vector{Float64}) = (arr.β .= β; arr)

# 스티어링 (broadside 기준): β_n = β0 - n k d sin(θ0) 
# 틀린 수식 수정 필요... β_gap 입력받아서 일정 간격으로 
function steer!(arr::ArrayAntenna, θ0::Float64; β0::Float64=0.0)
    k = 2π / arr.λ
    βinc = -k * arr.d * sin(θ0)
    arr.β .= β0 .+ (0:arr.n-1) .* βinc
    return arr
end

# ================= Binomial taper =================
binomial_weights(N::Integer) = Float64.([binomial(N-1, k) for k in 0:N-1])

function normalize_weights(w::AbstractVector{<:Real}; norm::Symbol=:max)
    w = Float64.(w)
    norm === :max && return w ./ maximum(w)
    norm === :sum && return w ./ sum(w)
    norm === :none && return w
    error("norm must be :max, :sum, or :none")
end

# ================= Chebyshev core =================
# 실수형 T_m(x): |x|≤1 -> cos(m*acos x), |x|>1 -> cosh(m*acosh |x|)에 부호 처리
@inline function chebT_real(m::Int, x::Float64)
    if abs(x) <= 1
        return cos(m * acos(x))
    else
        s = (x < 0 && isodd(m)) ? -1.0 : 1.0
        return s * cosh(m * acosh(abs(x)))
    end
end

# 스펙트럼 샘플 C_k = T_M(β cos(π k / N)),  k=0..N-1
function chebyshev_spectrum(N::Int; sll_dB::Real=30.0)
    @assert N ≥ 1
    if N == 1
        return [1.0]
    end
    R = 10.0^(sll_dB/20)       # 리플비
    M = N - 1
    β = cosh(acosh(R)/M)
    C = Vector{Float64}(undef, N)
    @inbounds for k in 0:N-1
        xk = β * cos(π * k / N)
        C[k+1] = chebT_real(M, xk)
    end
    println(C)
    return C
end

# IFFT(+fftshift)로 소자 가중치 도출
function chebyshev_weights(N::Int; sll_dB::Real=30.0, norm::Symbol=:max)
    C = chebyshev_spectrum(N; sll_dB)
    w = FFTW.fftshift(real(ifft(C)))  # 실수 대칭이 본질, 수치오차만 실수화
    # 수치 대칭 보정 & 미세 음수 클램프
    w = (w .+ reverse(w)) ./ 2
    w .= max.(w, 0.0)
    println("weights = ", normalize_weights(w; norm))
    return normalize_weights(w; norm)
end

# ================= Taper selector =================
"""
    set_taper!(arr; taper=:uniform|:binomial|:chebyshev, norm=:max, sll_dB=30)

taper 선택하여 I(가중치) 적용. Chebyshev는 목표 SLL(dB) 지정.
"""
function set_taper!(arr::ArrayAntenna; taper::Symbol=:uniform, norm::Symbol=:max,
                    sll_dB::Real=30.0, strict_even::Bool=false)
    if taper === :chebyshev && iseven(arr.n)
        if strict_even
            error("Chebyshev taper는 보통 홀수 소자(중앙 소자) 기준으로 설계됩니다. n=$(arr.n)")
        else
            @warn "Chebyshev taper는 보통 홀수 소자 기준. n=$(arr.n) — 결과가 기대와 다를 수 있음."
        end
    end
    w = taper === :uniform   ? ones(arr.n) :
        taper === :binomial  ? binomial_weights(arr.n) :
        taper === :chebyshev ? chebyshev_weights(arr.n; sll_dB, norm=:none) :
        error("unknown taper: $taper")
    set_currents!(arr, normalize_weights(w; norm))
end

# ================= Array Factor (Broadside convention) =================
function array_factor_pattern(arr::ArrayAntenna; θ_range_deg=0.0:0.1:360.0)
    k = 2π / arr.λ
    θdeg = collect(θ_range_deg)
    θ = deg2rad.(θdeg)
    n = collect(0:arr.n-1)
    kd_n = k * arr.d .* n
    w = arr.I .* exp.(1im .* arr.β)
    phase_mat = cos.(θ) * kd_n'        # (L×1)*(1×N) = L×N
    AF = exp.(1im .* phase_mat) * w    # (L×N)*(N) = L
    return θdeg, AF
end

function plot_and_save_patterns(arr::ArrayAntenna;
        θ_range_deg=0.0:0.1:360.0,
        cartesian_path="AF_cart.png",
        polar_path="AF_polar.png",
        dB_floor=-60.0)

    θdeg, AF = array_factor_pattern(arr; θ_range_deg)
    mag = abs.(AF); mag ./= maximum(mag)
    AFdB = 20 .* log10.(mag .+ eps())

    p1 = plot(θdeg, AFdB, xlabel="θ (deg)", ylabel="|AF| (dB)",
              title="Array Factor (broadside)", ylim=(dB_floor, 0),
              legend=false, grid=true)
    savefig(p1, cartesian_path)

    p2 = plot(deg2rad.(θdeg), mag, proj=:polar, legend=false,
              title="Array Factor (polar)")
    savefig(p2, polar_path)

    return cartesian_path, polar_path
end

# ================= (옵션) SLL 측정 헬퍼 =================
# 메인로브 근처 ±Δdeg 제외 후 사이드로브 최대 dB
function measure_SLL_dB(θdeg::Vector{Float64}, AF::Vector{ComplexF64}; exclude_deg::Float64=10.0)
    mag = abs.(AF); mag ./= maximum(mag)
    AFdB = 20 .* log10.(mag .+ eps())
    i0 = argmax(mag)
    θ0 = θdeg[i0]
    mask = trues(length(θdeg))
    for (i,θ) in enumerate(θdeg)
        if abs(θ - θ0) ≤ exclude_deg
            mask[i] = false
        end
    end
    SLL = maximum(AFdB[mask])
    return SLL
end

# ================= Example =================
# 1) Uniform vs Chebyshev(30 dB)
arr = ArrayAntenna(9, 0.5, 1.0)
set_taper!(arr; taper=:uniform, norm=:max)
set_phases!(arr, zeros(arr.n))
θu, AFu = array_factor_pattern(arr)
println("Uniform SLL ≈ ", measure_SLL_dB(θu, AFu), " dB")   # ~ -13 dB 주변 예상
plot_and_save_patterns(arr; cartesian_path="AF_uniform.png", polar_path="AF_uniform_polar.png")

arr2 = ArrayAntenna(9, 0.5, 1.0)
set_taper!(arr2; taper=:chebyshev, sll_dB=30, norm=:max)
set_phases!(arr2, zeros(arr2.n))
θc, AFc = array_factor_pattern(arr2)
println("Chebyshev(30 dB) SLL ≈ ", measure_SLL_dB(θc, AFc), " dB")  # ~ -30 dB 근처
plot_and_save_patterns(arr2; cartesian_path="AF_cheb30.png", polar_path="AF_cheb30_polar.png")
