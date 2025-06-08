
module SVDInverse
import FromFile: @from
using LinearAlgebra
using SpecialFunctions
using FFTW
using Optimization
using OptimizationOptimJL
using Zygote
using ForwardDiff
using LsqFit
using DrWatson
using Lux, Statistics
@from "$(srcdir("forward.jl"))" using MoMForward: ConstantParameter, get_Ei_in_rxs, get_Ei_in_domain, get_GS
export get_GD_Z, GD, get_α⁺, get_α⁺_LM, get_V⁺_V⁻, get_shared_variables, get_delta_dat_and_sta_phaseless, get_delta_dat_and_sta, shared_variables_for_save, get_ξ₀, born_phaseless, shared_variables_for_born, get_shared_variables_born

@doc raw"""
    get_GD_Z(config::ConstantParameter)
"""
function get_GD_Z(config::ConstantParameter)
    M = config.grid_number_x
    N = config.grid_number_y
    grid_size_x = config.doi_size_x / config.grid_number_x
    grid_size_y = config.doi_size_y / config.grid_number_y
    x_diff = ((1-M):1:(M-1)) .* grid_size_x
    y_diff = ((1-N):1:(N-1)) .* grid_size_y

    R = sqrt.(transpose(x_diff .^ 2) .+ y_diff .^ 2)
    ZZ = -config.η *
         π *
         (config.grid_radius / 2) *
         besselj(1, config.k₀ * config.grid_radius) .* hankelh1.(0, config.k₀ .* R)
    ZZ[N, M] = -config.η *
               pi *
               (config.grid_radius / 2) *
               hankelh1(1, config.k₀ * config.grid_radius) - 1im * config.η / (config.k₀)

    Z = zeros(ComplexF64, 2 * N - 1, 2 * M - 1)
    Z[1:N, 1:M] = ZZ[N:(2*N-1), M:(2*M-1)]
    Z[(N+1):(2*N-1), (M+1):(2*M-1)] = ZZ[1:(N-1), 1:(M-1)]
    Z[1:N, (M+1):(2*M-1)] = ZZ[N:(2*N-1), 1:(M-1)]
    Z[(N+1):(2*N-1), 1:M] = ZZ[1:(N-1), M:(2*M-1)]

    return Z
end

function get_J_now(N, M, Ni, J)
    J_now_buff = Zygote.Buffer(zeros(eltype(J), 2 * N - 1, 2 * M - 1, Ni))
    for i in 1:2*N-1
        for j in 1:2*M-1
            for k in 1:Ni
                J_now_buff[i, j, k] = convert(eltype(J), 0)
            end
        end
    end
    J_now_buff[1:N, 1:M, :] = reshape(J, N, M, Ni)
    return copy(J_now_buff)
end

@doc raw"""
    GD(J, Z)

We calculate the Green's function of the domain using the inverse Fourier transform. $G_D J$ is `GD(J, Z)`, where `Z` is given by [`get_GD_Z`](@ref).
"""
function GD(J, Z)
    Ni = size(J, 2)
    N, M = Int.((size(Z) .+ 1) ./ 2)
    J_now = zeros(eltype(J), 2 * N - 1, 2 * M - 1, Ni)
    J_now[1:N, 1:M, :] = reshape(J, N, M, Ni)
    # J_now = get_J_now(N, M, Ni, J)
    Z_fft = fft(Z)
    # opa = zeros(eltype(J), N * M, Ni)
    # for tx in 1:Ni
    #     temp_opa = ifft(Z_fft .* fft(J_now[:, :, tx]))
    #     opa[:, tx] = vec(temp_opa[1:N, 1:M])
    # end
    opa = hcat([vec(ifft(Z_fft .* fft(J_now[:, :, tx]))[1:N, 1:M]) for tx in 1:Ni]...)
    return opa
end

function get_α⁺(F, E_i, GS, V⁺;)
    F = F
    E_i = E_i
    GS = GS
    V⁺ = V⁺
    α⁺₀ = ones(ComplexF64, size(V⁺, 2), size(F, 2))
    rosenbrock = (α⁺, F) -> begin
        tmp = E_i .+ GS * V⁺ * α⁺
        obj = F .- conj.(tmp) .* tmp
        return sum(abs2.(obj))
    end
    optf = Optimization.OptimizationFunction(rosenbrock, Optimization.AutoZygote())
    prob = Optimization.OptimizationProblem(optf, α⁺₀, F)

    α⁺ = Optimization.solve(prob, LBFGS()).u
    return α⁺
end

function get_α⁺_LM(F, E_i, GS, V⁺)
    model = (p, α⁺) -> begin
        tmp = E_i .+ GS * V⁺ * α⁺
        return p .- real.(conj.(tmp) .* tmp)
    end

    α⁺₀ = ones(size(V⁺, 2))
    α⁺ = curve_fit(model, F, zeros(size(F)), α⁺₀)
    return α⁺.param
end
function get_V⁺_V⁻(GS, L)
    U, S, V = svd(GS, full=true)
    return V[:, 1:L], V[:, (L+1):end]
end

function get_delta_dat_and_sta_phaseless(
    xi, GS, GDZ, V⁻, F, Ei_in_domain, J⁺, Ei_in_rxs)
    A = V⁻ - Diagonal(xi) * (GD(V⁻, GDZ))
    B = repeat(xi, 1, size(Ei_in_domain, 2)) .* (Ei_in_domain + GD(J⁺, GDZ)) - J⁺
    α⁻ = (A' * A) \ (A' * B)
    C = Ei_in_rxs + GS * (J⁺ + V⁻ * α⁻)
    delta_dat = F - abs2.(C)
    delta_sta = A * α⁻ - B
    return delta_dat, delta_sta
end
function get_delta_dat_and_sta(
    xi, GS, GDZ, V⁻, F, Ei_in_domain, J⁺, Es_in_rxs)
    A = V⁻ - Diagonal(xi) * (GD(V⁻, GDZ))
    B = repeat(xi, 1, size(Ei_in_domain, 2)) .* (Ei_in_domain + GD(J⁺, GDZ)) - J⁺
    α⁻ = (A' * A) \ (A' * B)
    delta_dat = GS * (J⁺ + V⁻ * α⁻) - Es_in_rxs
    delta_sta = A * α⁻ - B
    return delta_dat, delta_sta
end

function get_shared_variables(parameters_all, L)
    num_freqs = length(parameters_all)
    Ei_in_rxs_all = Vector{Array{ComplexF64}}(undef, num_freqs)
    Ei_in_domain_all = Vector{Array{ComplexF64}}(undef, num_freqs)
    GS_all = Vector{Array{ComplexF64}}(undef, num_freqs)
    GDZ_all = Vector{Array{ComplexF64}}(undef, num_freqs)
    V⁺_all = Vector{Array{ComplexF64}}(undef, num_freqs)
    V⁻_all = Vector{Array{ComplexF64}}(undef, num_freqs)
    GDV_neg_all = Vector{Array{ComplexF64}}(undef, num_freqs)
    GSV_neg_all = Vector{Array{ComplexF64}}(undef, num_freqs)
    for sub_freq_ids in axes(parameters_all, 1)
        parameters = parameters_all[sub_freq_ids]
        Ei_in_rxs = get_Ei_in_rxs(parameters)
        Ei_in_domain = get_Ei_in_domain(parameters)
        GS = get_GS(parameters)
        GDZ = get_GD_Z(parameters)
        V⁺, V⁻ = get_V⁺_V⁻(GS, L)
        Ei_in_rxs_all[sub_freq_ids] = Ei_in_rxs
        Ei_in_domain_all[sub_freq_ids] = Ei_in_domain
        GS_all[sub_freq_ids] = GS
        GDZ_all[sub_freq_ids] = GDZ
        V⁺_all[sub_freq_ids] = V⁺
        V⁻_all[sub_freq_ids] = V⁻
        parameters_all[sub_freq_ids] = parameters
        GDV_neg_all[sub_freq_ids] = GD(V⁻, GDZ)
        GSV_neg_all[sub_freq_ids] = GS * V⁻
    end
    return Ei_in_rxs_all,
    Ei_in_domain_all,
    GS_all,
    GDZ_all,
    V⁺_all,
    V⁻_all,
    GDV_neg_all,
    GSV_neg_all
end


function get_shared_variables_born(parameters_all)
    num_freqs = length(parameters_all)
    Ei_in_rxs_all = Vector{Array{ComplexF64}}(undef, num_freqs)
    Ei_in_domain_all = Vector{Array{ComplexF64}}(undef, num_freqs)
    GS_all = Vector{Array{ComplexF64}}(undef, num_freqs)
    GDZ_all = Vector{Array{ComplexF64}}(undef, num_freqs)

    for sub_freq_ids in axes(parameters_all, 1)
        parameters = parameters_all[sub_freq_ids]
        Ei_in_rxs = get_Ei_in_rxs(parameters)
        Ei_in_domain = get_Ei_in_domain(parameters)
        GS = get_GS(parameters)
        GDZ = get_GD_Z(parameters)
        Ei_in_rxs_all[sub_freq_ids] = Ei_in_rxs
        Ei_in_domain_all[sub_freq_ids] = Ei_in_domain
        GS_all[sub_freq_ids] = GS
        GDZ_all[sub_freq_ids] = GDZ
        parameters_all[sub_freq_ids] = parameters

    end
    return (Ei_in_rxs_all,
        Ei_in_domain_all,
        GS_all,
        GDZ_all)

end
"""
    shared_variables_for_save(shared_variables_config)
"""
function shared_variables_for_save(shared_variables_config)
    @unpack parameters_all, L, config_path, selected_subcarriers, frequencies_all, Tx_pos, Rx_pos = shared_variables_config
    Ei_in_rxs_all, Ei_in_domain_all, GS_all, GDZ_all, V⁺_all, V⁻_all, GDV_neg_all, GSV_neg_all = get_shared_variables(
        parameters_all,
        L)
    dict_res = Dict("Ei_in_rxs_all" => Ei_in_rxs_all,
        "Ei_in_domain_all" => [ComplexF64.(a) for a in Ei_in_domain_all],
        "GS_all" => [ComplexF64.(a) for a in GS_all],
        "GDZ_all" => [ComplexF64.(a) for a in GDZ_all],
        "Vpos_all" => [ComplexF64.(a) for a in V⁺_all],
        "Vminus_all" => [ComplexF64.(a) for a in V⁻_all],
        "txs_pos" => Tx_pos,
        "rxs_pos" => Rx_pos,
        "config_path" => config_path,
        "selected_subcarriers" => selected_subcarriers,
        "GDV_neg_all" => GDV_neg_all,
        "GSV_neg_all" => GSV_neg_all,
        "frequencies" => frequencies_all)
    return dict_res
end

function shared_variables_for_born(shared_variables_config)
    @unpack parameters_all, L, config_path, selected_subcarriers, frequencies_all, Tx_pos, Rx_pos = shared_variables_config
    Ei_in_rxs_all, Ei_in_domain_all, GS_all, GDZ_all = get_shared_variables_born(
        parameters_all)
    dict_res = Dict("Ei_in_rxs_all" => Ei_in_rxs_all,
        "Ei_in_domain_all" => [ComplexF64.(a) for a in Ei_in_domain_all],
        "GS_all" => [ComplexF64.(a) for a in GS_all],
        "GDZ_all" => [ComplexF64.(a) for a in GDZ_all],
        "txs_pos" => Tx_pos,
        "rxs_pos" => Rx_pos,
        "config_path" => config_path,
        "selected_subcarriers" => selected_subcarriers,
        "frequencies" => frequencies_all)
    return dict_res
end
function get_ξ₀(scatter, ξ_values, air_permittivity=1.0)
    @assert (size(unique(scatter), 1) - 1) == (size(ξ_values, 1) ÷ 2) "The size of scatter and ξ_values must be the same!"
    material_labels = [element for element in unique(scatter) if element ≠ air_permittivity]
    ξ_array = zeros(ComplexF64, size(vec(scatter), 1), size(unique(scatter), 1) - 1)
    for i in axes(material_labels, 1)
        ξ_array[vec(scatter .== material_labels[i]), i] .= complex(ξ_values[2*i-1],
            ξ_values[2*i])
    end

    return vec(sum(ξ_array, dims=2))

end


function born_phaseless(; GS, F, Ei_in_domain, Ei_in_rxs, scatter, parameters, α=0.02)
    scatter_bce = zeros(size(scatter))
    scatter_bce[scatter.!=1] .= 1.0
    bce = BinaryCrossEntropyLoss()
    rosenbrock_phaseless = (ξ, F) -> begin
        temp = abs2.(GS * Diagonal(ξ) * Ei_in_domain .+ Ei_in_rxs)
        ξ_pwr = abs2.(ξ)
        image = (ξ_pwr .- minimum(ξ_pwr)) ./ (maximum(ξ_pwr) - minimum(ξ_pwr))
        # return sum(abs2.(temp ./ mean(temp) - F ./ mean(F))) + α_phaseless * sum(abs2.(ξ))
        return sum(abs2.(temp ./ mean(temp) - F ./ mean(F))) + α * bce(image, vec(scatter_bce))
    end
    ξ₀ = 0 .- 0.6im .* rand(size(Ei_in_domain, 1))
    optf_phaseless = Optimization.OptimizationFunction(rosenbrock_phaseless, Optimization.AutoZygote())
    prob_phaseless = Optimization.OptimizationProblem(optf_phaseless, ξ₀, F)
    ξ_born_phaseless = reshape(Optimization.solve(prob_phaseless, LBFGS()).u, parameters.grid_number_x, parameters.grid_number_y)
    return real.(ξ_born_phaseless ./ (-1im * parameters.k₀ / parameters.η)) .+ 1.0
end
end