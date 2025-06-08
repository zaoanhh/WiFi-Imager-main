doi_size = 1.0
grid_number = 40
tx_xs = [0.6, 0.0, -0.6, 0.0]
tx_ys = [0.0, -0.6, 0.0, 0.6]
all_esp_num = 44
cen_freq = 2.462e9 #channel 11
delta_f = 312.5e3 # The frequency interval between two adjacent sub-carriers is 312.5kHz
air_permittivity = 1.0 + 0.0im
dis_ant = 0.085
Ni = length(tx_xs)
tx_pp = [0]
tx_scale = Ni ÷ 4
Ns = all_esp_num - Ni
L = 1
rx_xs = vcat(repeat(tx_xs[(1*tx_scale):(1*tx_scale)], div(Ns, 4)),
    [x for x in div(all_esp_num - 4, 8):-1:div(-(all_esp_num - 4), 8) if x ∉ tx_pp] .*
    dis_ant,
    repeat(tx_xs[(3*tx_scale):(3*tx_scale)], div(Ns, 4)),
    [x for x in div(-(all_esp_num - 4), 8):1:div((all_esp_num - 4), 8) if x ∉ tx_pp] .*
    dis_ant)
rx_ys = vcat(
    [x
     for x in div((all_esp_num - 4), 8):-1:div(-(all_esp_num - 4), 8) if x ∉ tx_pp] .*
    dis_ant,
    repeat(tx_ys[(2*tx_scale):(2*tx_scale)], div(Ns, 4)),
    [x for x in div(-(all_esp_num - 4), 8):1:div((all_esp_num - 4), 8) if x ∉ tx_pp] .*
    dis_ant,
    repeat(tx_xs[(4*tx_scale):(4*tx_scale)], div(Ns, 4)))
Rx_pos = [rx_xs'; rx_ys']
Tx_pos = [tx_xs'; tx_ys']
selected_subcarriers = [-58, -54, -50, -46, -42, -38, -34, -30, -26, -22, -18, -14, -10, -6, -2, 2, 6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58][1:6:end]
freq_num = length(selected_subcarriers)
frequencies_all = cen_freq .+ selected_subcarriers .* delta_f
all_materials_labels = Dict("rectangle" => 2.0, "wood" => 2.0, "glass" => 3.0, "leather" => 4.0)
window_size = 200
sample_num = 5