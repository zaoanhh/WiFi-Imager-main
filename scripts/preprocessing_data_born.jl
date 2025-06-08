using DrWatson, Revise
@quickactivate "WiFi-Imager"
includet(srcdir("funcs.jl"))

using JLD2, ProgressMeter

config_path = projectdir("config", "solid_born.jl")
includet(config_path)
data_folder = datadir("solid_born_jld")
mkpath(data_folder)
parameters_all = map(
    frequency -> ConstantParameter(freq=frequency,
        doi_size_x=doi_size,
        doi_size_y=doi_size,
        grid_number_x=grid_number,
        grid_number_y=grid_number,
        txs_pos=Tx_pos,
        rxs_pos=Rx_pos),
    frequencies_all)
material_names, material_labels, rectangle_pos = read_label_xlsx(datadir("坐标对应.xlsx"), all_materials_labels)

# noting_csi_abs2_all = get_nothing_csiabs2(datadir("guti"), selected_subcarriers .+ 67; tx_num=4, rx_num=40, noting_file_name="test")

# jldsave(joinpath(data_folder, "noting_csi_abs2_all.jld2"); noting_csi_abs2_all)

p = Progress(
    size(material_names, 1),
    dt=0.5,
    barglyphs=BarGlyphs("[=> ]"),
    barlen=50,
    color=:yellow,
)
Threads.@threads for file_index in axes(material_names, 1)
    temp_file_name = joinpath(data_folder, "$(material_names[file_index])_winsize=$(window_size)_samnum=$(sample_num).jld2")
    next!(p)
    if isfile(temp_file_name)
        continue
    else
        F_all, scatter, files_is_empty = get_csiabs2_and_scatter(material_names[file_index], material_labels[file_index], rectangle_pos[file_index], datadir("guti"), parameters_all[1], selected_subcarriers .+ 67; tx_num=Ni, rx_num=Ns, sample_num=sample_num, window_size=window_size)
        jldsave(temp_file_name; F_all, scatter, files_is_empty)
    end
end
finish!(p)