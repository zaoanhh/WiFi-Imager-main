using DrWatson, Revise
@quickactivate "WiFi-Imager"
includet(srcdir("funcs.jl"))
using CairoMakie, Optimization, OptimizationOptimJL, LinearAlgebra, ThreadPinning, JLD2, ProgressMeter

if Threads.nthreads() > 1
    BLAS.set_num_threads(1)
    if Sys.islinux()
        pinthreads(:cores)
    end
end

function main()


    config_path = projectdir("config", "solid_born.jl")
    includet(config_path)
    @info "Config Load Finish!"

    save_dir = datadir("solid_res", "born")
    mkpath(save_dir)
    data_folder = datadir("solid_born_jld")
    share_dir = datadir("share_dir", "born")

    parameters_all = map(
        frequency -> ConstantParameter(freq=frequency,
            doi_size_x=doi_size,
            doi_size_y=doi_size,
            grid_number_x=grid_number,
            grid_number_y=grid_number,
            txs_pos=Tx_pos,
            rxs_pos=Rx_pos),
        frequencies_all)

    shared_variables_config = @strdict(parameters_all, L, config_path, selected_subcarriers, frequencies_all, Tx_pos, Rx_pos)
    dict_res, _ = produce_or_load(shared_variables_for_born, shared_variables_config, share_dir)
    @info "Shared Variables Load Finish!"

    material_names, material_labels, rectangle_pos = read_label_xlsx(datadir("solid_res", "坐标对应.xlsx"), all_materials_labels)



    Ei_in_rxs_all = dict_res["Ei_in_rxs_all"]
    Ei_in_domain_all = dict_res["Ei_in_domain_all"]


    task_num = length(material_names) * sample_num * length(selected_subcarriers)
    @info "Total Task Number: $task_num"

    @info "Start Optimization"

    p = Progress(
        task_num,
        dt=0.5,
        barglyphs=BarGlyphs("[=> ]"),
        barlen=50,
        color=:yellow,
    )
    Threads.@threads for file_index in axes(material_names, 1)
        # F_all, scatter, files_is_empty = get_csiabs2_and_scatter(material_names[file_index], material_labels[file_index], rectangle_pos[file_index], datadir("guti"), parameters_all[1], selected_subcarriers .+ 67; tx_num=Ni, rx_num=Ns, sample_num=sample_num, window_size=window_size)

        temp_file_name = joinpath(data_folder, "$(material_names[file_index])_winsize=$(window_size)_samnum=$(sample_num).jld2")
        F_all, scatter, files_is_empty = load(temp_file_name, "F_all", "scatter", "files_is_empty")
        if sum(files_is_empty) == length(files_is_empty)
            next!(p; step=sample_num * length(selected_subcarriers))
            continue
        end
        for sample_index in 1:sample_num
            scatter_res = zeros(ComplexF64, size(scatter, 1), size(scatter, 2), length(selected_subcarriers))
            scatter_save_file = joinpath(save_dir, "$(material_names[file_index])_$(sample_index).jld2")
            if isfile(scatter_save_file) == false
                for sub_freq_ids in axes(selected_subcarriers, 1)
                    F = F_all[files_is_empty.==false, :, sub_freq_ids, sample_index]
                    Ei_in_rxs = Ei_in_rxs_all[sub_freq_ids][files_is_empty.==false, :]
                    GS = dict_res["GS_all"][sub_freq_ids][files_is_empty.==false, :]
                    Ei_in_domain = Ei_in_domain_all[sub_freq_ids]

                    scatter_now = born_phaseless(GS=GS, F=F, Ei_in_domain=Ei_in_domain, Ei_in_rxs=Ei_in_rxs, scatter=scatter, parameters=parameters_all[sub_freq_ids], α=0.02)
                    scatter_res[:, :, sub_freq_ids] = scatter_now
                    next!(p)
                end
                jldsave(scatter_save_file; scatter_res)
            else
                next!(p; step=length(selected_subcarriers))
            end
        end
    end
    finish!(p)
    @info "Optimization Finish!"
end

main()