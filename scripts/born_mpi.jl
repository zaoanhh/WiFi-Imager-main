using DrWatson, Revise
includet(srcdir("funcs.jl"))
using CairoMakie, Optimization, OptimizationOptimJL, LinearAlgebra, ThreadPinning, JLD2, ProgressMeter, MPI
BLAS.set_num_threads(1)
if Threads.nthreads() > 1
  if Sys.islinux()
    pinthreads(:cores)
  end
end

function split_count(N::Integer, n::Integer)
  q, r = divrem(N, n)
  return [i <= r ? q + 1 : q for i = 1:n]
end

function main()

  config_path = projectdir("config", "solid_born.jl")
  includet(config_path)

  save_dir = datadir("solid_res", "born_mpi")
  mkpath(save_dir)
  data_folder = datadir("solid_born_jld")
  share_dir = datadir("share_dir", "born_mpi")
  MPI.Init()
  comm = MPI.COMM_WORLD
  rank = MPI.Comm_rank(comm)
  comm_size = MPI.Comm_size(comm)
  root = 0
  if rank == root
    @info "Config Load Finish!"
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
    shared_variables_kwargs = (prefix="solid_shared_variables", verbose=false, tag=false)
    dict_res, _ = produce_or_load(shared_variables_for_born, shared_variables_config, share_dir; filename=hash, shared_variables_kwargs...)
    @info "Shared Variables Load Finish!"

    material_names, material_labels, rectangle_pos = read_label_xlsx(datadir("坐标对应.xlsx"), all_materials_labels)

    counts = split_count(length(material_names), comm_size)
    ids = collect(1:length(material_names))
    ids_vbuf = VBuffer(ids, counts) # VBuffer for scatter
    size_ubuf = UBuffer(counts, 1)

  else
    parameters_all = nothing
    dict_res = nothing
    material_names = nothing
    size_ubuf = UBuffer(nothing)
    ids_vbuf = VBuffer(nothing)
    material_labels = nothing
    rectangle_pos = nothing
  end
  MPI.Barrier(comm)
  local_size = MPI.Scatter(size_ubuf, NTuple{1,Int}, root, comm)
  local_ids = MPI.Scatterv!(ids_vbuf, zeros(Int64, local_size), root, comm)
  material_names = MPI.bcast(material_names, comm)
  dict_res = MPI.bcast(dict_res, comm)
  parameters_all = MPI.bcast(parameters_all, comm)
  material_labels = MPI.bcast(material_labels, comm)
  rectangle_pos = MPI.bcast(rectangle_pos, comm)
  Ei_in_rxs_all = dict_res["Ei_in_rxs_all"]
  Ei_in_domain_all = dict_res["Ei_in_domain_all"]

  if rank == root
    task_num = size(local_ids, 1) * sample_num * length(selected_subcarriers)
    p = Progress(
      task_num,
      dt=0.5,
      barglyphs=BarGlyphs("[=> ]"),
      barlen=50,
      color=:yellow,
    )
  end

  for file_index in local_ids
    mark = 0
    temp_file_name = joinpath(data_folder, "$(material_names[file_index])_winsize=$(window_size)_samnum=$(sample_num).jld2")
    F_all, scatter, files_is_empty = load(temp_file_name, "F_all", "scatter", "files_is_empty")
    GC.gc()
    if sum(files_is_empty) == length(files_is_empty)
      if rank == root
        next!(p; step=sample_num * length(selected_subcarriers))
      end
      continue
    end
    for sample_index in 1:sample_num
      scatter_res = zeros(ComplexF64, size(scatter, 1), size(scatter, 2), length(selected_subcarriers))
      scatter_save_file = joinpath(save_dir, "$(material_names[file_index])_$(sample_index).jld2")
      if mark == 0 && rank == root
        @info "Start Optimization"
        mark = 1
      end
      if isfile(scatter_save_file) == false
        for sub_freq_ids in axes(selected_subcarriers, 1)
          F = F_all[files_is_empty.==false, :, sub_freq_ids, sample_index]
          Ei_in_rxs = Ei_in_rxs_all[sub_freq_ids][files_is_empty.==false, :]
          GS = dict_res["GS_all"][sub_freq_ids][files_is_empty.==false, :]
          Ei_in_domain = Ei_in_domain_all[sub_freq_ids]
          GC.gc()
          scatter_now = born_phaseless(GS=GS, F=F, Ei_in_domain=Ei_in_domain, Ei_in_rxs=Ei_in_rxs, scatter=scatter, parameters=parameters_all[sub_freq_ids], α=0.02)
          scatter_res[:, :, sub_freq_ids] = scatter_now
          if rank == root
            next!(p)
          end
          GC.gc()
        end
        jldsave(scatter_save_file; scatter_res)
      else
        if rank == root
          next!(p; step=length(selected_subcarriers))
        end
      end
    end
  end
  if rank == root
    finish!(p)
  end
  MPI.Barrier(comm)
end

main()
