using DrWatson, Revise
@quickactivate "WiFi-Imager"
includet(srcdir("funcs.jl"))
using CairoMakie, JLD2, ProgressMeter, LinearAlgebra, CairoMakie, TSne, MLJ, NearestNeighborModels, DataFrames


config_path = projectdir("config", "solid.jl")
includet(config_path)
@info "Config Load Finish!"

material_names, material_labels, rectangle_pos = read_label_xlsx(datadir("solid_res", "坐标对应.xlsx"), all_materials_labels)
@info "file_namse have been read!"

parameters_all = map(
    frequency -> ConstantParameter(freq=frequency,
        doi_size_x=doi_size,
        doi_size_y=doi_size,
        grid_number_x=grid_number,
        grid_number_y=grid_number,
        txs_pos=Tx_pos,
        rxs_pos=Rx_pos),
    frequencies_all)


labels_gt = Vector{Float64}()
features = Vector{Any}()
for file_index in axes(material_names, 1)
    for sample_index in 1:sample_num
        temp_file_name = joinpath(datadir("solid_res", "som"), "$(material_names[file_index])_$(sample_index).jld2")
        if isfile(temp_file_name)
            temp_data = JLD2.load(temp_file_name, "scatter_res")
            for material_index in axes(material_labels[file_index], 1)
                push!(labels_gt, material_labels[file_index][material_index])
                rectangle_now = zeros(ComplexF64, grid_number, grid_number)
                padding_a_rectangle!(rectangle_now, rectangle_pos[file_index][material_index], 1.0, parameters_all[1])
                features_now = zeros(Float64, length(selected_subcarriers) * 2)
                for freq_index in 1:length(selected_subcarriers)
                    scatter_now = (temp_data[:, :, freq_index] ./ (-1im * parameters_all[freq_index].k₀ / parameters_all[freq_index].η)) .* rectangle_now
                    uni_fea = unique(scatter_now)
                    features_now[2*freq_index-1] = real(uni_fea[findfirst(uni_fea .!= 0)])
                    features_now[2*freq_index] = imag(uni_fea[findfirst(uni_fea .!= 0)])
                end
                push!(features, features_now)
            end
        end
    end
end


@info "Start t-SNE"
X = hcat(features...)'
Y = tsne(X)
let
    f = Figure()
    ax = Axis(f[1, 1])
    unique_labels = unique(labels_gt)
    for i in axes(unique_labels, 1)
        index_now = findall(isequal(unique_labels[i]), labels_gt)
        scatter!(ax, Y[index_now, 1], Y[index_now, 2], label="$i", markersize=5)
    end
    axislegend(ax)
    f
end

X_df = DataFrame(X, :auto)
y_df = coerce(labels_gt, Union{Missing,Multiclass})

KNN = MLJ.@load KNNClassifier
knn = KNN(K=1)
evaluate(knn, X_df, y_df, resampling=CV(nfolds=5),
    measures=[accuracy])


## use MLP

using Lux, ADTypes, Optimisers, Statistics, Zygote, MLUtils, OneHotArrays, Random, Printf

model = Chain(Dense(size(X, 2) => 16, relu),Dense(16 => 16, relu), Dense(16 => size(unique(labels_gt), 1)))

train_split = 0.8
batchsize = 2
(x_train, y_train), (x_test, y_test) = splitobs((X', onehotbatch(labels_gt, unique(labels_gt))); at=train_split, shuffle=true)

train_dataloader = DataLoader(collect.((x_train, y_train)); batchsize, shuffle=true)
test_dataloader = DataLoader(collect.((x_test, y_test)); batchsize, shuffle=false)

loss = CrossEntropyLoss(; logits=Val(true))

function mlp_accuracy(model, ps, st, dataloader)
    total_correct, total = 0, 0
    st = Lux.testmode(st)
    for (x, y) in dataloader
        target_class = onecold(y)
        predicted_class = onecold(Array(first(model(x, ps, st))))
        total_correct += sum(target_class .== predicted_class)
        total += length(target_class)
    end
    return total_correct / total
end

rng = Xoshiro(0)
train_state = Lux.Experimental.TrainState(
    rng, model, Adam(3.0f-4); transform_variables=identity)

### Warmup the model
x_proto = randn(rng, Float32, size(X, 2), 1)
y_proto = onehotbatch([unique(labels_gt)[1]], unique(labels_gt))
Lux.Experimental.compute_gradients(AutoZygote(), loss, (x_proto, y_proto), train_state)

### Lets train the model
nepochs = 200
for epoch in 1:nepochs
    stime = time()
    for (x, y) in train_dataloader
        (gs, _, _, train_state) = Lux.Experimental.single_train_step!(
            AutoZygote(), loss, (x, y), train_state)
    end
    ttime = time() - stime

    tr_acc = mlp_accuracy(
        model, train_state.parameters, train_state.states, train_dataloader) * 100
    te_acc = mlp_accuracy(
        model, train_state.parameters, train_state.states, test_dataloader) * 100
    if epoch % 10 == 0


        @printf "[%2d/%2d] \t Time %.2fs \t Training Accuracy: %.2f%% \t Test Accuracy: \
                 %.2f%%\n" epoch nepochs ttime tr_acc te_acc
    end
end
