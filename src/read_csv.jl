module CSV2CSI
using LinearAlgebra, CSV, DataFrames, JSON
export read_csv, dbinv
# CSI_VAID_SUBCARRIER_INTERVAL = 1

# csi_vaid_subcarrier_index = vcat(collect(7:CSI_VAID_SUBCARRIER_INTERVAL:32), collect(34:CSI_VAID_SUBCARRIER_INTERVAL:59), collect(67:CSI_VAID_SUBCARRIER_INTERVAL:94),
#     collect(96:CSI_VAID_SUBCARRIER_INTERVAL:123))
csi_vaid_subcarrier_index = collect(1:128)
CSI_DATA_COLUMNS = length(csi_vaid_subcarrier_index)

function dbinv(x)
    ret = 10 .^ (x ./ 10)
    return ret
end

function raw2csi(csi_raw_data::Array)
    CSI_DATA_INDEX = size(csi_raw_data, 1)
    csi_data_array = zeros(ComplexF64, CSI_DATA_INDEX, CSI_DATA_COLUMNS)

    # 将CSI_DATA整数两两组合恢复成复数形式
    for row in 1:CSI_DATA_INDEX
        csi_data_array[row, :] .= complex.(csi_raw_data[row][ csi_vaid_subcarrier_index .* 2], csi_raw_data[row][ csi_vaid_subcarrier_index .* 2 .- 1])
    end

    return csi_data_array
end

function read_csv(file_path::String; window_size=100)
    df = CSV.read(file_path, DataFrame)

    if size(df,1)<window_size+1
        return true, ComplexF64[], Float64[], String[]
    end

    raw_data = df[!,"data"]
    raw_timestamp = df[!,"local_timestamp"]
    timestamp = (raw_timestamp .- raw_timestamp[1]) ./ 1e6
    mac_addr = unique(df[!,"mac"])
    row_rssi = df[!,"rssi"]
    row_noise_floor = df[!,"noise_floor"]
    #csi = Array{Any}(undef, length(raw_data))
    csi = Vector{Any}()
    rssi = Vector{Float64}()
    noise_floor = Vector{Float64}()
    for (i, entry) in enumerate(raw_data)
        #csi[i] = JSON.parse(entry)
        temp = JSON.parse(entry)
        if size(temp,1)>255
                push!(csi,temp)
                push!(rssi,row_rssi[i])
                push!(noise_floor,row_noise_floor[i])
        end
    end
    if size(csi,1)<window_size+1
            return true, ComplexF64[],Float64[],String[]
    end
    data = raw2csi(csi) # size is (packet,subcarrier)
    return false, data, timestamp, mac_addr
end
end