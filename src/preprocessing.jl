### PreProcessing

# The domain is as follows:
# ^
# |
# y
# |
# |---x--->

module PreProcessing
using Random
using DrWatson, XLSX, JSON, Statistics
import FromFile: @from
@from "$(srcdir("forward.jl"))" using MoMForward: ConstantParameter
@from "$(srcdir("read_csv.jl"))" using CSV2CSI
export get_scatterer, get_csi_and_label, read_label_xlsx, get_csiabs2_and_scatter, get_nothing_csiabs2, padding_a_rectangle!

mac_addr_all = ["1a:00:00:00:00:01", "1a:00:00:00:00:02", "1a:00:00:00:00:03", "1a:00:00:00:00:04"]



# function getLinearEquation(p1x, p1y, p2x, p2y)
#     sign = 1
#     a = p2y - p1y
#     if a < 0
#         sign = -1
#         a = sign * a
#     end
#     b = sign * (p1x - p2x)
#     c = sign * (p1y * p2x - p1x * p2y)
#     return [a, b, c]
# end
@doc raw"""
    getLinearEquation(p1x, p1y, p2x, p2y)

Returns the coefficients of a straight line.

# Arguments:

  - `p1x`: The x coordinate of the first point on the line.
  - `p1y`: The y coordinate of the first point on the line.
  - `p2x`: The x coordinate of the second point on the line.
  - `p2y`: The y coordinate of the second point on the line.

# Returns:

  - The coefficients of the straight line, ``[a, b, c]``, where ``ax + by + c = 0``.
"""
function getLinearEquation(x1, y1, x2, y2)
    # 如果直线平行于y轴，即x1==x2，则方程为x=k
    if x1 == x2
        return [1.0, 0.0, -x1]

    else
        # 斜率m=(y2-y1)/(x2-x1)，截距b=y1-m*x1
        m = (y2 - y1) / (x2 - x1)


        # 将方程转换为ax+by=c的形式
        a = m
        b = -1.0
        c = (x1 * y2 - x2 * y1) / (x1 - x2)

        return [a, b, c]
    end
end

@doc raw"""
    get_dis_between_a_point_and_a_line(point, line_point1, line_point2)

Returns the distance between a point and a straight line.
With a line is ``ax + by + c = 0``, the distance is ``|ax + by + c| \sqrt{a^2 + b^2}``.

# Arguments:
- `point`: The coordinates of the point, ``[x, y]``.
- `line_point1`: The coordinates of the first point on the line, ``[x, y]``.
- `line_point2`: The coordinates of the second point on the line, ``[x, y]``.

`line_point1` and `line_point2` is used to calculate straight-line equations, [`getLinearEquation(p1x, p1y, p2x, p2y)`](@ref).

# Returns:
- The distance between the point and the line.
"""
function get_dis_between_a_point_and_a_line(point, line_point1, line_point2)
    a, b, c =
        getLinearEquation(line_point1[1], line_point1[2], line_point2[1], line_point2[2])
    return abs(a * point[1] + b * point[2] + c) / sqrt(a^2 + b^2)
end

function padding_a_rectangle!(scatter_gt, pos, label, parameters::ConstantParameter)
    @assert length(pos) == 8 "The length of pos must be 8!"
    @assert length(label) == 1 "The length of label must be 1!"
    x1, y1, x2, y2, x3, y3, x4, y4 = pos

    dis_14 = sqrt((x1 - x4)^2 + (y1 - y4)^2)
    dis_23 = sqrt((x2 - x3)^2 + (y2 - y3)^2)
    dis_12 = sqrt((x1 - x2)^2 + (y1 - y2)^2)
    dis_34 = sqrt((x3 - x4)^2 + (y3 - y4)^2)

    @assert dis_14 ≈ dis_23 "The input is not a rectangle!"
    @assert dis_12 ≈ dis_34 "The input is not a rectangle!"

    for i in eachindex(parameters.grid_xs)
        for j in eachindex(parameters.grid_ys)
            point_now = [parameters.grid_xs[i], parameters.grid_ys[j]]
            dis_point_to_12 = get_dis_between_a_point_and_a_line(
                point_now,
                [x1, y1],
                [x2, y2],
            )
            dis_point_to_34 = get_dis_between_a_point_and_a_line(
                point_now,
                [x3, y3],
                [x4, y4],
            )

            dis_point_to_14 = get_dis_between_a_point_and_a_line(
                point_now,
                [x1, y1],
                [x4, y4],
            )
            dis_point_to_23 = get_dis_between_a_point_and_a_line(
                point_now,
                [x2, y2],
                [x3, y3],
            )

            if dis_point_to_12 < dis_14 &&
               dis_point_to_34 < dis_14 &&
               dis_point_to_14 < dis_12 &&
               dis_point_to_23 < dis_12
                scatter_gt[j, i] = label
            end
        end
    end

end


function get_scatterer(poses, labels, parameters::ConstantParameter)
    scatter_gt = ones(ComplexF64, length(parameters.grid_ys), length(parameters.grid_xs))
    @assert size(poses, 1) == size(labels, 1) "The number of poses and labels must be the same!"
    for i in eachindex(poses)
        padding_a_rectangle!(scatter_gt, poses[i], labels[i], parameters)
    end
    return scatter_gt
end

"""
    get_csi_and_label(material_names, material_labels, rectangle_pos, folder,parameters::ConstantParameter,selected_freq_index; tx_num=4, rx_num=40,sample_num=10,window_size=100)

Read the label from the xlsx file and get the CSI data. The size of CSI data is (rx_num, tx_num, freq_num, sample_num).
"""
function get_csiabs2_and_scatter(material_name, material_labels, rectangle_pos, folder, parameters::ConstantParameter, selected_freq_index; tx_num=4, rx_num=40, sample_num=10, window_size=100)

    scatter = real.(get_scatterer(rectangle_pos, material_labels, parameters))
    freq_num = length(selected_freq_index)
    if minimum(selected_freq_index) < 0
        selected_freq_index .= (selected_freq_index .+ 67)
    end
    csi_abs2 = zeros(Float64, rx_num, tx_num, freq_num, sample_num)
    files_is_empty = zeros(Bool, rx_num)

    for rx_index in 1:rx_num
        for tx_index in 1:tx_num
            file_now = joinpath(folder, "$(material_name)_tx$(tx_index)_rx$(rx_index).csv")
            empty_mark, csi_data, timestamp, mac_addr = CSV2CSI.read_csv(file_now) # the size of csi_data is (packet,subcarrier)
            files_is_empty[rx_index] = empty_mark

            if !empty_mark
                while size(csi_data, 1) < window_size * sample_num
                    csi_data = repeat(csi_data, 2, 1)
                end
                for sample_index in 1:sample_num
                    csi_abs2[rx_index, findfirst(isequal(mac_addr[1]), mac_addr_all), :, sample_index] = mean(abs2.(csi_data[(sample_index-1)*window_size+1:sample_index*window_size, selected_freq_index]), dims=1)
                end
            end
        end
    end
    return csi_abs2, scatter, files_is_empty
end
function get_nothing_csiabs2(folder, selected_freq_index; tx_num=4, rx_num=40, noting_file_name="test")

    freq_num = length(selected_freq_index)
    if minimum(selected_freq_index) < 0
        selected_freq_index .= (selected_freq_index .+ 67)
    end
    noting_csi_abs2 = zeros(Float64, rx_num, tx_num, freq_num)
    files_is_empty = zeros(Bool, rx_num, tx_num)

    for rx_index in 1:rx_num
        for tx_index in 1:tx_num
            file_now = joinpath(folder, "$(noting_file_name)_tx$(tx_index)_rx$(rx_index).csv")
            empty_mark, csi_data, timestamp, mac_addr = CSV2CSI.read_csv(file_now) # the size of csi_data is (packet,subcarrier)
            files_is_empty[rx_index, tx_index] = empty_mark

            if !empty_mark
                noting_csi_abs2[rx_index, findfirst(isequal(mac_addr[1]), mac_addr_all), :] = mean(abs2.(csi_data[:, selected_freq_index]), dims=1)
            end
        end
    end
    return noting_csi_abs2
end

function read_label_xlsx(filename, all_materials_labels)
    workbook = XLSX.openxlsx(filename)
    sheet = workbook[1]
    material_names = Vector{Any}()
    material_labels = Vector{Any}()
    rectangle_pos = Vector{Any}()
    for r in XLSX.eachrow(sheet)
        no_empty_cells = [r[i] for i in 1:10 if !ismissing(r[i])]
        push!(material_names, no_empty_cells[1])
        push!(material_labels, [all_materials_labels[key_now] for key_now in split(no_empty_cells[1], "-")[1:end-1]])
        push!(rectangle_pos, [(parse.(Float64, split(x, ","))) ./ 100 for x in no_empty_cells[2:end]])
    end
    close(workbook)
    return material_names, material_labels, rectangle_pos
end
end

