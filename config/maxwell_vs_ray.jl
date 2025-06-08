if !isdefined(Main, :Drwatson)
    using DrWatson
end

freq = 5e9
grids_y = 100
grids_x = 30
dis_target_rx = [0.1,0.2,0.3,0.4,0.5]
doi_size_x = 0.3
doi_size_y = 1.0
tx_num = 1
rx_num = 3
dis_antenna = 0.03

tag_size_x = 0.1
tx_x = -3.1
configs = @strdict freq grids_x grids_y doi_size_x doi_size_y tx_num rx_num dis_antenna tag_size_x  dis_target_rx tx_x