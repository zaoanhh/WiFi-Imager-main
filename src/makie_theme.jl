using CairoMakie
my_colors = ["#D43F3AFF", "#EEA236FF", "#5CB85CFF", "#46B8DAFF",
        "#357EBDFF", "#9632B8FF", "#B8B8B8FF"]
function new_cycle_theme()
    # https://nanx.me/ggsci/reference/pal_locuszoom.html
    cycle = Cycle([:color, :linestyle, :marker], covary=true) # alltogether
    my_markers = [:circle, :rect, :utriangle, :dtriangle, :diamond,
        :pentagon, :cross, :xcross]
    my_linestyle = [nothing, :dash, :dot, :dashdot, :dashdotdot]
    return Theme(
        fontsize = 24,
        labelsize = 26,
        linewidth = 2,
        # font="Sarasa Mono SC Nerd",
        colormap = :linear_bmy_10_95_c78_n256,
        palette = (
            color = my_colors,
            marker = my_markers,
            linestyle = my_linestyle,
        ),
        Axis = (
            backgroundcolor= (:white, 0.2),
            xgridstyle = :dash,
            ygridstyle = :dash
        ),
        Lines = (
            cycle= cycle,
        ),
        ScatterLines = (
            cycle = cycle,
        ),
        Scatter = (
            cycle = cycle,
        ),
        Legend = (
            backgroundcolor = (:grey, 0.05),
            framecolor = (:white, 0.2),
            labelsize = 18,
        )
    )
end;
set_theme!(merge(new_cycle_theme(),theme_latexfonts()))
CairoMakie.activate!(type="svg")
