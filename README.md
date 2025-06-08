# WiFi-Imager-code
This is the preprocessing code for the paper "The field-based model: a new perspective on RF-based material sensing". We employ the Born approximation to process CSI data, enabling pre-imaging using phaseless data. The image enhancement part is implemented based on a modified U-Net network. Since this portion of the code is highly mature and presents a "two-language" issue with this repository, it has not been open-sourced here.


```bibtex
@article{WiFiImager2025,
  author = "Shang Fei,Jiang Haocheng,Yang Panlong,Yan Dawei,Du Haohua,Li Xiang-Yang",
  title = "The field-based model: a new perspective on RF-based material sensing",
  journal = "SCIENCE CHINA Information Sciences",
  year = "2025",
 pages = "-",
  url = "http://www.sciengine.com/publisher/Science China Press/journal/SCIENCE CHINA Information Sciences///10.1007/s11432-024-4444-8,
  doi = "https://doi.org/10.1007/s11432-024-4444-8"
}
```

Some [documentation](https://zaoanhh.github.io/WiFi-Imager-main/dev/) has been set up for this project. It can be viewed by
running the file `docs/make.jl`, and then launching the generated file
`docs/build/index.html`.
Alternatively, the documentation may be already hosted online.



This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> WiFi-Imager

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "WiFi-Imager"
```
which auto-activate the project and enable local path handling from DrWatson.







