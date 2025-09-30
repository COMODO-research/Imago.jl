# Imago.jl
A Julia package of 3D medical image segmentation. 

# Description
This package enables basic semi-automated (but mostly manual) 3D medical image segmentation. The package uses [Comodo.jl](https://github.com/COMODO-research/Comodo.jl) tools, as well as [Makie]() based visualiation and [DICOM.jl](https://github.com/JuliaHealth/DICOM.jl) based DICOM file handling. Image data is rendered in 3D and the user can navigate slices. Next if the image contains intensity variations defining lines, then "contours" can be desired that define the boundary of objects of interest. An exampl in the video below is MRI data of the lower leg. The contours here are the cortical bone/muscle boundary and together the contours for all slices segment the tibia bone. With `Image` users can create automated contours, but also edit, cut, merge, smooth, and manually draw them. 

<video width="100%" controls>
  <source src="https://github.com/COMODO-research/Imago.jl/blob/main/assets/anim/segmentation_turkish_march.mp4" type="video/mp4">
</video>

# Getting started 

This package is currently not registered. Hence to install use: 
```
] add https://github.com/COMODO-research/Imago.jl.git
```

# License
https://github.com/COMODO-research/Imago.jl/blob/main/LICENSE

