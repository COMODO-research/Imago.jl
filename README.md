# Imago.jl
A Julia package of 3D medical image segmentation. 

# Description
This package enables basic semi-automated (but mostly manual) 3D medical image segmentation. The package uses [Comodo.jl](https://github.com/COMODO-research/Comodo.jl) tools, as well as [Makie]() based visualiation and [DICOM.jl](https://github.com/JuliaHealth/DICOM.jl) based DICOM file handling. Image data is rendered in 3D and the user can navigate slices. Next if the image contains intensity variations defining lines, then "contours" can be desired that define the boundary of objects of interest. An exampl in the video below is MRI data of the lower leg. The contours here are the cortical bone/muscle boundary and together the contours for all slices segment the tibia bone. With `Imago` users can create automated contours, but also edit, cut, merge, smooth, and manually draw them. 

![Animation shown bone segmentation process for MRI of the lower leg](https://github.com/COMODO-research/Imago.jl/blob/main/assets/anim/segmentation_tibia.gif)

# Getting started 
This package is currently not registered. Hence to install use: 
```
] add https://github.com/COMODO-research/Imago.jl.git
```
Next one could start by testing an example from the example folder. 

# Documentation
Under construction. 

# License
https://github.com/COMODO-research/Imago.jl/blob/main/LICENSE

