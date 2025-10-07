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

The following is a quick demo using included demo data: 
```julia
using Imago

# Example data 
testCase = 1 
if testCase == 1 # MRI data for a part of the lower leg (excludes knee and foot)
    dcmFolder = getdemodata("MRI_human_lower_leg")
elseif testCase == 2 # CT data Visible human female
    dcmFolder = getdemodata("VHF_Aligned_CT_DICOM")
elseif testCase == 3 # CT data Visible human male
    dcmFolder = getdemodata("VHM_Aligned_CT_DICOM")
end

contoursegment(dcmFolder)
```

# Documentation

The `contoursegment` function is meant to be used to segment contours (e.g. lines) for one particular feature of interest, such as a bone, at a time. The user should complete the contouring process for one such feature, save the result, and restart `contoursegment` to segment a different feature if desired. 

**View navigation**. The left shows a 3D slice view and the right a 2D slice view. The user can rotate/zoom/pan the 3D view and zoom/pan the 2D view. In addition the set sliders in the bottom can be used to alter the slices shown. Only the 3rd slider influences the current slice to be segmented. 

**Alter colorbar limits**. The colorbar limits only influence the visual appearance, not the contouring performance. The user may choose alternative colorbar limits e.g. in order to help improve contrast of features of interest. 

**Sampling contours**. This option is on by default and the program returns to this mode when other functionality is completed or exited. 
Click anywhere in the slice image on the right to attempt to draw a contour close to that click and for the intensity for the location clicked. It is often helpful to attempt to improve the contour by clicking at a region where it departs from the desired location. Sampled contours can be thought of as "draft" contours, i.e. they are not accepted and assigned to the current slice yet. Sampled contours have colours from the`viridis` colormap and have a thing white outline. 

**Cutting draft contours**. Press `c` to trigger cutting mode. Press any non-used key, such as ESC to exit. 
Click once to define a first corner of a cutting window (a preview for the cutting window will now be shown in red), then click a second time to define the second corner of the cutting window. All contour content within the window will be cut away. 

**Accepting contours**. Press `a` to trigger accepting mode. Press any non-used key, such as ESC to exit. 
Accepting contours refers to assigning them to keep for the current slice. Once a contour is accepted it is shown in a solid green in both the 3D and 2D slice views. When there is only one draft contour in the 2D view this contour will now immediately be accepted. If there are more draft contours in the view the user first needs to click to select the draft contour to accept (nearest to click is accepted). If this was the first contour to be accepted it is simply accepted. If the current slice already contains accepted contours the user now needs to A) left click close to an existing contour to merge the accepted contour with, or B) right click anywhere to accept the contour as a new entry (not attached to any previously selected contours). 

**Smooth accepted contours**. Press `s` to trigger smoothing mode. Press any non-used key, such as ESC, to exit. 
Smoothing is only for accepted contours and currently only for closed contours (e.g. loops). If there is only one accepted contour in the current slice then that contour is immediated smoothed using the current smoothing settings. If more than one accepted contour is assigned to the current slice then the user next needs to click close to a contour to smooth. 

**Auto-smooth draft contours** Swith the toggle buttom to turn on auto-smoothing of draft contours. 

**Change smoothing parameter $\lambda$** Edit the value for the smoothing parameter using the text input box on the right. 

**Saving accepted contours**. To save contours the user should first define the file name to same the contours to by specifying the path to an XML file in the text input box above the sliders (Note that depending on your operating system, e.g. for Ubuntu Linux, you may need to install additional software to enable pasting text into this text box). Once the path is defined the user may use the save button to save the contours to the XML file. It is recommended that users always check that the saved file has been created completely before exiting the program. 

**Loading contours**. To load a contour first define the path of the XML file to load using the text box above the sliders. Once defined user can click the load buttom to load the previously created contour set. Once loaded the contours should visualise like accepted contours in green. 

# License
Imago.jl is licensed using the [Apache 2.0 license](https://github.com/COMODO-research/Imago.jl/blob/main/LICENSE). 

