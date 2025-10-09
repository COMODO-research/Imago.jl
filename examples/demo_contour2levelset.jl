using Comodo
using Comodo.GLMakie
using Imago
using Geogram
using FileIO

GLMakie.closeall()

# Control parameters 
pointSpacing = 3.0 # Desired OBJ output mesh point spacing
saveObj = false # Turn on or off saving of the OBJ file 
testCase = 6 # Selects a data case below

# Files names to import 
segmentationPath = joinpath(imagodir(), "assets", "segmentation")
if testCase == 1 
    fileNameContours = "VHF_right_femur.xml"
elseif testCase == 2 
    fileNameContours = "VHF_right_tibia.xml"
elseif testCase == 3 
    fileNameContours = "VHF_right_patella.xml"
elseif testCase == 4 
    fileNameContours = "VHF_right_fibula.xml"
elseif testCase == 5 
    fileNameContours = "VHF_right_skin_leg.xml"
elseif testCase == 6 
    fileNameContours = "VHF_right_talus.xml"
end

# File name for saving OBJ file 
fileName_part, _ = splitext(fileNameContours)
fileName_OBJ = joinpath(segmentationPath,fileName_part * ".obj")

# Load contours 
contourSets, voxelSize = load_contours(joinpath(segmentationPath, fileNameContours)) 

# Create levelset image from contours 
M, xRange, yRange, zRange = contour2levelset(contourSets, voxelSize)

# Construct surface from levelset image 
F, V = getisosurface(M; x = xRange, y = yRange, z = zRange, level = 0.0, cap = true, padValue=1e8)      

# Remesh and resample surface using Geogram
numPoints = spacing2numvertices(F, V, pointSpacing)
Fg, Vg = ggremesh(F, V; nb_pts=numPoints)

# Saving OBJ file 
if saveObj == true
    M_export = GeometryBasics.Mesh(Vg,Fg)
    save(fileName_OBJ, M_export)
end

################################################################################
# Visualization
xMin, xMax, yMin, yMax, zMin, zMax = get_contour_extrema(contourSets) # To set axis limits
siz = size(M) # To set slider limits for level set image 
indStart = ceil(Int, siz[3]/2)

fig = Figure(size=(1200,800))
ax1 = AxisGeom(fig[1,1], limits=(xMin, xMax, yMin, yMax, zMin, zMax), title="Level set and reconstructed surface")
hs1 = heatmap!(ax1, xRange, yRange, M[:,:,indStart], colormap = :bluesreds, colorrange=(-1.0,1.0))
Makie.translate!(hs1, 0.0, 0.0, zMin+(indStart-1)*voxelSize[3])
hp1 = meshplot!(ax1, Fg, Vg, color=(:white,0.5), transparency=true, strokewidth=0.0)

Colorbar(fig[:,2], hs1)

# Slider for viewing levelset image 
stepRange = 1:siz[3]
hSlider = Slider(fig[2, 1], range = stepRange, startvalue = indStart, linewidth=30)
on(hSlider.value) do iSlice    
    hs1[3] =  M[:,:,iSlice]
    Makie.translate!(hs1, 0.0, 0.0, zMin+(iSlice-1)*voxelSize[3])
end

ax2 = AxisGeom(fig[1,3], limits=(xMin, xMax, yMin, yMax, zMin, zMax), title="Surface mesh")
hp2 = meshplot!(ax2, Fg, Vg, color=:white, strokewidth=0.5, strokecolor=:black)

# Slider for mesh density variation
# stepRange = range(pointSpacing/2.0,pointSpacing*2.0,20)
# hSlider2 = Slider(fig[2, 3], range = stepRange, startvalue = pointSpacing, linewidth=30)

# on(hSlider2.value) do pointSpacing    
#     numPoints = spacing2numvertices(F,V,pointSpacing)
#     Fg,Vg = ggremesh(F,V; nb_pts=numPoints)
#     hp2[1] = GeometryBasics.Mesh(Vg,Fg)    
#     hp1[1] = GeometryBasics.Mesh(Vg,Fg)    
# end

screen = display(GLMakie.Screen(), fig)