"""
    imagodir()

# Description 

This function simply returns the string for the Imago path. This is helpful for instance to load items, such as meshes, from the `assets` folder. 
"""
function imagodir()
    pkgdir(@__MODULE__)
end

function dicomdir2dicomvec(dcmFolder; extension_exclude = [".txt"])
    fileSet = get_image_file_set(dcmFolder; extension_exclude = extension_exclude)        
    
    # Read in first 
    dcm = dcm_parse(joinpath(dcmFolder,fileSet[1]))
    dicomData = Dict{Int, typeof(dcm)}()
    for f in fileSet # for each file 
        dcm = dcm_parse(dcmFolder*"/"*f)
        dicomData[dcm.InstanceNumber] = dcm
    end
    return dicomData
end

function get_image_file_set(imageFolder; extension_exclude = [".txt"])        
    fileNameSet = readdir(imageFolder)
    imageNameSet = Vector{String}()
    for fileName in fileNameSet
        _, fileExtension = splitext(fileName)
        if !in(lowercase(fileExtension),extension_exclude)
            push!(imageNameSet, fileName)
        end
    end
    return imageNameSet
end

function getvalrange(dicomData)    
    sliceKeySet = sort(collect(keys(dicomData)))    
    maxVal = -Inf
    minVal =  Inf
    for sliceKey in sliceKeySet        
        m = dicomData[sliceKey].PixelData
        maxVal = max(maxVal,maximum(m))
        minVal = min(minVal,minimum(m))
    end
    return maxVal, minVal
end

function getslice(dicomData,ind,dir=1)
    sliceKeySet = sort(collect(keys(dicomData)))
    dcm = dicomData[sliceKeySet[1]]        
    T = eltype(dcm.PixelData)
    siz = size(dcm.PixelData)
    n = length(dicomData)
    if dir == 1 
        A = Matrix{T}(undef,n,siz[2])
    else
        A = Matrix{T}(undef,siz[1],n)
    end
    for (iSlice, sliceKey) in enumerate(sliceKeySet)
        if dir == 1 
            A[iSlice,:] = dicomData[sliceKey].PixelData[ind,:]
        elseif dir == 2
            A[:,iSlice] = dicomData[sliceKey].PixelData[:,ind]
        end
    end
    return A
end

function getdemodata(dataName="MRI_human_lower_leg")    
    if dataName == "MRI_human_lower_leg"            
        dataDir = joinpath(imagodir(), "assets", "data", dataName)        
        dataDir_unzip = joinpath(dataDir,"data")
    elseif dataName == "VHF_Aligned_CT_DICOM"
        dataDir = joinpath(imagodir(), "assets", "data", "Visible_human")
        dataDir_unzip = joinpath(dataDir, "VHF_Aligned_CT_DICOM")        
    elseif dataName == "VHM_Aligned_CT_DICOM"
        dataDir = joinpath(imagodir(), "assets", "data", "Visible_human")
        dataDir_unzip = joinpath(dataDir, "VHM_Aligned_CT_DICOM")        
    end
    if isdir(dataDir_unzip)
        return dataDir_unzip        
    end
end

function xy2ij_click(xy_click, voxelSize)
    return ceil.(Int, xy_click./voxelSize[1:2])
end

function contour_at_click(xy_click, voxelSize, raw_image; c = NaN)
    ij_click = xy2ij_click(xy_click, voxelSize)
    if isnan(c)
        global contourLevel = raw_image[ij_click[1],ij_click[2]]
    else
        global contourLevel = c
    end

    c = Contour.contours(x[1:end-1],y[1:end-1],raw_image,[contourLevel])  
    xk = Vector{Vector{Float64}}()              
    yk = Vector{Vector{Float64}}()      
    dp = 1.0
    q = 1
    for (i,l) in enumerate(Contour.lines(Contour.levels(c)[1]))
        xc, yc = Contour.coordinates(l)    
        d = minimum(sqrt.((xc .- xy_click[1]).^2 .+ (yc .- xy_click[2]).^2))                            
        if i == 1                
            xk = xc
            yk = yc
            dp = d
        else
            if d < dp 
                xk = xc
                yk = yc
                q = i
                dp = d
            end
        end                                              
    end      
    V_contour_raw = [Point{2,Float64}(xk[i],yk[i]) for i in eachindex(xk)]
    V_contour_raw = V_contour_raw[1:end-1]
    numDigitsMerge = 3-mag(mean(voxelSize))

    # Create rounded coordinates to help obtain unique set
    # Note -0.0+0.0 = 0.0 so addition of zero points helps avoid 0.0 and -0.0 being seen as unique
    VR = [round.(v,digits = numDigitsMerge)+Point{2,Float64}(0.0,0.0) for v in V_contour_raw]

    # Get unique indices and reverse for rounded vertices
    _,indUnique,indMap = gunique(VR; return_index=Val(true), return_inverse=Val(true),sort_entries=false)
    V_contour_raw = V_contour_raw[sort(indUnique)] # The unique node set

    if toggle.active[]
        V_contour_raw = smooth_contour(V_contour_raw; λ = λ)
    end
    return [Point{3,Float64}(v[1],v[2],dicomData[sliceIndex].SliceLocation) for v in V_contour_raw]
end

function smooth_contour(V_contour_raw; λ = 1.0)
    if length(V_contour_raw)>4
        L = curve_length(V_contour_raw; close_loop=true)
        S = fit(BSplineOrder(4), L[1:end-1], V_contour_raw, λ, BSplineKit.Periodic(maximum(L)))
        L_fit = range(0.0,maximum(L),length(L))
        V_fit = S.(L_fit[1:end-1])
        V_fit_even = evenly_space(V_fit, min(voxelSize[1], voxelSize[2]); close_loop = true, niter=10, spline_order = 4) 
        return V_fit_even
    else
        # V_fit_even = evenly_space(V_contour_raw, min(voxelSize[1], voxelSize[2]); close_loop = true, niter=10, spline_order = 2) 
        return V_contour_raw
    end
end   

function makeplottable(contourSets)
    if make_2D
        n = 2
        V_nan = Point{2,Float64}(NaN, NaN)
    else
        n = 3
        V_nan = Point{3,Float64}(NaN, NaN, NaN)
    end

    V_plot = Vector{Point{n,Float64}}()
    C_plot = Vector{Int}()
    for contourSet in contourSets 
        for (i,V) in enumerate(contourSet)
            if make_2D
                V = [Point{2,Float64}(v[1], v[2]) for v in V]
            end
            push!(V_plot, V_nan)
            append!(V_plot, V)
            push!(V_plot, V_nan)
            append!(C_plot, fill(i,length(V)+2))
        end    
    end
    return V_plot, C_plot
end

function update_draft_plot(h, contourSet; use_color=false)
    V_nan = Point{3,Float64}(NaN, NaN, NaN)
    V_plot = Vector{Point{3,Float64}}()
    C_plot = Vector{Int}()    
    for (i,V) in enumerate(contourSet)                    
        append!(V_plot, V)
        push!(V_plot, V_nan)
        append!(C_plot, fill(i,length(V)+1))
    end    
    
    h[1] = V_plot
    if use_color == true 
        h.color = C_plot
    end
end

function update_accepted_plot(h1, h2, contourSets)
    V_nan = Point{3,Float64}(NaN, NaN, NaN)
    V_plot3D = Vector{Point{3,Float64}}()    
    V_plot3D_ax1 = Vector{Point{3,Float64}}()        
    for (i,contourSet) in enumerate(contourSets)
        for V in contourSet     
            if !isempty(V)   
                append!(V_plot3D, V)
                push!(V_plot3D, V_nan)
                if i == sliceIndex
                    append!(V_plot3D_ax1, V)
                    push!(V_plot3D_ax1, V_nan)
                end        
            end
        end    
    end    
    # if !isempty(V_plot3D)   
        h1[1] = V_plot3D
        h2[1] = V_plot3D_ax1 
    # end   
end

function nearestContourIndex(xy_pos, contourSet)    
    p = Point{3,Float64}(xy_pos[1], xy_pos[2], dicomData[sliceIndex].SliceLocation)
    d = Inf 
    indexClose = 0                
    for (i,P) in enumerate(contourSet)
        dNow = minimum(norm.(P.-p))
        if dNow<d
            indexClose = i # Switch to current nearest
            d = dNow
        end
    end
    return indexClose
end


function load_contours(file_name)    
    doc = read(file_name, Node)
    contour_segmentation_node = doc.children[2]

    info_node = contour_segmentation_node.children[1]
    numSlices = parse(Int,info_node.children[1].children[1].value)
    voxelSize = parse.(Float64,split(info_node.children[2].children[1].value,","))
    contours_node = contour_segmentation_node.children[2]
    contourSets_accepted = [ Vector{Vector{Point{3,Float64}}}()  for _ in 1:numSlices]  
    for slice_node in children(contours_node)
        iSlice = parse(Int,slice_node.attributes["id"])
        for contour_node in children(slice_node)
            V = Vector{Point{3,Float64}}()        
            for point_node in children(contour_node)                
                s = split(point_node.children[1].value,",")
                p = Point{3,Float64}(parse(Float64,s[1]), parse(Float64,s[2]), parse(Float64,s[3]))            
                push!(V,p)
            end
            push!(contourSets_accepted[iSlice],V)
        end
    end
    return contourSets_accepted, voxelSize
end


function save_contours(save_path, contourSets_accepted, voxelSize)
    numSlices = length(contourSets_accepted)
    io = open(save_path, "w")
    write(io, """<?xml version="1.0" encoding="UTF-8"?> \n""")    
    write(io, "<contour_segmentation> \n")
    write(io, "<info> \n")
    write(io, "    <number_of_slices>$numSlices</number_of_slices>\n") 
    write(io, "    <voxel_size>" * @sprintf("%.16e, %.16e, %.16e", voxelSize[1], voxelSize[2], voxelSize[3]) * "</voxel_size>\n")
    write(io, "</info> \n")
    write(io, "<contours> \n")
    for (iSlice,contourSet) in enumerate(contourSets_accepted)
        if !isempty(contourSet)
            write(io, """    <slice id="$iSlice"> \n""")            
            for (iContour,V) in enumerate(contourSets_accepted[iSlice])            
                write(io, """        <contour id="$iContour"> \n""")
                for (iPoint,v) in enumerate(V)
                    write(io, """            <point id="$iPoint">""" * join([@sprintf("%.16e", x) for x ∈ v],',') * "</point> \n")            
                end            
                write(io, "        </contour> \n")
            end
            write(io, "    </slice> \n")            
        end
    end    
    write(io, "</contours>")
    write(io, "</contour_segmentation> \n")
    close(io)
end

function get_contour_extrema(contourSets)
    xMin =  Inf
    xMax = -Inf
    yMin =  Inf
    yMax = -Inf
    zMin =  Inf
    zMax = -Inf
    for contourSet in contourSets
        for P in contourSet
            for p in P
                xMin = min(xMin,p[1])
                xMax = max(xMax,p[1])
                yMin = min(yMin,p[2])
                yMax = max(yMax,p[2])
                zMin = min(zMin,p[3])
                zMax = max(zMax,p[3])            
            end
        end
    end
    return xMin, xMax, yMin, yMax, zMin, zMax
end

function contour2bool(contourSets, voxelSize)
    xMin, xMax, yMin, yMax, zMin, zMax = get_contour_extrema(contourSets)

    # Create image coordinate grid
    nAdd = 2.0
    xRange = xMin-nAdd*voxelSize[1]:voxelSize[1]:xMax+nAdd*voxelSize[1]
    yRange = yMin-nAdd*voxelSize[2]:voxelSize[2]:yMax+nAdd*voxelSize[2]
    zRange = zMin:voxelSize[3]:zMax

    indicesSlices = findall(.!isempty.(contourSets))
    B = fill(false, (length(xRange), length(yRange), length(indicesSlices)))
    for (k,iSlice) in enumerate(indicesSlices)
        for P in contourSets[iSlice]                
            for (i,x) in enumerate(xRange)
                for (j,y) in enumerate(yRange)                
                    s = inpolygon(Point{3,Float64}(x, y, P[1][3]),P)
                    if s >= 0 # In or on 
                        B[i,j,k] = !B[i,j,k] # Flip bool (this ensures complex geometry e.g. featuring holes etc, are treated properly)
                    end
                end            
            end
        end
    end
    return B, xRange, yRange, zRange
end

function boolBoundary(B)
    # Cross shaped kernel
    kernelCart = [CartesianIndex(-1,  0,  0), # Previous row
                  CartesianIndex( 1,  0,  0), # Next row
                  CartesianIndex( 0, -1,  0), # Previous column
                  CartesianIndex( 0,  1,  0), # Next column
                  CartesianIndex( 0,  0, -1), # Previous slice
                  CartesianIndex( 0,  0,  1)] # Next slice

    Bb = zeros(Bool, size(B)) # False array initialisation
    for ijk in findall(B) # For all true element in B            
        for ijk_shift in kernelCart # For each Cartesian index shift in the kernel 
            ijk_check = ijk + ijk_shift # The proposed CartesianIndex to check
            if checkbounds(Bool, B, ijk_check) && B[ijk_check] == false            
                Bb[ijk] = true # Set true if we found a false in the kernel                                               
                break # Stop if at least one was found
            end                        
        end
    end
    return Bb
end

function boolGrow!(B; n=1)
    if any(B) 
        # Cross shaped kernel
        kernelCart = [CartesianIndex(-1,  0,  0), # Previous row
                      CartesianIndex( 1,  0,  0), # Next row
                      CartesianIndex( 0, -1,  0), # Previous column
                      CartesianIndex( 0,  1,  0), # Next column
                      CartesianIndex( 0,  0, -1), # Previous slice
                      CartesianIndex( 0,  0,  1)] # Next slice    
        for _ in 1:n                  
            for ijk in findall(B) # For all true elements in B            
                for ijk_shift in kernelCart # For each Cartesian index shift in the kernel 
                    ijk_check = ijk + ijk_shift # The proposed CartesianIndex to check
                    if checkbounds(Bool, B, ijk_check) && B[ijk_check] == false # If valid and currently false
                        B[ijk_check] = true # Switch to true
                    end                        
                end
            end
        end
        return B
    else # Non-true so just return B
        return B
    end
    
end

function boolGrow(B; n=1)
    return boolGrow!(deepcopy(B); n=n)
end

function boolShrink!(B; n=1)
    if any(B) 
        # Cross shaped kernel
        kernelCart = [CartesianIndex(-1,  0,  0), # Previous row
                      CartesianIndex( 1,  0,  0), # Next row
                      CartesianIndex( 0, -1,  0), # Previous column
                      CartesianIndex( 0,  1,  0), # Next column
                      CartesianIndex( 0,  0, -1), # Previous slice
                      CartesianIndex( 0,  0,  1)] # Next slice  

        for _ in 1:n                  
            for ijk in findall(.!B) # For all false elements in B            
                for ijk_shift in kernelCart # For each Cartesian index shift in the kernel 
                    ijk_check = ijk + ijk_shift # The proposed CartesianIndex to check
                    if checkbounds(Bool, B, ijk_check) && B[ijk_check] # If valid and currently true                
                        B[ijk_check] = false # Switch to false
                    end                        
                end
            end
        end
        return B
    else # Non-true so just return B
        return B
    end
end

function boolShrink(B; n=1)
    return boolShrink!(deepcopy(B); n=n)
end

function contour2levelset(contourSets, voxelSize)
    B, xRange, yRange, zRange = contour2bool(contourSets, voxelSize)
    B2 = boolBoundary(B) # Boundary boolean
    boolGrow!(B2; n=1)
    M = fill(Float64(1),size(B))
    M[B] .= Float64(-1)
    indicesSlices = findall(.!isempty.(contourSets))
    for ijk in findall(B2)    
        d = Inf    
        for P in contourSets[indicesSlices[ijk[3]]]        
            for p in P
                dn = sqrt((xRange[ijk[1]]-p[1])^2 + (yRange[ijk[2]]-p[2])^2)
                d = min(d,dn)
            end                
        end
        if B[ijk]
            d *= -1.0                    
        end    
        M[ijk] =  d          
    end
    return M, xRange, yRange, zRange
end

################################################################################

function contoursegment(dcmFolder)    
    cursor_sample = GLFW.CreateStandardCursor(GLFW.CROSSHAIR_CURSOR)
    cursor_cut    = GLFW.CreateStandardCursor(GLFW.ARROW_CURSOR)
    cursor_draw   = GLFW.CreateStandardCursor(GLFW.HAND_CURSOR)
    cursor_accept = GLFW.CreateStandardCursor(GLFW.HAND_CURSOR)
    cursor_smooth = GLFW.CreateStandardCursor(GLFW.HAND_CURSOR)
    cursor_demote = GLFW.CreateStandardCursor(GLFW.HAND_CURSOR)
    cursor_delete = GLFW.CreateStandardCursor(GLFW.HAND_CURSOR)
    sliderWidth = 20 # Width for the sliders in the slidergrid

    global dicomData = dicomdir2dicomvec(dcmFolder)
    sliceKeySet = sort(collect(keys(dicomData)))
    startSlice = sliceKeySet[1] # ceil(Int,length(dicomData)/2)
    firstSlice = sliceKeySet[1] 
    lastSlice = sliceKeySet[end] 
    numSlices =  length(dicomData)
    dcm = dicomData[startSlice]
    raw_image = dcm.PixelData
    siz = (size(raw_image,1),size(raw_image,2),numSlices)
    if !isnothing(dcm.SpacingBetweenSlices)
        sliceStepSize = dcm.SpacingBetweenSlices
    elseif !isnothing(dcm.SliceThickness)
        sliceStepSize = dcm.SliceThickness
    elseif length(dicomData)>1
        sliceStepSize = dicomData[sliceKeySet[2]].SliceLocation-dicomData[sliceKeySet[1]].SliceLocation
    else
        sliceStepSize = 1.0
    end

    maxVal, minVal = getvalrange(dicomData) 

    global voxelSize = (dcm.PixelSpacing[1], dcm.PixelSpacing[2], sliceStepSize)
    # sliceSpacing = dcm.SpacingBetweenSlices
    imagePositionPatient = dcm.ImagePositionPatient
    imageOrientationPatient = dcm.ImageOrientationPatient


    initialSliceIndices = [ceil.(Int,siz[1]./2), ceil.(Int,siz[2]./2), startSlice]
    global sliceIndex = initialSliceIndices[3]
    global contourSet_draft = Vector{Vector{Point{3,Float64}}}()
    global contourSets_accepted = [ Vector{Vector{Point{3,Float64}}}()  for _ in 1:numSlices]  
    global V_cut = [Point{2,Float64}(NaN,NaN), Point{2,Float64}(NaN,NaN), Point{2,Float64}(NaN,NaN), Point{2,Float64}(NaN,NaN)]
    global V_draw = Vector{Point{2,Float64}}()
    global indexNearest = 0
    global xy_click = Vector{Float64}()
    global contourLevel = zero(eltype(dicomData[startSlice].PixelData))
    global λ = 1.0
    global colorbarLimit_upper = maxVal
    global colorbarLimit_lower = minVal
    global save_path = joinpath(imagodir(),"temp.xml")
    global x = 0.0:voxelSize[1]:voxelSize[1]*siz[1]
    global y = 0.0:voxelSize[2]:voxelSize[2]*siz[2]
    global z = dicomData[firstSlice].SliceLocation:voxelSize[3]:(dicomData[firstSlice].SliceLocation + (voxelSize[3]*siz[3]))

    A = getslice(dicomData,initialSliceIndices[1],1)
    B = getslice(dicomData,initialSliceIndices[2],2)

    ## Visualization------------------------------------------
    cmap = :grays
    figSize = (1800, 1200)
    fig = Figure(size=figSize)

    Label(fig[1, 3][1, 1], "Smooth raw")
    global toggle = Toggle(fig[1, 3][1,2], active = false)

    Label(fig[1, 3][2, 1], "λ: ")
    h_textBox_smooth = Textbox(fig[1, 3][2, 2], placeholder = "$λ")

    h_textBox_save = Textbox(fig[2, :], placeholder = save_path)
    # h_textBox_save.displayed_string = save_path
    # h_textBox_save.stored_string = save_path
    h_button_save = Button(fig[1, 3][3, 2], label = "Save")

    h_button_load = Button(fig[1, 3][4, 2], label = "Load")

    Label(fig[1, 3][5, 1], "cMax: ")
    h_textBox_cMax = Textbox(fig[1, 3][5, 2], placeholder = "$maxVal")

    Label(fig[1, 3][6, 1], "cMin: ")
    h_textBox_cMin = Textbox(fig[1, 3][6, 2], placeholder = "$minVal")

    ax1 = LScene(fig[1,1]); 
    # ax1 = AxisGeom(fig[1, 1], limits=(0,siz[1]*voxelSize[1],0,siz[2]*voxelSize[2],0,siz[3]*voxelSize[3]))


    hs1 = heatmap!(ax1,z,y,A, colormap = cmap, colorrange=(colorbarLimit_lower, colorbarLimit_upper))
    α = -0.5*pi
    a = (0.0,1.0,0.0) # x-y-z
    k = Makie.Quaternion(sin(α/2)*a[1], sin(α/2)*a[2], sin(α/2)*a[3], cos(α/2))
    Makie.rotate!(hs1, k)
    Makie.translate!(hs1, initialSliceIndices[1]*voxelSize[1], 0.0, 0.0)

    hs2 = heatmap!(ax1,x,z,B, colormap = cmap, colorrange=(colorbarLimit_lower, colorbarLimit_upper))
    α = 0.5*pi
    a = (1.0,0.0,0.0) # x-y-z
    k = Makie.Quaternion(sin(α/2)*a[1], sin(α/2)*a[2], sin(α/2)*a[3], cos(α/2))
    Makie.rotate!(hs2, k)
    Makie.translate!(hs2, 0.0, initialSliceIndices[2]*voxelSize[2], 0.0)

    hs3 = heatmap!(ax1,x,y,raw_image, colormap = cmap, colorrange=(colorbarLimit_lower, colorbarLimit_upper))
    Makie.translate!(hs3, 0.0, 0.0, dicomData[startSlice].SliceLocation)


    # ax2 = LScene(fig[1,1]); 
    ax2 = Axis(fig[1, 2], aspect = DataAspect(), xlabel = "X", ylabel = "Y", title="Sample contour", limits=(0,siz[1]*voxelSize[1],0,siz[2]*voxelSize[2]),xrectzoom=false, yrectzoom=false)
    hs4 = heatmap!(ax2,x,y,raw_image, colormap = cmap, colorrange=(colorbarLimit_lower, colorbarLimit_upper))

    currentMode = Observable(:sample)

    # Initialise plots (not visible for now) 
    V_empty = Vector{Point{3,Float64}}()   

    h1_draft    = lines!(ax1, V_empty, color=:yellow, linewidth=6, linestyle=(:dot,:dense))
    h1_accepted = lines!(ax1, V_empty, color=:green, linewidth=8)

    h2_draft_bg = lines!(ax2, V_empty, color = :white, linewidth=4)
    h2_draft    = lines!(ax2, V_empty, color = ones(Int,length(V_empty)), linewidth=2, colormap=:viridis) #, linestyle=(:dot,:dense)
    h2_accepted = lines!(ax2, V_empty, color=:green, linewidth=4)

    hCut   = lines!(ax2, [Point{2,Float64}(NaN,NaN)], color = :red, linewidth=2, visible=false) #, linestyle=(:dot,:dense)
    hDraw   = lines!(ax2, [Point{2,Float64}(NaN,NaN)], color = :red, linewidth=2, visible=false) #, linestyle=(:dot,:dense)

    # update_plots((hl3), val; use_color=false, make_2D=false)

    on(h_textBox_smooth.stored_string) do s
        global λ = parse(Float64, s)
    end

    on(h_textBox_save.stored_string) do s
        global save_path = s#parse(String, s)
    end

    on(h_textBox_cMax.stored_string) do c
        cVal = parse(Float64, c)
        if cVal > colorbarLimit_lower 
            global colorbarLimit_upper = cVal
            hs1.colorrange = (colorbarLimit_lower, colorbarLimit_upper)
            hs2.colorrange = (colorbarLimit_lower, colorbarLimit_upper)
            hs3.colorrange = (colorbarLimit_lower, colorbarLimit_upper)
            hs4.colorrange = (colorbarLimit_lower, colorbarLimit_upper)
        else
            h_textBox_cMax.stored_string = string(colorbarLimit_upper)
            h_textBox_cMax.displayed_string = string(colorbarLimit_upper)
        end     
    end

    on(h_textBox_cMin.stored_string) do c
        cVal = parse(Float64, c)
        if cVal < colorbarLimit_upper
            global colorbarLimit_lower = cVal
            hs1.colorrange = (colorbarLimit_lower, colorbarLimit_upper)
            hs2.colorrange = (colorbarLimit_lower, colorbarLimit_upper)
            hs3.colorrange = (colorbarLimit_lower, colorbarLimit_upper)
            hs4.colorrange = (colorbarLimit_lower, colorbarLimit_upper)
        else
            h_textBox_cMin.stored_string = string(colorbarLimit_lower)
            h_textBox_cMin.displayed_string = string(colorbarLimit_lower)
        end
    end


    on(h_button_save.clicks) do n
        save_contours(save_path, contourSets_accepted, voxelSize)
    end

    on(h_button_load.clicks) do n
        contourSets_imported, voxelSize = load_contours(save_path)  
        global contourSets_accepted = contourSets_imported
        update_accepted_plot(h1_accepted, h2_accepted, contourSets_accepted)
    end

    register_interaction!(ax2, :my_interaction) do event::MouseEvent, axis        
        if currentMode[] == :sample
            if event.type === MouseEventTypes.leftdown
                raw_image = hs4[3][]
                V_contour_raw = contour_at_click(event.data, voxelSize, raw_image)      
                global contourSet_draft = [V_contour_raw]
                global xy_click = [V_contour_raw[1][1], V_contour_raw[1][2]]

                update_draft_plot(h1_draft, contourSet_draft; use_color=false)
                update_draft_plot(h2_draft_bg, contourSet_draft; use_color=false)
                update_draft_plot(h2_draft, contourSet_draft; use_color=true)
            end
        elseif currentMode[] == :cut        
            if !isempty(contourSet_draft)
                if event.type === MouseEventTypes.leftdown
                    hCut.visible = true
                    xy_pos = event.data     
                    if isnan(V_cut[1]) # First click
                        global V_cut[1] = Point{2,Float64}(xy_pos[1], xy_pos[2])                    
                    else # Second click
                        global V_cut[2] = Point{2,Float64}(xy_pos[1], xy_pos[2])
                        if V_cut[2][1]<V_cut[1][1] && V_cut[2][2]<V_cut[2][1]
                            reverse!(V_cut)
                        elseif V_cut[2][1]>V_cut[1][1] && V_cut[2][2]<V_cut[2][1]
                            reverse!(V_cut)
                        end
                        V_poly = [V_cut[1], 
                                Point{2,Float64}(V_cut[2][1],V_cut[1][2]),
                                V_cut[2],
                                Point{2,Float64}(V_cut[1][1],V_cut[2][2])]
                            
                        V_groups_all = Vector{Vector{Point{3,Float64}}}()    
                        for (i,V_contour_now) in enumerate(contourSet_draft)
                            V_group = Vector{Vector{Point{3,Float64}}}()    
                            V_part = Vector{Point{3,Float64}}()
                            firstTrue = false
                            lastTrue = false
                            cutMade = false
                            for (j,p) in enumerate(V_contour_now)
                                f = inpolygon(Point{2,Float64}(p[1], p[2]), V_poly)
                                if f!=-1 # Inside (or on)
                                    cutMade = true
                                    if !isempty(V_part)
                                        push!(V_group,V_part) # Add group
                                    end
                                    V_part = Vector{Point{3,Float64}}() # Empty group
                                else # Outside 
                                    push!(V_part,p) # Add point to group
                                    if j == 1
                                        firstTrue = true
                                    elseif j == length(V_contour_now)
                                        lastTrue = true
                                        push!(V_group,V_part) # Add group
                                    end
                                end
                            end
                            if cutMade && firstTrue && lastTrue                            
                                V_first = V_group[1]
                                V_last = V_group[end]
                                if norm(V_last[end]-V_first[1]) <= sqrt(sum(voxelSize[1:2].^2))
                                    deleteat!(V_group,length(V_group))
                                    prepend!(V_group[1],V_last)
                                end
                            end  
                            append!(V_groups_all,V_group)
                        end
                        contourSet_draft = V_groups_all
                        update_draft_plot(h1_draft, contourSet_draft; use_color=false)
                        update_draft_plot(h2_draft_bg, contourSet_draft; use_color=false)
                        update_draft_plot(h2_draft, contourSet_draft; use_color=true)            
                        global V_cut = [Point{2,Float64}(NaN,NaN), Point{2,Float64}(NaN,NaN)]                      
                        hCut[1] = V_cut
                        hCut.visible = false           
                    end                     
                elseif event.type === MouseEventTypes.over
                    if !isnan(V_cut[1])
                        xy_pos = event.data
                        p = Point{2,Float64}(xy_pos[1],xy_pos[2])

                        V_poly = [V_cut[1], 
                                Point{2,Float64}(p[1],V_cut[1][2]),
                                p,
                                Point{2,Float64}(V_cut[1][1],p[2]),
                                V_cut[1]]
                        hCut[1] = V_poly                    
                    end
                elseif event.type === MouseEventTypes.rightclick
                    currentMode[] = :sample  
                end
            else
                currentMode[] = :sample  
            end
        elseif currentMode[] == :draw                  
            if event.type === MouseEventTypes.leftdown
                xy_pos = event.data
                p = Point{2,Float64}(xy_pos[1],xy_pos[2])
                push!(V_draw,p)
                hDraw[1] = V_draw
                hDraw.visible = true
            elseif event.type === MouseEventTypes.rightclick
                push!(contourSet_draft,[Point{3,Float64}(v[1], v[2], dicomData[sliceIndex].SliceLocation) for v in V_draw])
                global V_draw = Vector{Point{2,Float64}}()
                hDraw.visible = false            
                hDraw[1] = [Point{2,Float64}(NaN,NaN)]
                update_draft_plot(h1_draft, contourSet_draft; use_color=false)
                update_draft_plot(h2_draft_bg, contourSet_draft; use_color=false)
                update_draft_plot(h2_draft, contourSet_draft; use_color=true)   
                currentMode[] = :sample  
            end
        elseif currentMode[] == :accept            
            if !isempty(contourSet_draft)
                if event.type === MouseEventTypes.leftdown                   
                    if indexNearest == 0 # First time                                      
                        global indexNearest = nearestContourIndex(event.data, contourSet_draft)                    
                        if isempty(contourSets_accepted[sliceIndex])
                            push!(contourSets_accepted[sliceIndex],contourSet_draft[indexNearest])          
                            update_accepted_plot(h1_accepted, h2_accepted, contourSets_accepted)

                            # Remove accepted from draft set and update plot 
                            deleteat!(contourSet_draft,indexNearest)
                            update_draft_plot(h1_draft, contourSet_draft; use_color=false)
                            update_draft_plot(h2_draft_bg, contourSet_draft; use_color=false)
                            update_draft_plot(h2_draft, contourSet_draft; use_color=true) 

                            global indexNearest = 0 # Reset to zero
                        else
                            ax2.title = "Accept draft contour, -> l.click to merge with nearest OR r.click to add as new."
                        end                    
                    else                    
                        indexNearestAccepted = nearestContourIndex(event.data, contourSets_accepted[sliceIndex])
                        pStart_contour = contourSet_draft[indexNearest][1]
                        pEnd_contour = contourSet_draft[indexNearest][end]
                        pStart_current = contourSets_accepted[sliceIndex][indexNearestAccepted][1]
                        pEnd_current = contourSets_accepted[sliceIndex][indexNearestAccepted][end]

                        d_start_start = norm(pStart_current - pStart_contour)
                        d_start_end = norm(pStart_current - pEnd_contour)
                        d_end_start = norm(pEnd_current - pStart_contour)
                        d_end_end = norm(pEnd_current - pEnd_contour)
                        _,indMin = findmin([d_start_start, d_start_end, d_end_start, d_end_end])
                        
                        if indMin == 1 # Contour start is closest to current start 
                            prepend!(contourSets_accepted[sliceIndex][indexNearestAccepted],reverse(contourSet_draft[indexNearest]))
                        elseif indMin == 2 # Contour end is closest to current start 
                            prepend!(contourSets_accepted[sliceIndex][indexNearestAccepted],contourSet_draft[indexNearest])
                        elseif indMin == 3 # Contour start is closest to current end
                            append!(contourSets_accepted[sliceIndex][indexNearestAccepted],contourSet_draft[indexNearest])
                        elseif indMin == 4 # Contour end is closest to current end                
                            append!(contourSets_accepted[sliceIndex][indexNearestAccepted],reverse(contourSet_draft[indexNearest]))
                        end
                        update_accepted_plot(h1_accepted, h2_accepted, contourSets_accepted)

                        # Remove accepted from draft set and update plot 
                        deleteat!(contourSet_draft,indexNearest)
                        update_draft_plot(h1_draft, contourSet_draft; use_color=false)
                        update_draft_plot(h2_draft_bg, contourSet_draft; use_color=false)
                        update_draft_plot(h2_draft, contourSet_draft; use_color=true) 
                        global indexNearest = 0 # Reset to zero
                    end
                elseif event.type === MouseEventTypes.rightclick
                    if indexNearest == 0 # First time        
                        currentMode[] = :sample         
                    else
                        push!(contourSets_accepted[sliceIndex],contourSet_draft[indexNearest])
                        update_accepted_plot(h1_accepted, h2_accepted, contourSets_accepted)

                        # Remove accepted from draft set and update plot 
                        deleteat!(contourSet_draft,indexNearest)
                        update_draft_plot(h1_draft, contourSet_draft; use_color=false)
                        update_draft_plot(h2_draft_bg, contourSet_draft; use_color=false)
                        update_draft_plot(h2_draft, contourSet_draft; use_color=true) 

                        global indexNearest = 0 # Reset to zero
                    end
                end
            else
                currentMode[] = :sample  
            end
        elseif currentMode[] == :smooth
            if event.type === MouseEventTypes.leftdown  
                if !isempty(contourSets_accepted[sliceIndex])                
                    indexNearestAccepted = nearestContourIndex(event.data, contourSets_accepted[sliceIndex])
                    contourSets_accepted[sliceIndex][indexNearestAccepted] = smooth_contour(contourSets_accepted[sliceIndex][indexNearestAccepted]; λ = λ)
                    update_accepted_plot(h1_accepted, h2_accepted, contourSets_accepted)
                end
            elseif event.type === MouseEventTypes.rightclick
                currentMode[] = :sample  
            end
        elseif currentMode[] == :demote            
            if event.type === MouseEventTypes.leftdown                  
                if !isempty(contourSets_accepted[sliceIndex])
                    contourSet = contourSets_accepted[sliceIndex]
                    indexNearestAccepted = nearestContourIndex(event.data, contourSet)
                    push!(contourSet_draft, contourSet[indexNearestAccepted])
                    deleteat!(contourSet,indexNearestAccepted)
                    contourSets_accepted[sliceIndex] = contourSet
                    update_accepted_plot(h1_accepted, h2_accepted, contourSets_accepted)
                    update_draft_plot(h1_draft, contourSet_draft; use_color=false)
                    update_draft_plot(h2_draft_bg, contourSet_draft; use_color=false)
                    update_draft_plot(h2_draft, contourSet_draft; use_color=true)                 
                end
            elseif event.type === MouseEventTypes.rightclick
                currentMode[] = :sample  
            end
        elseif currentMode[] == :delete
            if !isempty(contourSet_draft)
                if event.type === MouseEventTypes.leftdown  
                    indexDelete = nearestContourIndex(event.data, contourSet_draft)
                    deleteat!(contourSet_draft,indexDelete)
                    update_draft_plot(h1_draft, contourSet_draft; use_color=false)
                    update_draft_plot(h2_draft_bg, contourSet_draft; use_color=false)
                    update_draft_plot(h2_draft, contourSet_draft; use_color=true) 
                elseif event.type === MouseEventTypes.rightclick
                    currentMode[] = :sample  
                end
            else
                currentMode[] = :sample  
            end
        end
    end

    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press || event.action == Keyboard.repeat        
            if event.key == Keyboard.c                
                currentMode[] = :cut
            elseif event.key == Keyboard.d                
                currentMode[] = :draw
            elseif event.key == Keyboard.a                
                currentMode[] = :accept
            elseif event.key == Keyboard.s                
                currentMode[] = :smooth
            elseif event.key == Keyboard.backspace                
                currentMode[] = :demote            
            elseif event.key == Keyboard.delete                
                currentMode[] = :delete            
            else
                hDraw.visible = false            
                hDraw[1] = [Point{2,Float64}(NaN,NaN)]                
                currentMode[] = :sample             
            end        
        end
    end

    on(currentMode) do val
        if val == :sample        
            ax2.title = "Sample draft contour"
            GLFW.SetCursor(window, cursor_sample)
        elseif val == :cut
            ax2.title = "Cut draft contour, l.click first and second corner"
            GLFW.SetCursor(window, cursor_cut)
        elseif val == :draw
            ax2.title = "Draw draft contour, l.click to add points"        
            GLFW.SetCursor(window, cursor_draw)
        elseif val == :accept
            ax2.title = "Accept draft contour, l.click to accept nearest."
            if length(contourSet_draft)==1 && isempty(contourSets_accepted[sliceIndex])
                # Only one draft contour so no need to select, just accept this one
                push!(contourSets_accepted[sliceIndex],contourSet_draft[1])
                update_accepted_plot(h1_accepted, h2_accepted, contourSets_accepted)

                # Empty draft contours now and update plots
                global contourSet_draft = Vector{Vector{Point{3,Float64}}}()
                update_draft_plot(h1_draft, contourSet_draft; use_color=false)
                update_draft_plot(h2_draft_bg, contourSet_draft; use_color=false)
                update_draft_plot(h2_draft, contourSet_draft; use_color=true) 
            end        
            GLFW.SetCursor(window, cursor_accept)
        elseif val == :smooth
            ax2.title = "Smooth accepted contour"
            if length(contourSets_accepted[sliceIndex]) == 1# Just one so just smooth that one
                contourSets_accepted[sliceIndex][1] = smooth_contour(contourSets_accepted[sliceIndex][1]; λ = λ)
                update_accepted_plot(h1_accepted, h2_accepted, contourSets_accepted)
            end        
            GLFW.SetCursor(window, cursor_smooth)
        elseif val == :demote
            ax2.title = "Demote accepted contour"
            if length(contourSets_accepted[sliceIndex]) == 1# Just one so just smooth that one
                push!(contourSet_draft, contourSets_accepted[sliceIndex][1])
                contourSets_accepted[sliceIndex]=Vector{Vector{Point{3,Float64}}}()
                update_accepted_plot(h1_accepted, h2_accepted, contourSets_accepted)
                update_draft_plot(h1_draft, contourSet_draft; use_color=false)
                update_draft_plot(h2_draft_bg, contourSet_draft; use_color=false)
                update_draft_plot(h2_draft, contourSet_draft; use_color=true) 
            end        
            GLFW.SetCursor(window, cursor_demote)
        elseif val == :delete
            ax2.title = "Delete draft contour"
            if length(contourSet_draft) == 1# Just one so just smooth that one
                # Empty draft contours now and update plots
                global contourSet_draft = Vector{Vector{Point{3,Float64}}}()
                update_draft_plot(h1_draft, contourSet_draft; use_color=false)
                update_draft_plot(h2_draft_bg, contourSet_draft; use_color=false)
                update_draft_plot(h2_draft, contourSet_draft; use_color=true) 
            end        
            GLFW.SetCursor(window, cursor_delete)
        end
    end

    sg = SliderGrid(
        fig[3, 1:3],
        (label = "X-slice", range = 1:1:siz[1], startvalue = ceil(Int,siz[1]/2)),
        (label = "Y-slice", range = 1:1:siz[2], startvalue = ceil(Int,siz[2]/2)),
        (label = "Z-slice", range = sliceKeySet[1]:1:sliceKeySet[end], startvalue = ceil(Int,siz[3]/2)),
        tellheight = true, valign=:bottom)

    sg.sliders[1].linewidth = sliderWidth
    sg.sliders[2].linewidth = sliderWidth
    sg.sliders[3].linewidth = sliderWidth

    hSlider1 = sg.sliders[1] #Slider(fig[3, :], range = 1:1:siz[1], startvalue = ceil(Int,siz[1]/2), linewidth=30)
    on(hSlider1.value) do stepIndex 
        hs1[3] = getslice(dicomData,stepIndex,1)
        Makie.translate!(hs1, stepIndex*voxelSize[1], 0.0, 0.0)
    end

    hSlider2 = sg.sliders[2] # Slider(fig[4, :], range = 1:1:siz[2], startvalue = ceil(Int,siz[2]/2), linewidth=30)
    on(hSlider2.value) do stepIndex 
        hs2[3] = getslice(dicomData,stepIndex,2)
        Makie.translate!(hs2, 0.0, stepIndex*voxelSize[2], 0.0)
    end

    hSlider3 = sg.sliders[3]#Slider(fig[5, :], range = sliceKeySet, startvalue = startSlice, linewidth=30)
    on(hSlider3.value) do stepIndex 
        global sliceIndex = stepIndex

        global contourSet_draft = Vector{Vector{Point{3,Float64}}}()
        
        update_accepted_plot(h1_accepted, h2_accepted, contourSets_accepted)
        
        hs3[3] = hs4[3] = dicomData[sliceIndex].PixelData       
        Makie.translate!(hs3, 0.0, 0.0, dicomData[sliceIndex].SliceLocation)    
        
        if !isempty(xy_click)
            V_contour_raw = contour_at_click(xy_click, voxelSize, dicomData[sliceIndex].PixelData, c=contourLevel)
            if !isempty(V_contour_raw)
                global contourSet_draft = [V_contour_raw]
                update_draft_plot(h1_draft, contourSet_draft; use_color=false)
                update_draft_plot(h2_draft_bg, contourSet_draft; use_color=false)
                update_draft_plot(h2_draft, contourSet_draft; use_color=true) 
            end  
        end  
    end

    slidercontrol(hSlider3,fig) 

    Colorbar(fig[:,4],hs4)
    rowgap!(fig.layout,2)
    colgap!(fig.layout,2)
    colsize!(fig.layout, 1, Relative(0.45))
    colsize!(fig.layout, 2, Relative(0.45))
    rowsize!(fig.layout, 1, Relative(0.9))
    # rowsize!(fig.layout, 2, Relative(0.9))
    screen = display(fig)
    window = GLMakie.to_native(screen)
    GLFW.SetCursor(window, cursor_sample)

    return fig, ax1, ax2
end