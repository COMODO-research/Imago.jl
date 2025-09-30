"""
    imagodir()

# Description 

This function simply returns the string for the Imago path. This is helpful for instance to load items, such as meshes, from the `assets` folder. 
"""
function imagodir()
    pkgdir(@__MODULE__)
end

function dicomdir2dicomvec(dcmFolder)
    fileSet = readdir(dcmFolder)

    # Read in first 
    dcm = dcm_parse(dcmFolder*"/"*fileSet[1])
    dicomData = Dict{Int, typeof(dcm)}()
    # dicomData = Vector{typeof(dcm)}(undef,max_instanceNumber)
    for f in fileSet # for each file 
        dcm = dcm_parse(dcmFolder*"/"*f)
        dicomData[dcm.InstanceNumber] = dcm
    end
    return dicomData
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
    if dir ==1 
        A = Matrix{T}(undef,n,siz[dir])
    else
        A = Matrix{T}(undef,siz[dir],n)
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
    println(joinpath(imagodir(), "assets", "data", dataName))    
    dataDir = joinpath(imagodir(), "assets", "data", dataName)
    dataDir_unzip = joinpath(dataDir,"data")
    if isdir(dataDir_unzip)
        return dataDir_unzip
    # else
    #     mkdir(dataDir_unzip) # Create directory 
    #     fileName_zip = joinpath(dataDir,"data.zip")
    #     if isfile(fileName_zip)
    #         iobuffer = IOBuffer(readstring(fileName_zip))
    #         r = ZipFile.Reader(iobuffer)
    #     end        
    end
end
