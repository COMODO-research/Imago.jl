module Imago

using Comodo
using Comodo.GeometryBasics
using Comodo.BSplineKit
using Comodo.GLMakie
using Comodo.GLMakie: GLFW, Colors
using Comodo.Statistics
using Comodo.LinearAlgebra

using FileIO
using DICOM
import Contour
using Printf
using XML

# Include functions
include("functions.jl")

# Export imported modules for later possible use
export Comodo
export BSplineKit
export GLMakie
export GLFW
export Colors
export DICOM
export XML
export Printf
export Contour
export FileIO

# Export functions
export imagodir, dicomdir2dicomvec, getvalrange, getslice, getdemodata, get_image_file_set, contoursegment

end # module Imago
