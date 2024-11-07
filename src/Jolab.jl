module Jolab

using CoordinateTransformations, Rotations, StaticArrays, ArgCheck, StructArrays, FillArrays, Rotations, LinearAlgebra
using Bessels, HCubature

export light_interaction, intensity
export translate_referenceframe, rotate_referenceframe
export Forward, Backward
export ScatteringMatrix
export MonochromaticAngularSpectrum, MonochromaticSpatialBeam, MeshedPlaneWaveScalar, MonochromaticAngularSpectrum_gaussian, MonochromaticSpatialBeam_gaussian
export MonochromaticAngularSpectrumRadialSymmetric, MonochromaticAngularSpectrumRadialSymmetric_gaussian
export MonochromaticSpatialBeamRadialSymmetric_gaussian, MonochromaticSpatialBeamRadialSymmetric
export Propagation
export Fourier

import StructArrays: component

abstract type AbstractField{T,D} end
abstract type AbstractMode{T} end
abstract type AbstractFieldMode{T,D} <: AbstractField{T,D} end
abstract type AbstractPlaneWave{T,D} <: AbstractFieldMode{T,D} end
abstract type AbstractPointSource{T,D} <: AbstractFieldMode{T,D} end
abstract type AbstractMedium{T} end
abstract type AbstractOpticalElement{T} end

abstract type AbstractDirection end
struct Forward <: AbstractDirection end
struct Backward <: AbstractDirection end
struct Bothway <: AbstractDirection end

const RealOrComplex{T} = Union{T, Complex{T}}

include("errors.jl")
include("arrays.jl")
include("meshes.jl")
include("Medium.jl")
include("ReferenceFrame.jl")
include("Beam.jl")
include("MeshedBeam.jl")
include("ScatteringMatrix.jl")
include("Propagation.jl")

include("AbstractOpticalElement.jl")
include("DielectricStack.jl")
include("Lens.jl")
include("Fourier.jl")
include("Axicon.jl")
include("Fibre.jl")
include("SingleModeFibre.jl")

include("auxiliary_functions.jl")

function __init__()
end
end