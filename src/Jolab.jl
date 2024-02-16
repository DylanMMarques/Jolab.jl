module Jolab

using CoordinateTransformations, Rotations, StaticArrays, ArgCheck, StructArrays, FillArrays, Rotations, LinearAlgebra

export light_interaction, intensity
export translate_referenceframe, rotate_referenceframe

abstract type AbstractMode{T} end
abstract type AbstractPlaneWave{T} <: AbstractMode{T} end
abstract type AbstractMedium{T} end
abstract type AbstractOpticalElement{T} end

include("Medium.jl")
include("ReferenceFrame.jl")
include("PlaneWave.jl")
include("Beam.jl")
include("DielectricStack.jl")
include("SingleModeFibre.jl")
end