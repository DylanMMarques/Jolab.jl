module Jolab

using CoordinateTransformations, Rotations, StaticArrays, ArgCheck, StructArrays, FillArrays, Rotations, LinearAlgebra

export light_interaction, intensity
export translate_referenceframe, rotate_referenceframe
export Forward, Backward

abstract type AbstractMode{T} end
abstract type AbstractFieldMode{T,D} end
abstract type AbstractPlaneWave{T,D} <: AbstractFieldMode{T,D} end
abstract type AbstractMedium{T} end
abstract type AbstractOpticalElement{T} end

abstract type AbstractDirection end
struct Forward <: AbstractDirection end
struct Backward <: AbstractDirection end

include("Medium.jl")
include("ReferenceFrame.jl")
include("PlaneWave.jl")
include("Beam.jl")
include("DielectricStack.jl")
include("SingleModeFibre.jl")

include("auxiliary_functions.jl")
end