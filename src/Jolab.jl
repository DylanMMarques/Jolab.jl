module Jolab

using CoordinateTransformations, Rotations, StaticArrays, ArgCheck, StructArrays, FillArrays, Rotations, LinearAlgebra
using Bessels, Roots, HCubature

export light_interaction, intensity
export translate_referenceframe, rotate_referenceframe
export Forward, Backward

import StructArrays: component

abstract type AbstractField{T,D} end
abstract type AbstractMode{T} end
abstract type AbstractFieldMode{T,D} <: AbstractField{T,D} end
abstract type AbstractPlaneWave{T,D} <: AbstractFieldMode{T,D} end
abstract type AbstractMedium{T} end
abstract type AbstractOpticalElement{T} end

abstract type AbstractDirection end
struct Forward <: AbstractDirection end
struct Backward <: AbstractDirection end

const RealOrComplex{T} = Union{T, Complex{T}}

include("Medium.jl")
include("ReferenceFrame.jl")
include("PlaneWave.jl")
include("Beam.jl")
include("ScatteringMatrix.jl")

include("DielectricStack.jl")
# include("SingleModeFibre.jl")
include("Fibre.jl")

include("auxiliary_functions.jl")
end