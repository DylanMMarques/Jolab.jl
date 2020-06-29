module Jolab
	using LinearAlgebra
	using SpecialFunctions
	using Roots
	using Interpolations
	using StaticArrays
	using DSP
	using FFTW
	using HCubature
	import FunctionWrappers
	import FunctionWrappers: FunctionWrapper
	using Random

	export FieldAngularSpectrum, FieldAngularSpectrum_fromspace, FieldAngularSpectrum_fromspacefft, FieldAngularSpectrum_gaussian, intensity, changereferenceframe!, scalartovectorial!, vectorialtoscalar!, scalartovectorial, vectorialtoscalar, FieldAngularSpectrumSymmetric_gaussian, FieldAngularSpectrumSymmetric
	export FieldSpace, FieldSpace_gaussian, FieldSpace_uniform, FieldSpace_fromangspe, FieldSpace_fromangspefft, FieldSpaceSymmetric, FieldSpaceSymmetric_gaussian
	export angspeto3Dspace
	export PointDetector, signal
	export ReferenceFrame
	export lightinteraction, MultilayerStructure
	export SingleModeFibre, signal, FieldSpace_fromfibre, FieldAngularSpectrum_fromfibre, FieldSpaceSymmetric_fromfibre, FieldAngularSpectrumSymmetric_fromfibre
	export numberofmodes
	export adddims
	export Lens
	export RoughInterface, lightinteraction_recursive
	export PropagationCoefficientScalar
	export Mirror, rotatestructure, translatestructure
	export SpatialLightModulator

	export CircularStepIndexModes, CircularStepIndexFibre, findmodes!


	macro tol(T=1.)
    	return convert(typeof(T), 1E-15)
	end
	function tobedone()
		error("Sorry, function not implemented yet.")
	end

	abstract type AbstractCoefficient{T<:Real} end
	abstract type AbstractCoefficientCache{T<:Real} end
	abstract type AbstractOpticalComponent{T<:Real} end
	abstract type AbstractDetector{T<:Real} <: AbstractOpticalComponent{T} end

	abstract type AbstractFieldMonochromatic{T<:Real} end
	abstract type AbstractModes{T<:Real} end
	abstract type AbstractPropagationComponent{T<:Real} <: AbstractOpticalComponent{T} end
	abstract type AbstractFieldAngularSpectrum{T} <: AbstractFieldMonochromatic{T} end
	abstract type AbstractFieldSpace{T} <: AbstractFieldMonochromatic{T} end

	include("ReferenceFrame.jl")

	include("AbstractFieldMonochromatic.jl")
	include("AbstractCoefficient.jl")
	include("AbstractOpticalComponent.jl")
	include("AbstractDetector.jl")

	include("numericalFunctions.jl");
	include("JolabFunction.jl")
	include("RefractiveIndexDatabase.jl")
	include("AngularSpectrum.jl");
	include("AngularSpectrumSymmetric.jl");
	include("FieldSpace.jl")
	include("FieldSpaceSymmetric.jl")
	include("FieldModes.jl")
	include("PointDetector.jl")
	include("PropagationCoefficient.jl")
	include("ScatteringCoefficient.jl")
	include("MultilayerStructure.jl")
	include("RoughInterface.jl")
	include("Lens.jl")
	include("SingleModeFibre.jl")
	include("Mirror.jl")
	include("interpolation.jl")
	include("CircularStepIndexFibre.jl")
	include("SpatialLightModulator.jl")
	include("ScatteringMatrix.jl")
	include("RecursiveAlgorithms.jl")
	include("PlotsRecipes.jl")
end
