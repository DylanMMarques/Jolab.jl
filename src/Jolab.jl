module Jolab
	using LinearAlgebra, SparseArrays, ProgressMeter
	using SpecialFunctions
	using Roots
	using Interpolations
	using StaticArrays
	using DSP
	using FFTW
	using HCubature
	using Random, Statistics

	export FieldAngularSpectrumScalar, FieldAngularSpectrumScalar_gaussian, FieldAngularSpectrumScalar_uniform
	export intensity, changereferenceframe!, changereferenceframe
	export FieldSpaceScalar, FieldSpaceScalar_gaussian, FieldSpaceScalar_uniform
	export PointDetector, signal
	export ReferenceFrame
	export lightinteraction, MultilayerStructure
	export SingleModeFibre, signal, FieldSpace_fromfibre, FieldAngularSpectrum_fromfibre
	export numberofmodes
	export Lens
	export RoughMultilayerStructure, lightinteraction_recursive, lightinteraction_recursivegridded
	export PropagationCoefficientScalar
	export Mirror, rotatestructure, translatestructure
	export SpatialLightModulator, SpatialLightModulator_slm, SpatialLightModulator_mask, SpatialLightModulator_reflectivemask
	export FourierTransform, PhantomLayer
	export Axicon, AxiconFourier

	export CircularStepIndexModes, CircularStepIndexFibre, findmodes!

	macro tol(T=1.)
    	return convert(typeof(T), 1E-15)
	end
	function tobedone()
		error("Sorry, function not implemented yet.")
	end

	abstract type AbstractField{T<:Real, D} end
	abstract type AbstractFieldMonochromatic{T<:Real, D} end
	abstract type AbstractFieldAngularSpectrum{T,D} <: AbstractFieldMonochromatic{T,D} end
	abstract type AbstractFieldSpace{T,D} <: AbstractFieldMonochromatic{T,D} end

	abstract type AbstractCoefficient{T<:Real, L<:Union{AbstractFieldMonochromatic{T}, Number}, R<:Union{AbstractFieldMonochromatic{T}, Number}} end
	abstract type AbstractCoefficientCache{T<:Real} end
	abstract type AbstractOpticalComponent{T<:Real} end
	abstract type AbstractDetector{T<:Real} <: AbstractOpticalComponent{T} end

	abstract type AbstractModes{T<:Real} end
	abstract type AbstractPropagationComponent{T<:Real} <: AbstractOpticalComponent{T} end

	include("JolabFunction.jl")
	include("Material.jl")
	include("ReferenceFrame.jl")

	include("AbstractFieldMonochromatic.jl")
	include("AbstractCoefficient.jl")
	include("AbstractOpticalComponent.jl")
	include("AbstractDetector.jl")

	include("FieldSpaceScalar.jl")
	include("AngularSpectrumScalar.jl")
	include("iterators.jl")

	include("numericalFunctions.jl");
	include("RefractiveIndexDatabase.jl")
	include("FieldModes.jl")
	include("PolichromaticField.jl")
	include("PointDetector.jl")
	include("ScatteringMatrix.jl")
	include("ScatteringCoefficient.jl")
	include("MultilayerStructure.jl")
	include("Lens.jl")
	include("SingleModeFibre.jl")
	include("Mirror.jl")
	include("FourierTransform.jl")
	include("interpolation.jl")
	include("CircularStepIndexFibre.jl")
	include("SpatialLightModulator.jl")
	include("RoughMultilayerStructure.jl")
	include("PhantomLayer.jl")
	include("Axicon.jl")
	include("AxiconFourier.jl")
	include("lightinteraction_recursive.jl")
	include("lightinteraction_recursivegridded.jl")
	include("PlotsRecipes.jl")
end
