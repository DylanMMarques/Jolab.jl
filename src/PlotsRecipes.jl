using RecipesBase

function complextoplot(z::Array{<:Number, 2}; type = :abs2)
	if type == :abs2
		zp = abs2.(z)
	elseif type == :abs
		zp = abs.(z)
	elseif type == :real
		zp = real.(z)
	elseif type == :imag
		zp = imag.(z)
	elseif type == :angle
		zp = angle.(z)
	else
		error("Type not known. Use abs2, abs, real, imag, angle")
	end
	return zp;
end
@recipe f(angspe::AbstractFieldAngularSpectrum; type= :abs2) = (angspe.nsx_X, angspe.nsy_Y, complextoplot(reshape(angspe.e_SXY,length(angspe.nsx_X), length(angspe.nsy_Y)), type=type)')

@recipe f(space::FieldSpaceScalar; type= :abs2) = (space.x_X, space.y_Y, complextoplot(reshape(space.e_SXY,length(space.x_X), length(space.y_Y)), type=type)')
@recipe function f(modes::CircularStepIndexModes{T}; type= :abs2) where {T<:Real}
	sizeX = 200;
	sizeY = 200;
	x_X = range(- 2 * modes.r, 2 * modes.r, length = sizeX)
	y_Y = range(- 2 * modes.r, 2 * modes.r, length = sizeY)
	e_SXY = Array{Complex{T},3}(undef, 1, sizeX, sizeY)
	circularstepindex_modefield!(e_SXY, modes.r, modes.ncore, modes.na, modes.λ, modes.m[1], modes.β[1], modes.C[1], modes.D[1], 1, x_X, y_Y, 0)
	(x_X, y_Y, complextoplot(e_SXY, type = type)')
end

@recipe function f(ref::ReferenceFrame{T}; length = 1::Real) where {T<:Real}
	seriestype = :path3d
	(x1 , y1, z1) = rotatecoordinatesto(length, 0, 0, ref.θ, ref.ϕ)
	(x2 , y2, z2) = rotatecoordinatesto(0, length, 0, ref.θ, ref.ϕ)
	(x3 , y3, z3) = rotatecoordinatesto(0, 0, length, ref.θ, ref.ϕ)
	[([ref.x; ref.x + x1], [ref.y; ref.y + y1], [ref.z; ref.z + z1]),
	([ref.x; ref.x + x2], [ref.y; ref.y + y2], [ref.z; ref.z + z2]),
	([ref.x; ref.x + x3], [ref.y; ref.y + y3], [ref.z; ref.z + z3])]
end
