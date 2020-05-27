using RecipesBase

function complextoplot(z::Array{<:Number, 3}; type = :abs2)
	if type == :abs2
		if size(z, 1) == 1
			zp = abs2.(z[1,:,:]);
		else
			zp = dotdim(z, z, 1);
		end
	elseif type == :abs
		if size(z, 1) == 1
			zp = norm.(z[1,:,:]);
		else
			zp = .√(dotdim(z, z, 1));
		end
	elseif type == :realx
		zp = real.(z[1,:,:]);
	elseif type == :realy
		if size(z, 1) == 1
			zp = real.(z[1,:,:]);
		else
			zp = real.(z[2,:,:]);
		end
	elseif type == :realz
		if size(z, 1) == 1
			zp = real.(z[1,:,:]);
		else
			zp = real.(z[3,:,:]);
		end
	elseif type == :imagx
		zp = imag.(z[1,:,:]);
	elseif type == :imagy
		if size(z, 1) == 1
			zp = imag.(z[1,:,:]);
		else
			zp = imag.(z[2,:,:]);
		end
	elseif type == :imagz
		if size(z, 1) == 1
			zp = imag.(z[1,:,:]);
		else
			zp = imag.(z[3,:,:]);
		end
	elseif type == :anglex
		zp = angle.(z[1,:,:]);
	elseif type == :angley
		if size(z, 1) == 1
			zp = angle.(z[1,:,:]);
		else
			zp = angle.(z[2,:,:]);
		end

	elseif type == :anglez
		if size(z, 1) == 1
			zp = angle.(z[1,:,:]);
		else
			zp = angle.(z[3,:,:]);
		end
	else
		error("Type not known. Use abs2, abs, realx, realy, realz, imagx, imagy, imagz, anglex, angley, anglez")
	end
	return zp;
end
@recipe f(angspe::FieldAngularSpectrum; type= :abs2) = (angspe.nsx_X, angspe.nsy_Y, complextoplot(angspe.e_SXY, type=type)')

@recipe f(space::FieldSpace; type= :abs2) = (space.x_X, space.y_Y, complextoplot(space.e_SXY, type=type)')
@recipe f(modes::CircularStepIndexModes{T}; type= :abs2) where {T<:Real} = begin
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
