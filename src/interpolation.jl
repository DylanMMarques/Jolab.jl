@inline function fitrectangletoplane!(coefs_4::AbstractArray{T}, matrixtemp_44::AbstractArray{T}, x_4::AbstractArray{A}, y_4::AbstractArray{A}, z_4::AbstractArray{T}) where {T,A}
	@inbounds @simd for i in 1:4
		matrixtemp_44[i] = x_4[i] * y_4[i];
		matrixtemp_44[4 + i] = x_4[i];
		matrixtemp_44[8 + i] = y_4[i];
		matrixtemp_44[12 + i] = one(T);
	end
	matrixtemp_44 .= inv(matrixtemp_44);
	mul!(coefs_4, matrixtemp_44, vec(z_4));
end

@inline function fittriangletoplane!(coefs_3::AbstractVector{T}, matrixtemp_33::AbstractArray{T}, x_3::AbstractVector{A}, y_3::AbstractVector{A}, z_3::AbstractVector{T}) where {T,A}
	@inbounds @simd for i in 1:3
		matrixtemp_33[i] = x_3[i];
		matrixtemp_33[3 + i] = y_3[i];
		matrixtemp_33[6 + i] = one(T);
	end
	matrixtemp_33 .= inv(matrixtemp_33);
	mul!(coefs_3, matrixtemp_33, vec(z_3));
end

@inline function pointinsiderectangle(x_4::AbstractArray{<:Real}, y_4::AbstractArray{<:Real}, x::Real, y::Real)
	vx1 = x_4[1] - x;
	vy1 = y_4[1] - y;
	#abs(vx1) < MIN_VAL && (vx1 = vx1 / abs(vx1) * MIN_VAL)
	#abs(vy1) < MIN_VAL && (vy1 = vy1 / abs(vy1) * MIN_VAL)
	ind = [1;2;4;3;1]
	boolClockwise = true
	boolAntiClockwise = true
	@inbounds for i in eachindex(ind)
		vx2 = vx1;
		vy2 = vy1;
		vx1 = x_4[ind[i]] - x;
		vy1 = y_4[ind[i]] - y;

		# Check if all position are the same
		if abs(vx1) < @tol
			if abs(vy1) < @tol
				return true;
			end
		end

		# Check if position are inline with the square
		Δx = x_4[ind[i%5+1]] - x_4[ind[i]];
		Δy = y_4[ind[i%5+1]] - y_4[ind[i]];
		ty = (y - y_4[ind[i]]) / Δy
		if Δx < @tol
			((abs(vx1) < @tol) && (0 < ty < 1)) && (return true)
		end
		tx = (x - x_4[ind[i]]) / Δx
		if Δy < @tol
			((abs(vy1) < @tol) && (0 < tx < 1)) &&  (return true)
		end
		((abs(tx - ty) < @tol) && (0 < tx < 1)) && (return true)

		(vx1 * vy2 - vx2 * vy1 < 0) && (boolClockwise = false)
		(vx1 * vy2 - vx2 * vy1 > 0) && (boolAntiClockwise = false)
	end
	return boolAntiClockwise || boolClockwise
end

@inline @inbounds function pointinsidetriangle(x_3::AbstractArray{<:Real}, y_3::AbstractArray{<:Real}, x::Real, y::Real)
	MIN_VAL = 1E-8;
	vx1 = x_3[1] - x;
	vy1 = y_3[1] - y;
	#abs(vx1) < MIN_VAL && (vx1 = vx1 / abs(vx1) * MIN_VAL)
	#abs(vy1) < MIN_VAL && (vy1 = vy1 / abs(vy1) * MIN_VAL)

	@inbounds for i in [1; 2; 3; 1]
		vx2 = vx1;
		vy2 = vy1;
		vx1 = x_3[i] - x;
		vy1 = y_3[i] - y;

		if abs(vx1) < MIN_VAL
			if abs(vy1) < MIN_VAL
				return true;
			else
			#	vx1 = vx1 / abs(vx1) * MIN_VAL;
			end
		end
		#abs(vy1) < MIN_VAL && (vy1 = vy1 / abs(vy1) * MIN_VAL)

		(vx1 * vy2 - vx2 * vy1 < 0) && return false
	end
	return true
end

function interpolateirregulargridtogrid!(z_SCD::AbstractArray{Complex{T}, 3}, x_XY::AbstractArray{N}, y_XY::AbstractArray{N}, z_SXY::AbstractArray{Complex{T},3}, xI_C::AbstractArray{<:Real}, yI_D::AbstractArray{<:Real}) where {T<:Real, N<:Real}
	size(x_XY) == size(y_XY) == size(z_SXY)[2:3] || error("size(x_XY), size(y_XY) and size(z_SZY)[2:3] must be the same");
	sizeS, sizeX, sizeY = size(z_SXY);
	sizeC = length(xI_C);
	sizeD = length(yI_D);

	sizeS == size(z_SCD, 1) || error("size(z_SCD,1) and size(z_SXY,1) must be the same")
	size(z_SCD,2) == sizeC && size(z_SCD,3) == sizeD || error("size(z_SCD, 2) must be equal to length(xI_C) and size(z_SCD,3) must be equal to lenght(yI_D)");

	@inbounds @simd for i in eachindex(z_SCD)
		z_SCD[i] = zero(Complex{T});
	end
	xIstep = xI_C[2] - xI_C[1];
	yIstep = yI_D[2] - yI_D[1];
	coefsAbs_4S = MVector{4*sizeS,T}(undef);
	coefsAngle_4S = MVector{4*sizeS,T}(undef);
	xRect_4 = MVector{4,N}(undef);
	yRect_4 = MVector{4,N}(undef);
	zRect_4 = MVector{4,T}(undef);
	tempMatrix_44 = MArray{Tuple{4,4}, T, 2}(undef);
	@inbounds for iY in 2:sizeY
		for iX in 2:sizeX
			auxMin = min(x_XY[iX-1, iY-1], x_XY[iX-1,iY], x_XY[iX,iY-1], x_XY[iX,iY]);
			auxMax = max(x_XY[iX-1, iY-1], x_XY[iX-1,iY], x_XY[iX, iY-1], x_XY[iX, iY]);
			if xIstep > 0
				iCmin = floor(Int64, (auxMin - xI_C[1]) / xIstep)
				iCmax = ceil(Int64, (auxMax - xI_C[1]) / xIstep);
			else
				iCmax = ceil(Int64, (auxMin - xI_C[1]) / xIstep)
				iCmin = floor(Int64, (auxMax - xI_C[1]) / xIstep);
			end

			iCmin > sizeC || iCmax < 1 && continue
			iCmax > sizeC && (iCmax = sizeC)
			iCmin < 1 && (iCmin = 1)

			auxMin = min(y_XY[iX-1, iY-1], y_XY[iX-1,iY], y_XY[iX,iY-1], y_XY[iX,iY]);
			auxMax = max(y_XY[iX-1, iY-1], y_XY[iX-1,iY], y_XY[iX, iY-1], y_XY[iX, iY]);

			if yIstep > 0
				iDmin = floor(Int64, (auxMin - yI_D[1]) / yIstep)
				iDmax = ceil(Int64, (auxMax - yI_D[1]) / yIstep);
			else
				iDmax = ceil(Int64, (auxMin - yI_D[1]) / yIstep)
				iDmin = floor(Int64, (auxMax - yI_D[1]) / yIstep);
			end

			iDmin > sizeD || iDmax < 1 && continue;
			iDmax > sizeD && (iDmax = sizeD)
			iDmin < 1 && (iDmin = 1)

			xRect_4[1] = x_XY[iX-1, iY-1];
			xRect_4[2] = x_XY[iX-1,iY];
			xRect_4[3] = x_XY[iX,iY-1];
			xRect_4[4] = x_XY[iX,iY];
			yRect_4[1] = y_XY[iX-1,iY-1];
			yRect_4[2] = y_XY[iX-1,iY];
			yRect_4[3] = y_XY[iX,iY-1];
			yRect_4[4] = y_XY[iX,iY];

			for iS in 0:(sizeS-1)
				zRect_4[1] = abs(z_SXY[iS+1,iX-1,iY-1]);
				zRect_4[2] = abs(z_SXY[iS+1,iX-1,iY]);
				zRect_4[3] = abs(z_SXY[iS+1,iX,iY-1]);
				zRect_4[4] = abs(z_SXY[iS+1,iX,iY]);
				coefs1_4 = @view coefsAbs_4S[1+iS*4:4+iS*4];
				fitrectangletoplane!(coefs1_4, tempMatrix_44, xRect_4, yRect_4, zRect_4);

				zRect_4[1] = angle(z_SXY[iS+1,iX-1,iY-1]);
				zRect_4[2] = angle(z_SXY[iS+1,iX-1,iY]);
				zRect_4[3] = angle(z_SXY[iS+1,iX,iY-1]);
				zRect_4[4] = angle(z_SXY[iS+1,iX,iY]);
				unwrap!(zRect_4);

				coefs1_4 = @view coefsAngle_4S[1+iS*4:4+iS*4];
				fitrectangletoplane!(coefs1_4, tempMatrix_44, xRect_4, yRect_4, zRect_4);
			end
			for iD in iDmin:iDmax
				for iC in iCmin:iCmax
					if pointinsiderectangle(xRect_4, yRect_4, xI_C[iC], yI_D[iD])
						# z_SCD[1,iC,iD] += one(T);
						for iS in 0:(sizeS-1)
							auxAbs = coefsAbs_4S[1+4*iS] * xI_C[iC] * yI_D[iD] + coefsAbs_4S[2+4*iS] * xI_C[iC] + coefsAbs_4S[3+4*iS] * yI_D[iD] + coefsAbs_4S[4 + 4 * iS]
							auxAngle = coefsAngle_4S[1+4*iS] * xI_C[iC] * yI_D[iD] + coefsAngle_4S[2+4*iS] * xI_C[iC] + coefsAngle_4S[3+4*iS] * yI_D[iD] + coefsAngle_4S[4 + 4 * iS]
							z_SCD[iS+1,iC,iD] = auxAbs * exp(im * auxAngle)
						end
					end
				end
			end
		end
	end
end

function interpolateirregulargridtogrid!(z_SCD::AbstractArray{T, 3}, x_XY::AbstractArray{N}, y_XY::AbstractArray{N}, z_SXY::AbstractArray{T,3}, xI_C::AbstractArray{<:Real}, yI_D::AbstractArray{<:Real}) where {T<:Number, N<:Real}
	size(x_XY) == size(y_XY) == size(z_SXY)[2:3] || error("size(x_XY), size(y_XY) and size(z_SZY)[2:3] must be the same");
	sizeS, sizeX, sizeY = size(z_SXY);
	sizeC = length(xI_C);
	sizeD = length(yI_D);

	sizeS == size(z_SCD, 1) || error("size(z_SCD,1) and size(z_SXY,1) must be the same")
	size(z_SCD,2) == sizeC && size(z_SCD,3) == sizeD || error("size(z_SCD, 2) must be equal to length(xI_C) and size(z_SCD,3) must be equal to lenght(yI_D)");

	@inbounds @simd for i in eachindex(z_SCD)
		z_SCD[i] = zero(T);
	end
	xIstep = xI_C[2] - xI_C[1];
	yIstep = yI_D[2] - yI_D[1];
	coefs1_4S = MVector{4*sizeS,T}(undef);
	coefs2_4S = MVector{4*sizeS,T}(undef);
	xRect_4 = MVector{4,N}(undef);
	yRect_4 = MVector{4,N}(undef);
	zRect_4 = MVector{4,T}(undef);
	tempMatrix_44 = MArray{Tuple{4,4}, T, 2}(undef);
	@inbounds for iX in 2:sizeX
		for iY in 2:sizeY
			aux = min(x_XY[iX-1, iY-1], x_XY[iX-1,iY]);
			iCmin = floor(Int64, (aux - xI_C[1]) / xIstep)
			aux = max(x_XY[iX, iY-1], x_XY[iX, iY]);
			iCmax = ceil(Int64, (aux - xI_C[1]) / xIstep);

			iCmin > sizeC || iCmax < 1 && continue
			iCmax > sizeC && (iCmax = sizeC)
			iCmin < 1 && (iCmin = 1)

			aux = min(y_XY[iX-1, iY-1], y_XY[iX,iY-1]);
			iDmin = floor(Int64, (aux - yI_D[1]) / yIstep)
			aux = max(y_XY[iX-1, iY], y_XY[iX, iY]);
			iDmax = ceil(Int64, (aux - yI_D[1]) / yIstep);

			iDmin > sizeD || iDmax < 1 && continue;
			iDmax > sizeD && (iDmax = sizeD)
			iDmin < 1 && (iDmin = 1)

			xRect_4[1] = x_XY[iX-1, iY-1];
			xRect_4[2] = x_XY[iX-1,iY];
			xRect_4[3] = x_XY[iX,iY-1];
			xRect_4[4] = x_XY[iX,iY];
			yRect_4[1] = y_XY[iX-1, iY-1];
			yRect_4[2] = y_XY[iX-1,iY];
			yRect_4[3] = y_XY[iX,iY-1];
			yRect_4[4] = y_XY[iX,iY];

	#		xTrian1 = @view xRect_4[3:-1:1]; yTrian1 = @view yRect_4[3:-1:1];
	#		xTrian2 = @view xRect_4[2:4]; yTrian2 = @view yRect_4[2:4];
			for iS in 0:(sizeS-1)
				zRect_4[1] = z_SXY[iS+1,iX-1,iY-1];
				zRect_4[2] = z_SXY[iS+1,iX-1,iY];
				zRect_4[3] = z_SXY[iS+1,iX,iY-1];
				zRect_4[4] = z_SXY[iS+1,iX,iY];
	#			zTrian1 = @view zRect_4[1:3];
	#			zTrian2 = @view zRect_4[4:-1:2];
				coefs1_4 = @view coefs1_4S[1+iS*4:4+iS*4];
				#coefs2_4 = @view coefs2_3S[1+iS*4:4+iS*4];

				fitrectangletoplane!(coefs1_4, tempMatrix_44, xRect_4, yRect_4, zRect_4);
				#fittriangletoplane!(coefs2_4, tempMatrix_44, xTrian2, yTrian2, zTrian2);
			end
			for iC in iCmin:iCmax
				for iD in iDmin:iDmax
					if pointinsiderectangle(xRect_4, yRect_4, xI_C[iC], yI_D[iD])
						for iS in 0:(sizeS-1)
							z_SCD[iS+1,iC,iD] = coefs1_4S[1+4*iS] * xI_C[iC] * yI_D[iD] + coefs1_4S[2+4*iS] * xI_C[iC] + coefs1_4S[3+4*iS] * yI_D[iD] + coefs1_4S[4 + 4 * iS];
						end
					end
				end
			end
		end
	end
end

function interpolatetogrid!(z_SCD::AbstractArray{T,3}, x_XY::AbstractArray{<:Real,2}, y_XY::AbstractArray{<:Real,2}, z_SXY::AbstractArray{T,3}, xI_C::AbstractVector{<:Real}, yI_D::AbstractVector{<:Real}) where T
	interpolateirregulargridtogrid!(z_SCD, x_XY, y_XY, z_SXY, xI_C, yI_D);
	return z_SCD;
end

function interpolatetogrid!(z_CD::AbstractArray{T,2}, x_XY::AbstractArray{<:Real,2}, y_XY::AbstractArray{<:Real,2}, z_XY::AbstractArray{T,2}, xI_C::AbstractVector{<:Real}, yI_D::AbstractVector{<:Real}) where T
	z_SCD = reshape(z_CD, (1, size(z_CD)...));
	z_SXY = reshape(z_XY, (1, size(z_XY)...));
	interpolatetogrid!(z_SCD, x_XY, y_XY, z_SXY, xI_C, yI_D);
	return z_SCD;
end

function interpolatetogrid(x_XY::AbstractArray{<:Real,2}, y_XY::AbstractArray{<:Real,2}, z_XY::AbstractArray{T,2}, xI_C::AbstractVector{<:Real}, yI_D::AbstractVector{<:Real}) where {T}
	z_CD = zeros(T, length(xI_C), length(yI_D));
	interpolatetogrid!(z_CD, x_XY, y_XY, z_XY, xI_C, yI_D);
	return z_CD;
end

function interpolatetogrid(x_XY::AbstractArray{<:Real,2}, y_XY::AbstractArray{<:Real,2}, z_SXY::AbstractArray{T,3}, xI_C::AbstractVector{<:Real}, yI_D::AbstractVector{<:Real}) where {T}
	z_SCD = zeros(T, size(z_SXY,1), length(xI_C), length(yI_D));
	interpolatetogrid!(z_SCD, x_XY, y_XY, z_SXY, xI_C, yI_D);
	return z_SCD
end


function interpolatetogrid(x_XY::AbstractArray{<:Real,2}, y_XY::AbstractArray{<:Real,2}, z_SXY::AbstractArray{<:Number,3})
	(sizeX, sizeY) = size(x_XY)
	(sizeX, sizeY) != size(y_XY) && error("wrong sizes")

	if (x_XY[2] - x_XY[1]) > 0
		x_X = range(minimum(view(x_XY,1,:)), stop = maximum(view(x_XY,sizeX,:)), length = sizeX)
	else
		x_X = range(maximum(view(x_XY,1,:)), stop = minimum(view(x_XY,sizeX,:)), length = sizeX)
	end
	if (y_XY[1,2] - y_XY[1,1]) > 0
		y_Y = range(minimum(view(y_XY,:,1)), stop = maximum(view(y_XY,:,sizeY)), length = sizeY)
	else
		y_Y = range(maximum(view(y_XY,:,1)), stop = minimum(view(y_XY,:,sizeY)), length = sizeY)
	end
	return (x_X, y_Y, interpolatetogrid(x_XY, y_XY, z_SXY, x_X, y_Y))
end
