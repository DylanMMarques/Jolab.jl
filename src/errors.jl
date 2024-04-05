@enum ERROR_CODES::UInt64 begin
    INVALID_MEDIUM = 0
    INVALID_FRAME = 1
    INVALID_MEDIUM_COMPLEX = 2
    INVALID_WAVELENGTH = 3
    INVALID_FIELD_SAMPLING = 4
    INVALID_BEAM_TYPE = 5
end

import Base: <<

(<<)(int::Integer, code::ERROR_CODES) = int << UInt64(code)

function throw_error_msg(error_code)
    msg = "\n"
    n = 0
    while error_code != 0
        if isodd(error_code)
            code_n = 1 << (n - 1)
            msg *= "- " * errors_dict[ERROR_CODES(code_n), Forward] * '\n'
        end
        error_code >>= 1
        n += 1
    end
    throw(ArgumentError(msg))
end

const errors_dict = Dict(
    (INVALID_MEDIUM, Forward) => "The medium is invalid",
    (INVALID_MEDIUM, Backward) => "The medium is invalid",
    (INVALID_FRAME, Forward) => "The frame is invalid",
    (INVALID_FRAME, Backward) => "The frame is invalid",
    (INVALID_MEDIUM_COMPLEX, Forward) => "The medium cannot have a complex refractive index",
    (INVALID_MEDIUM_COMPLEX, Backward) => "The medium cannot have a complex refractive index",
)