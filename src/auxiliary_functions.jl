reverse_if_backward(::Type{Forward}, x) = x
reverse_if_backward(::Type{Backward}, x) = reverse(x)

import Base.!
Base.length(beam::Beam) = length(beam.modes)
Base.size(beam::Beam) = size(beam.modes)
(!)(::Type{Forward}) = Backward
(!)(::Type{Backward}) = Forward