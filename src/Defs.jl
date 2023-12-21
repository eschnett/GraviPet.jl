module Defs

function lincom(x0::S, y0::T, x1::S, y1::T, x::S) where {S,T}
    return T(y0 .* (T(x - x1) ./ T(x0 - x1)) + y1 .* (T(x - x0) ./ T(x1 - x0)))
end
export lincom

end
