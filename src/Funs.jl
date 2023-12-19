module Funs

"""
    struct Fun{Tuple{Ts...},R}

A function that takes inputs with the types `Ts` and returns a result of type `R`.
"""
struct Fun{T,R,X}
    fun::X
    Fun{T,R,X}(fun::X) where {T<:Tuple,R,X} = new(fun)
end
export Fun

Fun{T,R}(fun) where {T,R} = Fun{T,R,typeof(fun)}(fun)

(fun::Fun{T,R})(xs...) where {T,R} = fun.fun(T(xs...))::R

end
