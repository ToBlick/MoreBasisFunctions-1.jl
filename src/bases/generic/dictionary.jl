
"""
This function is exactly like `eval_element_derivative`, but it evaluates
the antierivative of the element instead.
"""
function eval_element_antiderivative(dict::Dictionary, idx, x)
    idxn = native_index(dict, idx)
    @boundscheck checkbounds(dict, idxn)
    unsafe_eval_element_antiderivative1(dict, idxn, x)
end

function eval_element_antiderivative(dict::Dictionary, idx::LinearIndex, x)
    @boundscheck checkbounds(dict, idx)
    unsafe_eval_element_antiderivative1(dict, native_index(dict, idx), x)
end

function eval_element_extension_antiderivative(dict::Dictionary, idx, x)
    @boundscheck checkbounds(dict, idx)
    unsafe_eval_element_antiderivative(dict, native_index(dict, idx), x)
end

function unsafe_eval_element_antiderivative1(dict::Dictionary{S,T}, idx, x) where {S,T}
    in_support(dict, idx, x) ? unsafe_eval_element_antiderivative(dict, idx, x) : zero(T)
end

unsafe_eval_element_antiderivative(dict::Dictionary, idx, x) =
    unsafe_eval_element_antiderivative(dict, native_index(dict, idx), x)
