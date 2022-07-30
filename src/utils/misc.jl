# make a list of parameters into a list of vectors of the same length
function consistent_params(ps...)
    lengths = Set(length.(ps))
    if length(lengths) == 1 # they're already the same length
        return ps
    elseif (length(lengths)==2 && 1 in lengths)
        len = [l for l in lengths if l!=1][1]
        new_ps = []
        for p in ps
            if length(p) == len
                push!(new_ps, p)
            else
                push!(new_ps, fill(p, len))
            end
        end
        return new_ps
    else
        throw(DimensionMismatch("Mixed length fields...they must be all the same length, or a combination of one length and length one ;)"))
    end
end

