module AllParameters

import Base: getproperty, setproperty!

export AllParameter

struct AllParameter
    content::Dict{Symbol, Any}
end
AllParameter() = AllParameter(Dict{Symbol, Any}())

getproperty(p::AllParameter, name::Symbol) =
    if name == :content
        getfield(p, name)
    else haskey(p.content, name)
        p.content[name]
    end
function setproperty!(p::AllParameter, name::Symbol, value)
    @assert name != :content ":content is a reserved name."
    p.content[name] = value
end

end
