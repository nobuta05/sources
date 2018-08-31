function MAD(x::Array{Float64,1})::Float64
  median(abs.(x.-median(x)))*1.4826
end