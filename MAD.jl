function MAD(xs)
  v=zeros(length(xs));
  v.=abs.(xs.-median(xs));
  # median(abs.(xs.-median(xs)))*1.4826
  median(v)*1.4826
end
