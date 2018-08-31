# Half-Sample Mode
function HSM(xs::Array{Float64,1}; issort=false)::Float64
  N = length(xs);
  if N==1
    return xs[1];
  elseif N==2
    return 0.5*(xs[1]+xs[2]);
  else
    ys = issort ? xs : sort(xs);
    M = cld(N,2);
    w = ys[N] - ys[1];
    stid = 1;
    for id in 1:fld(N,2)+1
      if (ys[id+M-1] - ys[id]) < w
        stid = id;
        w = ys[id+M-1] - ys[id];
      end
    end
    return HSM(ys[stid:stid+M-1],issort=true);
  end
end