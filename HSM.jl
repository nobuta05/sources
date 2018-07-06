# Half-Sample Mode
function HSM(xs)
  if ndims(xs) != 1
    println("Error");
  return nothing;
  end

  N = length(xs);
  if N==1
    return xs[1];
  elseif N==2
    return mean(xs);
  else
    ys = sort(xs);
    M = Int(ceil(N/2));
    w = ys[N] - ys[1];
    stid = 0;
    for id in 1:Int(floor(N/2))+1
      if (ys[id+M-1] - ys[id]) < w
        stid = id;
        w = ys[id+M-1] - ys[id];
      end
    end
    return HSM(ys[stid:stid+M-1]);
  end
end