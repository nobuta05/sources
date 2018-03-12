using Distributions;
include("mem.jl");

function mlr(y, X::Array{Float64,2}, β=nothing, h=nothing)
  N = length(y); p = 0;
  if β == nothing
    β = inv(X'*X)*X'*y;
  end
  if h == nothing
    h = 1.144 * std(y-X*β) /(N^(1/5));
  end
  Q = rand(N);

  function Φ(z)
    return exp(-0.5*z^2)/sqrt(2π)
  end
  function Φh(z,h)
    return Φ(z/h)/h;
  end
  function dΦh(z,h)
    return -z/(h^2) * Φh(z,h);
  end
  
  # get number of data and variables
  if ndims(X) > 1
    (N, p) = size(X);
  else
    N = length(X);
    p = 1;
  end
  EPS = 10e-8/N;
  Loop = 300;

  for l in 1:Loop
    if l == Loop
      println("Last");
    end
    ε = y - X*β;

    # e-step
    Q = Φh.(ε, h) |> vec;
    Q = Q ./ sum(Q);

    # m-step
    W = diagm(Q);
    ## optimize β
    βnxt = inv(X'*W*X)*X'*W*y;
    ε = y - X*βnxt;
    ## optimize h
    h = 1.144 * std(ε) /(N^(1/5));

    # judge
    if norm(βnxt-β) < EPS
      β = βnxt;
      break;
    else
      β = βnxt;
    end
  end

  return (β, Q, h);
end