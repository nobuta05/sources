using Distributions;
include("mem.jl");

function mlr(y, X)
  EPS = 10e-4;
  Loop = 100;
  N = length(y); p = 0;
  β = inv(X'*X)*X'*y;
  h = 1.144 * std(y-X*β) /(N^(1/5));
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
  function g(z,h)
    return Φh(z,h);
  end
  function dddg(z,h)
    return (3*z - z^3) * Φh(z,h) / (h^3);
  end
  
  # get number of data and variables
  if ndims(X) > 1
    (N, p) = size(X);
  else
    N = length(X);
    p = 1;
  end

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
    errdist = MixtureModel(
      Normal,
      map(x->(x,h), ε),
      map(i->1/N, 1:N)
    );
    m = mem(errdist);
    # z_i = \frac{ε_i - m}{h}
    z = (ε .- m) ./ h;
    G = g.(z, h);
    DDDG = dddg.(z,h);
    GW = mapreduce(i->X[i,:] .* G[i], hcat, 1:N);
    L = (GW * GW') ./ N;
    K = (X' * DDDG) ./ N;
    v2 = 1;
    h = ( (3*v2*(p+1))/(dot(K, inv(L)*K) .* N) )^(1/7);
    # h = sqrt(dot(ε, Q .* ε));
    # println("h: " * string(h));
    

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