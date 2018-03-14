function mlr(y::Array{Float64,1}, X::Array{Float64,2}, β=nothing)
  (N,p) = size(X);
  if β == nothing
    β = inv(X'*X)*X'*y;
  end
  h = 1.144 * std(y-X*β) / (N^(0.2));

  Q = zeros(N);
  ε = zeros(N);
  W = zeros(N,N);
  βnxt = zeros(p);
  const EPS = 10e-7/N;
  const Loop = 100;

  function Φ(z)
    return exp(-0.5*z^2)/sqrt(2π)
  end
  function Φh(z,h)
    return Φ(z/h)/h;
  end

  for l in 1:Loop
    ε .= y .- X*β;

    # e-step
    Q .= Φh.(ε, h);
    Q .= Q ./ sum(Q);

    # m-step
    # W .= diagm(Q);
    ## optimize β
    # βnxt .= (X'*W*X)\(X'*W*y);
    βnxt .= (X'*(X.*Q)) \ (X'*(Q.*y));
    ε .= y .- X*βnxt;
    ## optimize h
    h = 1.144 * std(ε) / (N^(0.2));

    # judge
    if norm(βnxt-β) < EPS
      β .= βnxt;
      break;
    else
      β .= βnxt;
    end
  end

  return (β, Q, h);
end