using StatsBase, Distributions;
include("/home/nobuta05/Gits/sources/diffusion.jl");
include("/home/nobuta05/Gits/sources/mem.jl");

function generate(N=500)
  β = [0;2];
  xd = Uniform(-5,5);
  mode_ed = 0.02426779526511704;
  ed = MixtureModel(Normal,[(0.0-mode_ed,1.0),(5.0-mode_ed,5.0)],[0.3,0.7]);
  X = hcat(ones(N), rand(xd,N));
  return (β,X,ed);
end

function estimated(εs,h)
  n = length(εs);
  return MixtureModel(Normal, map(ε->(ε,h),εs),map(i->1/n, 1:n))
end

function K(z)
  return exp(-z*z*0.5)/sqrt(2π);
end
function g(ε,m,hε)
  return K((ε-m)/hε)/hε;
end
function d2g(ε,m,hε)
  return K((ε-m)/hε)*( ((ε-m)/hε)^2-1 )/(hε^3);
end
function d3g(ε,m,hε)
  return K((ε-m)/hε)*( 3-((ε-m)/hε)^2 )*((ε-m)/hε)/(hε^4);
end
function hopt(εs,m,hε,X,p)
  N = length(εs);
  K = mapreduce(i->X[i,:].*d3g(εs[i],m,hε), +, 1:N)./N;
  L = mapreduce(i->X[i,:]*X[i,:]'.*g(εs[i],m,hε), +, 1:N)./N;
  return ((3*(p+1))/(N*dot(K,inv(L)*K)))^(1/7);
  end

function mlr_origin(y::Array{Float64,1}, X::Array{Float64,2}, βinit=nothing)
  (N,p) = size(X);
  β=zeros(p);
  if βinit == nothing
    β = inv(X'*X)*X'*y;
  else
    β += βinit;
  end
  ε = y-X*β;
  hε = 1.144 * std(ε) / (N^(0.2));
  m = mem(estimated(ε,hε));

  h = hopt(ε,m,hε,X,p);
  Q = zeros(N);
  W = zeros(N,N);
  βnxt = zeros(p);
  const EPS = 10e-8/N;
  const Loop = 500;

  function Φ(z)
    return exp(-0.5*z^2)/sqrt(2π)
  end
  function Φh(z,h)
    return Φ(z/h)/h;
  end

  for l in 1:Loop
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
    # hε = 1.144 * std(ε) / (N^(0.2));
    hε = diffusion_h(ε);
    m = mem(estimated(ε,hε));
    h = hopt(ε,m,hε,X,p);
    
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

function mlr(y::Array{Float64,1}, X::Array{Float64,2}, βinit=nothing)
  (N,p) = size(X);
  β=zeros(p);
  if βinit == nothing
    β = inv(X'*X)*X'*y;
  else
    β += βinit;
  end
  h = 1.144 * std(y-X*β) / (N^(0.2));

  Q = zeros(N);
  ε = zeros(N);
  W = zeros(N,N);
  βnxt = zeros(p);
  const EPS = 10e-8/N;
  const Loop = 500;

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

function mlr_mad(y::Array{Float64,1}, X::Array{Float64,2}, βinit=nothing)
  (N,p) = size(X);
  β=zeros(p);
  if βinit == nothing
    β = inv(X'*X)*X'*y;
  else
    β = βinit[:];
  end
  h = 1.144 * mad(y-X*β) / (N^(0.2));

  Q = zeros(N);
  ε = zeros(N);
  W = zeros(N,N);
  βnxt = zeros(p);
  const EPS = 10e-8/N;
  const Loop = 500;

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
    h = 1.144 * mad(ε) / (N^(0.2));

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


function mlr_hfix(y::Array{Float64,1}, X::Array{Float64,2}, h, βinit=nothing)
  (N,p) = size(X);
  β=zeros(p);
  if βinit == nothing
    β = inv(X'*X)*X'*y;
  else
    β = βinit[:];
  end
  
  Q = zeros(N);
  ε = zeros(N);
  W = zeros(N,N);
  βnxt = zeros(p);
  const EPS = 10e-8/N;
  const Loop = 500;

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

function mlr_hs(y::Array{Float64,1}, X::Array{Float64,2}, hsinit::Array{Float64,1}, βinit=nothing)
  hs = sort(hsinit,rev=true);
  (N,p) = size(X);
  β=zeros(p);
  if βinit == nothing
    β = inv(X'*X)*X'*y;
  else
    β = βinit[:];
  end
  
  Q = zeros(N);
  ε = zeros(N);
  W = zeros(N,N);
  βnxt = zeros(p);
  const EPS = 10e-8/N;
  const Loop = 500;

  function Φ(z)
    return exp(-0.5*z^2)/sqrt(2π)
  end
  function Φh(z,h)
    return Φ(z/h)/h;
  end

  for h in hs
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
    
      # judge
      if norm(βnxt-β) < EPS
        β .= βnxt;
        break;
      else
        β .= βnxt;
      end
    end
  end

  return (β, Q, hs);
end

function mlr_hs_comp(y::Array{Float64,1}, X::Array{Float64,2}, hsinit::Array{Float64,1}, βinit=nothing)

  hs = sort(hsinit,rev=true);
  (N,p) = size(X);
  β=zeros(p);
  if βinit == nothing
    β = inv(X'*X)*X'*y;
  else
    β = βinit[:];
  end
  
  Q = zeros(N);
  ε = zeros(N);
  W = zeros(N,N);
  βnxt = zeros(p);
  const EPS = 10e-8/N;
  const Loop = 500;

  function Φ(z)
    return exp(-0.5*z^2)/sqrt(2π)
  end
  function Φh(z,h)
    return Φ(z/h)/h;
  end

  for h in hs
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
    
      # judge
      if norm(βnxt-β) < EPS
        β .= βnxt;
        break;
      else
        β .= βnxt;
      end
    end

  end

  return (β, Q, hs);
end