using LinearAlgebra, Statistics, StatsBase

function MLR(y::Array{Float64,1}, X::Array{Float64,2})
  ϵ=1e-10;
  δ=1e-10;
  Loop=100;
  function ϕ(z::Float64)::Float64
    return(exp(-0.5*z^2)/sqrt(2π));
  end
  function ϕh(z::Float64, h::Float64)::Float64
    return(ϕ(z/h)/h);
  end

  (N,d)=size(X);
  βⁿ=zeros(d);
  ϵs=zeros(N);

  βⁿ.=(X'*X) \ (X'*y);
  ϵs.=y-X*βⁿ;
  h=1.144*mad(ϵs, normalize=true)*N^(-0.2);

  for l in 1:Loop
    qs=ϕh.(ϵs, h);
    qs.=qs./sum(qs);

    βˢ=(X'*(X.*qs)) \ (X'*(qs.*y));
    ϵs.=y-X*βˢ;
    h=1.144*mad(ϵs, normalize=true)*N^(-0.2);
    if norm(βˢ-βⁿ)<ϵ
      βⁿ.=βˢ;
      break;
    else
      βⁿ.=βˢ;
    end
  end
  return(Dict(
    :β => βⁿ
  ));
end
