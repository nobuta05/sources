using Random
Random.seed!(100);
include(homedir()*"/Gits/sources/HSM.jl");

function NewtonforKDE(xs; x=nothing)::Float64
  Loop = 100;
  ϵ = 10e-10;
  BORDER=maximum(abs.(xs));
  N=length(xs);
  h=max(1.144*mad(xs,normalize=true)*N^(-0.2), 1e-10);
  xⁿ= x==nothing ? HSM(xs) : x;
  xˢ=0.0;

  function ϕ(x)::Float64
    return exp(-0.5*x*x)/sqrt(2π);
  end
  function ϕh(x)::Float64
    return ϕ(x/h::Float64)/h;
  end
  function dLdx(x0::Float64)::Float64
    -dot((-xs.+x0)./h, ϕh.(-xs.+x0))/(N*h)
    #=
    s = 0.0;
    for i in 1:N
      @inbounds s -= (x0-xs[i])/h * ϕh(x0-xs[i]);
    end
    s/(N*h)
    =#
  end
  function ddLdxx(x0::Float64)::Float64
    -dot((-((-xs.+x0)./h).^2).+1, ϕh.(-xs.+x0) )/(N*h*h)
    #=
    s = 0.0;
    for i in 1:N
      @inbounds s -= (1-((x0-xs[i])/h)^2) * ϕh(x0-xs[i]);
    end
    s/(N*h*h)
    =#
  end

  for i in 1:Loop
    xˢ=xⁿ-dLdx(xⁿ)/ddLdxx(xⁿ);
    if abs(xˢ-xⁿ) < ϵ
      xⁿ = xˢ;
      break;
    elseif isnan(xˢ) || abs(xˢ)>BORDER
      xⁿ=0.;
      break;
    else
      xⁿ = xˢ;
    end
  end
  xⁿ
end
