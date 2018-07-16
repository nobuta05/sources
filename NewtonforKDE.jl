srand(100);
using StatsBase;
include("/home/nobuta05/Gits/sources/HSM.jl");

function NewtonforKDE2(xs)
  const Loop = 100;
  const ϵ = 10e-7;
  N=length(xs);
  h=1.144*std(xs)*N^(-0.2);
  xⁿ=HSM(xs);
  xˢ=0.0;

  function ϕ(x)
    return exp(-0.5*x*x)/sqrt(2π);
  end
  function ϕh(x)
    return ϕ(x/h)/h;
  end
  function dLdx(x0)
    # -dot((-xs.+x0)./h, ϕh.(-xs.+x0))/(N*h)
    s = 0.0;
    @simd for i in 1:N
      @inbounds s -= (x0-xs[i])/h * ϕh(x0-xs[i]);
    end
    s/(N*h)
  end
  function ddLdxx(x0)
    # -dot((-((-xs.+x0)./h).^2).+1, ϕh.(-xs.+x0) )/(N*h*h)
    s = 0.0;
    @simd for i in 1:N
      @inbounds s -= (1-((x0-xs[i])/h)^2) * ϕh(x0-xs[i]);
    end
    s/(N*h*h)
  end

  for i in 1:Loop
    xˢ=xⁿ-dLdx(xⁿ)/ddLdxx(xⁿ);
    if abs(xˢ-xⁿ) < ϵ
      xⁿ = xˢ;
      break;
    else
      xⁿ = xˢ;
    end
  end
  

  return xⁿ;
end