srand(100);
include(homedir()*"/Gits/sources/HSM.jl");
include(homedir()*"/Gits/sources/MAD.jl");

function NewtonforKDE(xs::Array{Float64,1})::Float64
  const Loop = 100;
  const ϵ = 10e-7;
  N=length(xs);
  h=1.144*MAD(xs)*N^(-0.2);
  xⁿ=HSM(xs);
  xˢ=0.0;

  function ϕ(x::Float64)::Float64
    return exp(-0.5*x*x)/sqrt(2π);
  end
  function ϕh(x::Float64)::Float64
    return ϕ(x/h::Float64)/h;
  end
  function dLdx(x0::Float64)::Float64
    # -dot((-xs.+x0)./h, ϕh.(-xs.+x0))/(N*h)
    s = 0.0;
    @simd for i in 1:N
      @inbounds s -= (x0-xs[i])/h * ϕh(x0-xs[i]);
    end
    s/(N*h)
  end
  function ddLdxx(x0::Float64)::Float64
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
    elseif xˢ>maximum(abs.(xs))
      xⁿ=0.;
      break;
    else
      xⁿ = xˢ;
    end
  end
  

  return xⁿ;
end