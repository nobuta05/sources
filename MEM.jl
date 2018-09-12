include(homedir()*"/Gits/sources/MAD.jl");
include(homedir()*"/Gits/sources/HSM.jl");

function MEM(xs;x=nothing)::Float64
  function ϕ(z::Float64)::Float64
    exp(-z*z*0.5) / sqrt(2π)
  end
  function ϕh(z::Float64,h::Float64)::Float64
    ϕ(z/h)/h
  end
  const ϵ=1e-8;
  const Loop=100;

  if x==nothing
    x=HSM(xs);
  end

  const N=length(xs);
  const h=1.144*MAD(xs)*N^(-0.2);
  xⁿ=x;
  for l in 1:Loop
    q=ϕh.(xs.-xⁿ,h);
    q=q./sum(q);
    xˢ=dot(xs,q);
    if norm(xˢ-xⁿ)<ϵ
      break;
    else
      xⁿ=xˢ;
    end
  end

  xⁿ
end
