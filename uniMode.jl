include(homedir()*"/Gits/sources/HSM.jl");
include(homedir()*"/Gits/sources/NewtonforKDE.jl");
include(homedir()*"/Gits/sources/MEM.jl");

function uniMode(xs)::Float64
  function ϕ(z::Float64)::Float64
    exp(-z*z*0.5) / sqrt(2π)
  end
  function ϕh(z::Float64,h::Float64)::Float64
    ϕ(z/h)/h
  end

  N=length(xs);
  h=max(1.144*mad(xs,normalize=true)*N^(-0.2), 1e-10);
  hsm=HSM(xs);
  ntn=NewtonforKDE(xs,x=hsm);
  if sum(ϕh.(xs.-hsm,h)) > sum(ϕh.(xs.-ntn,h))
    MEM(xs,x=hsm)
  else
    ntn
  end
end
