include(homedir()*"/Gits/sources/MAD.jl");
include(homedir()*"/Gits/sources/HSM.jl");
include(homedir()*"/Gits/sources/NewtonforKDE.jl");
include(homedir()*"/Gits/sources/MEM.jl");

function uniMode(xs::Array{Float64,1})::Float64
  function ϕ(z::Float64)::Float64
    exp(-z*z*0.5) / sqrt(2π)
  end
  function ϕh(z::Float64,h::Float64)::Float64
    ϕ(z/h)/h
  end
  
  N=length(xs);
  h=1.144*MAD(xs)*N^(-0.2);
  hsm=HSM(xs);
  ntn=NewtonforKDE(xs);
  if sum(ϕh.(xs.-hsm,h)) > sum(ϕh.(xs.-ntn,h))
    MEM(xs,x=hsm)
  else
    ntn
  end
end