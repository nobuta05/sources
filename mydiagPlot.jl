include(homedir()*"/Gits/sources/NewtonforKDE.jl");
include(homedir()*"/Gits/sources/MAD.jl");

# switch...1の時. modvは最頻値
# 2の時. 中央値
# 3の時．平均値
function mydiagPlot(inpX::Array{Float64,2},W::Array{Float64,2}; k=1,switch=1,d=nothing,v₀=nothing)::Tuple{Array{Float64,1},Array{Float64,1}}
  (N,p)=size(inpX);
  v=zeros(p);
  y=inpX*W;
  mods=map(i->NewtonforKDE(y[:,i]),1:p);
  if v₀!=nothing
    v=v₀;
  elseif switch==1
    v=W*mods;
  elseif switch==3
    v=mean(inpX,1)'|>vec;
  else
  end
  X=inpX.-v';
  T=X*W[:,p-k+1:p];
  # σ²s=map(l->dot(T[:,l],T[:,l])/N, 1:k);
  if switch==1
    σ²s=map(l->MAD(T[:,l])^2,1:k);
  elseif d==nothing
    σ²s=map(l->dot(T[:,l].-mean(T[:,l]),T[:,l].-mean(T[:,l]))/N,1:k);
  else
    σ²s=d;
  end
  ods=map(i->norm(X[i,:]-W[:,p-k+1:p]*T[i,:]),1:N);
  sds=map(i->sqrt(sum(T[i,l]*T[i,l]/σ²s[l] for l in 1:k)),1:N);

  (sds,ods)
end
