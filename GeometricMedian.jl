using Statistics, LinearAlgebra

function GM(X)
  δ=1e-10;
  ϵ=1e-10;
  Loop=100;
  (N,d)=size(X);
  cⁿ=zeros(d);
  cˢ=zeros(d);
  
  cⁿ.=mean(X,dims=1)|>vec;
  for l in 1:Loop
    ds=map(i->1/max(δ, norm(X[i,:]-cⁿ)), 1:N);
    cˢ.=(mean(X .* ds, dims=1) |>vec) / sum(ds);
    
    if(norm(cˢ-cⁿ) < ϵ)
      cⁿ.=cˢ;
      break;
    else
      cⁿ=cˢ;
    end
  end
  return(cⁿ);
end