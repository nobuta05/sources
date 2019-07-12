using LinearAlgebra
include(homedir()*"/Gits/sources/GeometricMedian.jl");

function cPCA(inpX::Array{Float64,2}; isCentered=false, scale=false)::Tuple{Array{Float64,1}, Array{Float64,2}}
  (N, d) = size(inpX);
  gm = zeros(d);
  X = zeros(N, d);
  if !isCentered
    gm .= GM(inpX);
    X .= inpX .- gm';
  else
    X .= inpX;
  end
  if scale
    X .= (map(i->1/norm(X[i,:]), 1:N) |> Diagonal) * X;
  end

  eigret = eigen(Symmetric((X'*X)./N), 1:d);
  inds = sortperm(eigret.values, rev=true);

  return (gm, eigret.vectors[:,inds]);
end
