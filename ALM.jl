srand(100);

A = rand(4,4);
A = (A+A')./2;
(_,U) = eig(A);
Σ = U*diagm([1,2,3,4])*U';

function ALM(v0=nothing)
  λ = 1;
  μ = 1;
  const α = 1.5;
end
function ALM(λ,μ,v0=nothing)
end