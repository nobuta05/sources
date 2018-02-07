# X = readcsv("X.dat");
# y = readcsv("y.dat");

function LR(y,X)
  function Pi(x, β)
    return 1/(1+exp(-dot(x|>vec ,β|>vec)));
  end
  
  function gradlogL(y, X, β)
    return X'*(y - map(i -> Pi(X[i,1:end], β), 1:size(X,1)));
  end
  
  function H(X,β)
    return -X' * diagm(map(i -> Pi(X[i,1:end], β), 1:size(X,1))) * X;
  end

  β = zeros(size(X,2),1);
  EPS = 10e-5;

  for i in 1:100
    nxtβ = β - inv(H(X,β)) * gradlogL(y, X, β);
    if norm(nxtβ-β) < EPS
      β = nxtβ
      break
    end
    β = nxtβ
  end
  pred = map(i -> Pi(X[i,1:end], β), 1:size(X,1))
  ind = map(i -> pred[i]>=0.5?1:0, 1:length(pred))
  return (β, pred);
end