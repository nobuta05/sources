# Conjugate-Gradient method
# minimize f(x) = \frac{1}{2} x^{\top} Q x + b^{\top} x
# given (Q,b), then return x^{*}

function CG(Q,b, xinit=nothing)
  p = length(b);
  x = zeros(size(b));
  const Loop = 100;
  const EPS = 10e-4;
  d = zeros(p,Loop);

  if xinit != nothing
    x = xinit[:];
  end
  d[:,1] = Q*x+b;
  
  for k in 1:Loop
    gk = Q*x+b;
    # gkをもとにdkを作成する．
    # k=0ならばdk = gkとすればよい
    # k>0ならば，dk
    if k == 1
      d[:,k] = gk;
    else
      # そのうち確認すべきこと．
      # 共役勾配法では，Q-直交な基底を1ステップ毎に1つ追加していくはずである．
      # その際gkをもとにグラムシュミットの直交化から，dkを作成すると考えるのが自然だと思われる．
      # しかし参考文献を見ても，Q-直交化というには要素が不十分のように思われる．
      # 以下は自分が正しいと考える，gkのQ-直交化によるdkの作成
      beta = map(i -> -dot(d[:,i],Q*gk)/dot(d[:,i],Q*d[:,i]), 1:k-1);
      d[:,k] = gk + d[:,1:k-1]*beta;
    end
    alpha = -dot(gk,d[:,k])/dot(d[:,k],Q*d[:,k]);
    x = x + d[:,k].*alpha;
    if norm(Q*x+b) < EPS
      break;
    end
  end
  
  return x;
end