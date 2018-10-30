# 線分を描画するための関数．X=(x1 x2 ... xN)ᵀ, Y=(y1 y2 ... yN)ᵀ を引数とする．
# x1とy1, x2とy2，... xNとyNで線分を結ぶ．
# 返り値として(l×2)行列Zを返す．scatter(x=Z[:,1], y=Z[:,2], mode="lines")とすれば描画できるようにする．

function mylineseg(X,Y)
  N=size(X)[1];
  Z=mapreduce(i->vcat(X[i,:]',Y[i,:]',[nothing nothing]), vcat, 1:N)
end
