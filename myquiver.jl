using Plotly
# scatterを使って無理やりベクトル場を表示する
# (xs[k],ys[k],vx[k],vy[k])が，(xs[k],ys[k])から伸びるベクトル(vx[k],vy[k])を表す．
# 1つ目の返り値が矢印を表示するscatterで
# 2つ目の返り値が終点を表示するscatter
function myquiver(xrng,yrng,f)
  mX = minimum(xrng);
  MX = maximum(xrng);
  mY = minimum(yrng);
  MY = maximum(yrng);
  # 最大長W
  W = min(MX-mX, MY-mY)*0.3;
  # ベクトルの最大長L
  L = maximum(mapreduce(y->map(x->norm(f(x,y)),xrng),vcat,yrng));

  retx=Array{Float64,1}(0);
  rety=Array{Float64,1}(0);
  px=Array{Float64,1}(0);
  py=Array{Float64,1}(0);
  for x in xrng
    for y in yrng
      v = f(x,y).*(W/L);
      retx = vcat(retx, [x,x+v[1],NaN]);
      rety = vcat(rety, [y,y+v[2],NaN]);
      px = vcat(px, [x+v[1]]);
      py = vcat(py, [y+v[2]]);
    end
  end
  
  return (scatter(x=retx,y=rety,mode="lines", lines=attr(width=1), showlegend=false), scatter(x=px,y=py,mode="markers",markers=attr(size=1), showlegend=false));
end

function myquiver(px,py,vx,vy; lengths=nothing)
  retx=Array{Float64,1}(0);
  rety=Array{Float64,1}(0);
  epx=Array{Float64,1}(0);
  epy=Array{Float64,1}(0);
  for i in 1:length(px)
    v = [vx[i],vy[i]];
    if lengths != nothing
      v = (v./norm(v)).*sqrt(lengths[i]);
    end
    retx = vcat(retx, [px[i], px[i]+v[1], NaN]);
    rety = vcat(rety, [py[i], py[i]+v[2], NaN]);
    epx = vcat(epx, [px[i]+v[1]]);
    epy = vcat(epy, [py[i]+v[2]]);
  end
  
  (scatter(x=retx,y=rety,mode="lines", lines=attr(width=1), showlegend=false), scatter(x=epx,y=epy,mode="markers",markers=attr(size=1), showlegend=false));
end

function myplotcomp(μ,U; lengths=nothing, name="")
  (p,K) = size(U);
  ret = Array{Float64,2}(0,2);
  for i in 1:K
    if lengths != nothing
      ret = vcat(ret, (μ-(U[:,i].*lengths[i]))');
      ret = vcat(ret, (μ+(U[:,i].*lengths[i]))');
      ret = vcat(ret, [NaN,NaN]');
    else
      ret = vcat(ret, (μ-U[:,i])');
      ret = vcat(ret, (μ+U[:,i])');
      ret = vcat(ret, [NaN,NaN]');
    end
  end
  scatter(x=ret[:,1], y=ret[:,2], mode="lines", name=name)
end