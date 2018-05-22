using Distributions;

function phi(x,t)
  return exp(-0.5*x^2/t)/sqrt(2π*t)
end

dΦ = Array{Function,1}(10);
dΦ[1] = (x,t) -> -x*phi(x,t)/t
dΦ[2] = (x,t) -> ( (x^2)/t-1 )*phi(x,t)/t
dΦ[3] = (x,t) -> x*( 3-(x^2)/t )*phi(x,t)/(t^2)
dΦ[4] = (x,t) -> ( 3-6*(x^2)/t+((x^2)/t)^2 )*phi(x,t)/(t^2)
dΦ[5] = (x,t) -> x*( -15+10*(x^2/t)-((x^2)/t)^2 )*phi(x,t)/(t^3)
dΦ[6] = (x,t) -> ( -15+45*x^2/t-15*((x^2)/t)^2+((x^2)/t)^3 )*phi(x,t)/(t^3)
dΦ[7] = (x,t) -> x*( 105-105*x^2/t+21*((x^2)/t)^2-((x^2)/t)^3 )*phi(x,t)/(t^4)
dΦ[8] = (x,t) -> ( 105-420*x^2/t+210*((x^2)/t)^2-28*((x^2)/t)^3+((x^2)/t)^4 )*phi(x,t)/(t^4)
dΦ[9] = (x,t) -> x*( -945+1260*x^2/t-378*((x^2)/t)^2+36*((x^2)/t)^3-((x^2)/t)^4 )*phi(x,t)/(t^5)
dΦ[10] = (x,t) -> ( -945+4725*x^2/t-3150*((x^2)/t)^2+630*((x^2)/t)^3-45*((x^2)/t)^4+((x^2)/t)^5 )*phi(x,t)/(t^5)

function diffusion_h(xs,ε=0.000001)
  const loop = 100;
  const L = 5;
  t = ε;
  ξ = 0.9;
  N = length(xs);
  for k in 1:loop
    if k==loop
      println("Last");
    end
    ts = zeros(L);
    ts[L] = t;
    ints = zeros(L);
    for l in L:-1:2
      # calc (integral of f^(l))^2
      ints[l] = sum(dΦ[2*l].(xs.-xs', 2*ts[l]))*((-1)^l)/(N^2);
      ts[l-1] = ( ( sqrt(2)+0.5^(l-1) )*prod(1:2:(2l-3)) / ( 3*N*sqrt(π)*ints[l] ) )^( 2/( 3+2(l-1) ) );
    end
    tnxt = ξ*ts[1];
    if norm(tnxt - t) < ε && k > 10
      t = tnxt;
      break;
    else
      t = tnxt;
    end
  end
  return sqrt(t);
end

function example()
  d = MixtureModel(Normal,[(0.0,1), (7.0,7.0)],[0.2,0.8]);
  return d;
end