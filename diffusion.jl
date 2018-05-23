using Roots,StatsBase;

function mydct(xs)
  N = length(xs);

  weights = [1;exp.(collect(1:(N-1)).*(-(2*π*im)/(2*N)))];
  data = vcat(xs[1:2:N],xs[N:-2:2]);
  return real.(weights.*fft(data));
end

function fixed_point(t,N,I,a2)
  l = 7;
  f = 2*(π^(2*l))*sum( (I.^l) .* a2 .* exp.(-I.*(2*π*t)) );
  for s in l-1:-1:2
    K0 = prod(1:2:2*s-1)/sqrt(2*π);
    c = (1+(1/2)^(s+1/2))/3;
    time = (2*c*K0/N/f)^(2/(3+2*s));
    f = 2*(π^(2*s))*sum( (I.^s) .* a2 .* exp.(-I.*(2*π*time)) );
  end
  return t-(2*N*sqrt(π)*f)^(-2/5);
end

function diffusion_h(xs)
  const n = 2^14;
  m = minimum(xs);
  M = maximum(xs);
  range = M-m;
  MIN = m - range/10;
  MAX = M + range/10;
  R = MAX - MIN;
  dx = R/n;
  xmesh = collect(0:dx:R).+MIN;
  N = length(xs);
  fitting = fit(Histogram,xs,xmesh,closed=:left);
  initial_data = fitting.weights;
  initial_data = initial_data./sum(initial_data);
  a = mydct(initial_data);
  I = collect(1:n-1).^2;
  a2 = (a[2:end]./2).^2;
  t_star = find_zero(t->fixed_point(t,N,I,a2), 0.1);
  return sqrt(t_star)*R;
end