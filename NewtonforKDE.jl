using Distributions
srand(100);

function NewtonforKDE(d, x=nothing, N::Int64=5)
  comps = components(d);
  N = length(comps);
  xs = map(comp->mean(comp), comps);
  h = std(comps[1]);

  function Φ(x)
    return exp(-0.5*x*x)/sqrt(2π);
  end
  function Φh(x)
    return Φ(x/h)/h;
  end
  function dL(x)
    return -mapreduce(xi->(x-xi)*Φh(x-xi),+,xs)/(N*h*h*pdf(d,x));
  end
  function ddL(x)
    a = mapreduce(
          xi->( 1-((x-xi)/h)^2 )*Φh(x-xi),
          +,
          xs
        );
    b = mapreduce(
          xi->((x-xi)/h) * Φh(x-xi),
          +,
          xs
    );
    return -(pdf(d,x)*a+b*b/N)/(N*h*h*pdf(d,x)^2);
  end

  function search(xinit)
    const Loop = 1000;
    const EPS = 10e-6;
    x = xinit;

    for l in 1:Loop
      xnxt = x - dL(x)/ddL(x);
      if norm(xnxt-x)<EPS
        x = xnxt;
        break;
      else
        x = xnxt;
      end
    end
    return x;
  end

  if x == nothing
    xs = rand(d, N);
    modes = search.(xs);
    pdfs = pdf.(d, modes);
    ind = sortperm(pdfs, rev=true)[1];
    mode = modes[ind];
  else
    mode = search(x);
  end

  return mode;
end