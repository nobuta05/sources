using Distributions, ArrayFire;
srand(100);

# ひとまずGMMを想定して実装
# 詳しい計算はpdf『MEMAlgorithm.pdf』を参照

# function mem(d::Sampleable, x=nothing, N::Int64=10)
function mem2(d, x=nothing, N::Int64=5)
  comps = components(d);
  prior = probs(d);
  K = length(prior);
  μs = AFArray(map(comp->mean(comp), comps));
  σs = AFArray(map(comp->std(comp), comps));
  σ2s = σs.^2;
  
  function mem_search(d, xinit)
    const Loop = 100;
    const EPS = 10e-5;
    x = 0.0;

    for i in 1:Loop
      # E-Step
      ps = AFArray(map(l->pdf(comps[l], x), 1:K) ./ pdf(d, x));
      
      # M-Step
      numerator = sum( (ps.*μs)./σ2s );
      denominator = sum(ps./σ2s);
      xnxt = numerator / denominator;

      if norm(xnxt - x) < EPS
        x = xnxt;
        println("loop:\t"*string(i));
        break;
      else
        x = xnxt;
      end
    end
    return x;
  end

  if x == nothing
    xs = rand(d, N);
    modes = mem_search.(d, xs);
    pdfs = pdf.(d, modes);
    ind = sortperm(pdfs, rev=true)[1];
    mode = modes[ind];
  else
    mode = mem_search(d,x);
  end

  return mode;
end