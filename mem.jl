using Distributions;
srand(100);

# ひとまずGMMを想定して実装
# 詳しい計算はpdf『MEMAlgorithm.pdf』を参照
function mem(d::Sampleable, x=nothing, N::Int64=10)
  comps = components(d);
  prior = probs(d);
  K = length(prior);
  μs = map(comp->mean(comp), comps);
  σs = map(comp->std(comp), comps);
  
  function mem_search(d, x)
    Loop = 100;
    EPS = 10e-7;

    for i in 1:Loop
      # E-Step
      ps = map(l->pdf(comps[l], x), 1:K) ./ pdf(d, x);
      # M-Step
      numerator = sum(map(l->ps[l]*μs[l]/(σs[l]^2), 1:K));
      denominator = sum(map(l->ps[l]/(σs[l]^2), 1:K));
      x_nxt = numerator / denominator;

      #= Output
      sum_ps_logcomps(x) = dot(ps, log.(pdf.(comps, x)));
      sum_ps_logπs = dot(ps, log.(prior));
      sum_ps_logps = dot(ps, log.(ps));
      surrogate(x) = sum_ps_logcomps(x) + sum_ps_logπs - sum_ps_logps;
      xrng = linspace(-5,5,500);
      plot(xrng, log.(pdf.(d, xrng)), label="log pdf", ylim=(-10,-1));
      plot!(xrng, surrogate.(xrng), label="surrogate", ylim=(-10,-1));
      vline!([x], label="current", lw=2);
      vline!([x_nxt], label="optimum", lw=2);
      savefig("i_"*string(i)*".pdf")
      =#

      if norm(x_nxt - x) < EPS
        break;
      else
        x = x_nxt;
      end
    end
    return x;
  end

  if x == nothing
    xs = rand(d, N);
    modes = zeros(N);
    modes = map(x->mem_search(d,x), xs);
    pdfs = pdf.(d, modes);
    ind = sortperm(pdfs, rev=true)[1];
    mode = modes[ind];
  else
    mode = mem_search(d,x);
  end

  return mode;
end