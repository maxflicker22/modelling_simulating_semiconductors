function [Rn, Rp, Gn, Gp, II] = generation_recombination_model (x, n, p, E, Jn, Jp, tn, tp, ...
                                                                theta, Cn, Cp, an, ap, Ecritn, ... 
                                                                Ecritp)
  
  denomsrh   = tn .* (p + theta) + tp .* (n + theta);
  factauger  = Cn .* n + Cp .* p;
  fact       = (1 ./ denomsrh + factauger);

  Rn = p .* fact;
  Rp = n .* fact;

  Gn = theta .^ 2 .* fact;
  Gp = Gn;

  %II = an * exp(-Ecritn./abs(E)) .* abs (Jn) + ap * exp(-Ecritp./abs(E)) .* abs (Jp);

  II = 0;
end