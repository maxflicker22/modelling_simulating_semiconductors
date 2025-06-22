function u = mobility_model (x, Na, Nd, Nref, E, u0, umin, vsat, beta)

  Neff = Na + Nd;
  Neff = (Neff(1:end-1) + Neff(2:end)) / 2;
  
  ubar = umin + (u0 - umin) ./ (1 + (Neff ./ Nref) .^ beta);
  u    = 2 * ubar ./ (1 + sqrt (1 + 4 * (ubar .* abs (E) ./ vsat) .^ 2));

end