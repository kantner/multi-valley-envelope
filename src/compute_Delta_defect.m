function R = compute_Delta_defect(psi0,par)

  R = 0;
  for in = 1 : length(par.n_range)

    n = par.n_range(in);
    R = R + par.C2(in) * sum(conj(psi0.^2) .* exp(-1i*(2*par.k0 + n * par.G0(3))*par.z ));

  end

end