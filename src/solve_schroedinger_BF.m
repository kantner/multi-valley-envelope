function [S, E] = solve_schroedinger_BF(x, par)
% solve Schrödinger eigenvalue problem in longitudinal direction (growth direction)
% for prescribed epitaxial profile X = X_QW + x
% S the slowly varying envelopes f

  % select sector
    sector = +1; % always solve +k0 sector (-k0 from complex conjugation)
    switch(sector)
      case +1
        idx_SBZ = par.idx_SBZ_p;
        k0      = +par.k0;
      case -1
        idx_SBZ = par.idx_SBZ_m;
        k0      = -par.k0;
      otherwise
        error('invalid option')
    end

    idx   = find(idx_SBZ==1); % convert logical vector into small index set
    N_idx = length(idx);

  % number of eigenvalues to export
    neigs = par.neigs;


    fprintf(1,'construct operator ...')
    tic

  %%%%%%%%%  
  % kinetic energy operator 
    T = par.const.hbar*par.const.hbar/(2*par.mass_l) * diag((par.k(idx) - k0).^2) * 1/par.energy_scale;

  % The operator is hbar^2 /(2 m_l) * K^2 with K \in S_+ 
  % Here we choose K = k-k0 to properly offset everything by k0

  %%%%%%%%%  
  % potential energy operator
    U_x   = potential_modification(x, par);
    U_tot = par.U_QW + par.U_F + U_x;

  % compute Fourier transform of the potential operator (matrix Fourier transform)
    fft_U = 1/par.N * fft(fft(diag(U_tot))')';

  % make Hermitian
    fft_U = 0.5*(fft_U + fft_U');
  
    %figure(2435);clf;hold all;
    %  plot(par.z/par.units.nm, U/par.units.meV,'ko-')

 

    %U_op = zeros(par.N_FBZ/2, par.N_FBZ/2);
    U_op = zeros(N_idx,N_idx);

  % add constant potential
    %U_op  = U_op + eye(N_idx) * 0.25*par.const.hbar*(par.omega_x + par.omega_y);

    for in = 1 : length(par.n_range)
      n = par.n_range(in);


      if abs(par.B(in)) > 0

    % shift the index set to select high-frequency components of the potential folded back to the FBZ    
      idx_shift = find(circshift(idx_SBZ, n * par.N_FBZ) == 1);


      %figure(2436);clf;hold all;
      %imagesc(abs(fft_U(idx,idx_shift))/par.units.meV)
      %colorbar

      %[par.k(idx),par.k(idx_shift)]*par.a0/(2*pi);
      %pause

      %{
      arg_x = n * par.G0(1)*par.l_x/2;
      arg_y = n * par.G0(2)*par.l_y/2;

      arg_x_sq = arg_x*arg_x;
      arg_y_sq = arg_y*arg_y;

      gaussian_x = exp(-arg_x_sq);
      gaussian_y = exp(-arg_y_sq);
      gaussian   = gaussian_x * gaussian_y;
      %}

      %aux = 0.5* par.const.hbar*par.omega_x*(0.5-arg_x_sq) + 0.5* par.const.hbar*par.omega_y*(0.5-arg_y_sq);
      %U_QD_contribution = 0*eye(N_idx) * aux;

      % assemble operator
      U_op = U_op ...
             + par.B(in) * (fft_U(idx, idx_shift))/par.energy_scale;

      % where is the nG0,z ???
      end


    end

  % construct operator (Fourier space)
    H = T + U_op;

  % make Hermitian
    H = 0.5 * (H + H');

    toc

  % solve: eigendecomposition
    fprintf(1,'diagonalize ...')
    tic
    [fft_S,E] = eig(H);
    E = diag(E)*par.energy_scale;
    toc

  % sort (eigenvalues returned by eig are not sorted!)
    [E,ord] = sort(E);

  % apply sorting to eigenvectors
    fft_S = fft_S(:,ord);

  % keep only lowest neigs states
    E = E(1:par.neigs);

  % compute matrix with real space eigenfunctions F as columns
    S = zeros(par.N, par.neigs);

    fprintf(1,'compute real-space wave functions ...')
    tic
    
    %figure(2435);hold all;
    for i = 1 : par.neigs
      phi = zeros(par.N,1);
      phi(idx) = fft_S(:,i);
      S(:,i)   = sqrt(par.N) .* ifft(phi);


      %plot(par.z/par.units.nm, E(i)/par.units.meV + par.N * abs(S(:,i)).^2,'r-')

    end

    toc


  % compute slowly varying envelopes f = exp(-i k0 z) * F
  % k0 here contains the sign corresponding to the sector selected above
    for i = 1 : par.neigs
      S(:,i)   = exp(-1i*k0*par.z) .* S(:,i);
    end





end

