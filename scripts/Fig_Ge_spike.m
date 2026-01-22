function [] = Fig_Ge_spike(par)
% compute Delta for Ge spike

  % filename
    filename = 'Fig_Ge_spike';

  % spike position
    param_range = linspace(-1,0,201)*par.h_QW;

  % Ge amplitude  
    xGe = 0.6;

  % spike width
    sigma_spike = 2 * par.units.ML;

  % potential offsets
    Uoff_range = par.U0_sweep_range;

  % allocate memory  
    Delta_loc  = zeros(length(param_range), length(Uoff_range));
    Delta_nloc = zeros(length(param_range), length(Uoff_range));
    Delta_filt = zeros(length(param_range), length(Uoff_range));

    R_loc  = zeros(length(param_range), length(Uoff_range));
    R_nloc = zeros(length(param_range), length(Uoff_range));
    R_filt = zeros(length(param_range), length(Uoff_range));


  % set electric field
    par.F    = 3.0 * 1E6;
    par.U_F  = -par.const.e0 * par.F * par.z;   

  % discretization method (kinetic energy operator, only for LK-type)
    par.discretization = 2; % 1 = finite difference | 2 = spectral

  % neigs
    par.neigs = 3;

  % set QW indicator function  
    par.z0           = -par.h_QW/2; % set mid point
    par.QW_indicator = 0.5 * tanh((par.z-par.z0+0.5*par.h_QW)/par.sigma_l) + 0.5*tanh((-par.z+par.z0+0.5*par.h_QW)/par.sigma_u);  % QW at [a-h/2, a+h/2]

  % QW alloy profile
    par.X_QW = par.X_barrier * (1 - par.QW_indicator);  

  % QW potential from nominal (smoothed) profile
    par.U_QW = par.dEc * par.X_QW;
    
  %%%%%%%%%%%%%%  
  % sweep
    for i_U = 1 : length(Uoff_range)
      for i_param = 1 : length(param_range)

    % apply offset to QW
      par.U_QW = par.dEc * par.X_QW + Uoff_range(i_U);

    % profile modification      
      x      = xGe * exp(-0.5*((par.z-param_range(i_param))/sigma_spike).^2);

    %%%%%%%%%%%%%%  
    % compute valley splitting: local
      par.envelope_model = 1; % 1 = LK-type | 2 = BF-type
      out_loc = compute_valley_splitting(x, 0, par);
      psi_loc = out_loc.S(:,out_loc.idx_gnd);

      Delta_loc(i_param, i_U)  = out_loc.Delta_det;

    %%%%%%%%%%%%%%  
    % compute valley splitting: non-local
      par.envelope_model = 2; % 1 = LK-type | 2 = BF-type
      out_nloc = compute_valley_splitting(x, 0, par);
      psi_nloc = out_nloc.S(:,out_nloc.idx_gnd);

      Delta_nloc(i_param, i_U) = out_nloc.Delta_det;      
      
    %%%%%%%%%%%%%%  
    % compute valley splitting: filtered local

    % compute filtered slowly varying envelope
      psi_filt = exp(-1i*par.k0*par.z)  .* ifft(par.idx_SBZ_p .* fft(exp(1i*par.k0*par.z) .* psi_loc)); 

    % renormalize
      psi_filt = psi_filt/sqrt(psi_filt'*psi_filt);

    % compute valley coupling parameter
      Delta_filt(i_param, i_U) = compute_Delta_det(psi_filt, x, par);


    %%%%%%%%%%%  
    % compute parameter R (measure for gauge invariance violation)

    % store potentials
      par.U_F_backup  = par.U_F;
      par.U_QW_backup = par.U_QW;      

    % replace by unity potential  
      par.U_F  = ones(size(par.z));
      par.U_QW = zeros(size(par.z));

    % compute  
      R_loc(i_param, i_U)  = compute_Delta_det(psi_loc,0*x,par);
      R_nloc(i_param, i_U) = compute_Delta_det(psi_nloc,0*x,par);
      R_filt(i_param, i_U) = compute_Delta_det(psi_filt,0*x,par);

    % swap back  
      par.U_F  = par.U_F_backup;
      par.U_QW = par.U_QW_backup;


    %%%%%%%%  
    % plot wave functions
      plot_wave_functions_loc_nloc(out_loc, out_nloc, x, par)

      end
    end

  %%%%%%%%%%%%  
  % save data  
    save([filename,'.mat'],'-v7.3')

  %%%%%%%%%%%%
  % plot    
    param_scale = par.units.nm;
    xaxis_label = 'spike position (nm)';
    plot_results(Delta_loc, Delta_nloc, Delta_filt, R_loc, R_nloc, R_filt, Uoff_range, param_range, param_scale, xaxis_label, par)

  %%%%%%%%%%%%%%%%
  % export
    height = 800;
    width  = 1600;
    set(gcf,'Position',[100 100 width height]);   
    exportgraphics(gcf, [filename,'.pdf'],'ContentType','vector')

    savefig(gcf, [filename,'.fig'])

end
 