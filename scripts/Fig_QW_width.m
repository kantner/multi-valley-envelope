function [] = Fig_QW_width(par)
% compute Delta for varying QW depth

  % filename
    filename = 'Fig_QW_width';

  % QW width
    param_range = linspace(5,105,201) * par.units.ML;

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
    par.F    = 3 * 1E6;
    par.U_F  = -par.const.e0 * par.F * par.z;   

  % discretization method (kinetic energy operator, only for LK-type)
    par.discretization = 2; % 1 = finite difference | 2 = spectral

  % neigs
    par.neigs = 10;

    
  %%%%%%%%%%%%%%  
  % sweep
    for i_U = 1 : length(Uoff_range)
      for i_param = 1 : length(param_range)
        
    % interface width
      par.sigma_u = 1 * par.units.ML;
      par.sigma_l = 1 * par.units.ML;

    % set QW indicator function  
      h                = param_range(i_param);
      par.z0           = -h/2; % set mid point      
      par.QW_indicator = 0.5 * tanh((par.z-par.z0+0.5*h)/par.sigma_l) + 0.5*tanh((-par.z+par.z0+0.5*h)/par.sigma_u);  % QW at [a-h/2, a+h/2]

    % QW alloy profile
      par.X_QW = par.X_barrier * (1 - par.QW_indicator);  

    % QW potential from nominal (smoothed) profile
      par.U_QW = par.dEc * par.X_QW;
          
    % apply offset to QW
      par.U_QW = par.dEc * par.X_QW + Uoff_range(i_U);
        
    % profile modification      
      x =  param_range(i_param) * par.QW_indicator;  

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
      figure(1); clf; hold all;

      subplot(1,2,1); hold all;
  
      plot_func = @(x) abs(x).^2;    
      
      plot(par.z/par.units.nm, (par.U_F + par.U_QW + par.dEc*x)/par.units.meV, 'k-','LineWidth',2,'DisplayName','potential U')
      for i = 1 : par.neigs
        out = out_loc;
        plot(par.z/par.units.nm, out.E(i)/par.units.meV + par.N * plot_func(out.S(:,i)), '-','Color',[1 1 1]*0.9,'LineWidth',2,'DisplayName','local')
  
        out = out_nloc;      
        plot(par.z/par.units.nm, out.E(i)/par.units.meV + par.N * plot_func(out.S(:,i)), '--','Color',[1 1 1]*0.7,'LineWidth',2,'DisplayName','non-local')
      end
      
    % highlight ground state  
      out = out_loc;
      plot(par.z/par.units.nm, out.E(out.idx_gnd)/par.units.meV + par.N * plot_func(psi_loc), '-','Color',[0 0 1]*0.9,'LineWidth',2,'DisplayName','local')
  
      plot(par.z/par.units.nm, out.E(out.idx_gnd)/par.units.meV + par.N * plot_func(psi_filt), ':','Color',[0 1 0]*0.9,'LineWidth',2,'DisplayName','local filtered')
            
      out = out_nloc;    
      plot(par.z/par.units.nm, out.E(out.idx_gnd)/par.units.meV + par.N * plot_func(psi_nloc), '--','Color',[1 0 0]*0.9,'LineWidth',2,'DisplayName','non-local')
    
      title('solid = loc | dashed = nloc')
      xlabel('z (nm)')
      ylabel('energy (ueV)')
      box on
      legend()
  
    %%%%%%%%  
    % plot wave functions power spectral density
      subplot(1,2,2); hold all;
  

      xscale = max(abs(fftshift(fft(out.S(:,out.idx_gnd)))).^2);

      out = out_loc;
      plot(fftshift(par.k) * par.a0/(2*pi), abs(fftshift(fft(out.S(:,out.idx_gnd)))).^2 / xscale, '-','Color',[0 0 1]*0.9,'LineWidth',2,'DisplayName','local')
  
      plot(fftshift(par.k) * par.a0/(2*pi), abs(fftshift(fft(psi_filt))).^2 / xscale, '-','Color',[0 1 0]*0.9,'LineWidth',2,'DisplayName','local filtered')

      out = out_nloc;
      plot(fftshift(par.k) * par.a0/(2*pi), abs(fftshift(fft(out.S(:,out.idx_gnd)))).^2 / xscale, '--','Color',[1 0 0]*0.9,'LineWidth',2,'DisplayName','non-local')
  
      xlabel('wave number k (2\pi/a_0)')  
      ylabel('PSD (a.u.)')

      xline((par.G0(3)/2-par.k0)* par.a0/(2*pi),'k--','DisplayName','max(k)')
      xline((0-par.k0)* par.a0/(2*pi),'k--','DisplayName','min(k)')
      
      box on
      set(gca,'YScale','log'  )
      legend()
  
      drawnow

      end
    end

  %%%%%%%%%%%%  
  % save data  
    save([filename,'.mat'],'-v7.3')

  %%%%%%%%%%%%
  % plot
    fig = figure(10006);clf;hold all;

    cmap = brighten(1-lajolla(length(Uoff_range)),0.3);
    cmap = flipud(cmap);

    param_scale = par.units.ML;


    for i_U = 1 : length(Uoff_range)

      color = [0.7 0.7 1]* i_U/length(Uoff_range);

    %%%%%%%%%%%%%%%%%%%%%%
    % valley splitting
      subplot(2,2,1); hold all;

      
      plot(param_range/param_scale, 2*abs(Delta_loc(:,i_U))/par.units.ueV ,'-','Color',cmap(i_U,:),'LineWidth',2,'DisplayName',['local (U_0 = ',num2str(Uoff_range(i_U)/par.units.eV),' eV)'])

      if i_U == length(Uoff_range)
      plot(param_range/param_scale, 2*abs(Delta_nloc(:,i_U))/par.units.ueV,'r-','LineWidth',2,'DisplayName','non-local')
      plot(param_range/param_scale, 2*abs(Delta_filt(:,i_U))/par.units.ueV,'g-.','LineWidth',2,'DisplayName','local filtered')      
      end


      legend('Location','best')
      box on
      xlabel('QW width (ML)')
      ylabel('valley splitting 2 |\Delta| (ueV)')
    

    %%%%%%%%%%%%%%%%%%%%%%
    % absolute error
    subplot(2,2,2); hold all;

      if i_U == 1
      plot(param_range/param_scale, 2*abs(abs(Delta_filt(:,i_U))-abs(Delta_nloc(:,i_U)))/par.units.ueV,'k-','LineWidth',2,'DisplayName','abs error')
      end

      legend()
      box on
      xlabel('QW width (ML)')
      ylabel('abs difference |E_{VS}^{filt}- E_{VS}^{nloc}| (ueV)')
      
    %%%%%%%%%%%%%%%%%%%%%%
    % arbitrariness (violation of gauge invariance)  
    subplot(2,2,3); hold all;
    
      if i_U == 1
      plot(param_range/param_scale, 2*abs(R_loc(:,i_U)),'b-','DisplayName','R_{loc}','LineWidth',2)
      end
      
      %plot(param_range/param_scale, 2*abs(R_nloc(:,i_U)),'r-','DisplayName','R_{nloc}','LineWidth',2)
      %plot(param_range/param_scale, 2*abs(R_filt(:,i_U)),'g--','DisplayName','R_{loc}^{filt}','LineWidth',2)

      legend()
      box on
      xlabel('QW width (ML)')
      ylabel('arbitrariness 2|R|')
    

    %%%%%%%%%%%%%%%%%%%%%%
    % relative error
    subplot(2,2,4); hold all;

      if i_U == 1
      plot(param_range/param_scale, abs(abs(Delta_filt(:,i_U))-abs(Delta_nloc(:,i_U)))./abs(Delta_nloc(:,i_U)),'k-','LineWidth',2,'DisplayName','rel error')
      end

      legend()
      box on
      xlabel('QW width (ML)')
      ylabel('relative error |E_{VS}^{filt}- E_{VS}^{nloc}|/E_{VS}^{nloc}')      
      set(gca,'YScale','log')
    end

    %%%%%%%%%%%%%%
    % export
    height = 800;
    width  = 1600;
    set(fig,'Position',[100 100 width height]);   
    exportgraphics(fig, [filename,'.pdf'],'ContentType','vector')

    savefig(fig, [filename,'.fig'])

end
 