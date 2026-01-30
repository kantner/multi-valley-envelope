function [] = Fig_wave_function(par)
  % solve eigenvalue problem and plot wave function
  

  % discretization method (kinetic energy operator, only for LK-type)
    par.discretization = 2; % 1 = finite difference | 2 = spectral

  % grid
    par.N_BZ  = 2^4; % number of Brillouin zones (determines spatial resolution dz)
    par.N_FBZ = 2^10; % number of points in 1st Brillouin zone (determines domain size L and how well k0 is matched on the grid)
    par.N     = par.N_BZ * par.N_FBZ;

    par.L  = 2*pi*par.N_FBZ/par.G0(3);
    par.dz = par.L/par.N; % = 2*pi/(G0(3) * par.N_BZ)
    
  % make position and k-space grid  
    par.z = [0:par.N-1]'*par.dz - par.L + 3*par.h_QW;
    par.k = 2*pi/par.L*[0:par.N/2-1, -par.N/2:-1]';
    
  % SBZ plane wave indicators: omit the origin (+k0) and the lower boundary (-k0) to ensure F_+^* = F_
    par.idx_SBZ_p = 0 < par.k & par.k < par.G0(3)/2;
    par.idx_SBZ_m = -par.G0(3)/2 < par.k & par.k< 0;
    
  % substitute k0 with closest match on grid
    par.k0_backup = par.k0;
    par.k0_target = par.k0;
    idx_k0 = find(abs(par.k-par.k0_target) == min(abs(par.k-par.k0_target)));
    par.k0 = par.k(idx_k0);
    par.k1 = 2*pi/par.a0 * (1-par.eps(3,3)) - par.k0;

    fprintf(1,'original    k0 = %.5f 2*pi/a0\n',par.k0_backup*par.a0/(2*pi))
    fprintf(1,'substitute  k0 = %.5f 2*pi/a0\n',par.k0*par.a0/(2*pi))
    fprintf(1,'difference dk0 = %.5f 2*pi/a0\n',(par.k0_backup-par.k0)*par.a0/(2*pi))

  % neigs
    par.neigs = 20;
   
  % interface width
    par.sigma_l = 0.5 * par.units.nm;
    par.sigma_u = 0.5 * par.units.nm;

  % set QW indicator function  
    par.h_QW         = 72*par.units.ML;
    par.z0           = -par.h_QW/2; % set mid point
    par.QW_indicator = 0.5 * tanh((par.z-par.z0+0.5*par.h_QW)/par.sigma_l) + 0.5*tanh((-par.z+par.z0+0.5*par.h_QW)/par.sigma_u);  % QW at [a-h/2, a+h/2]

  % QW alloy profile
    par.X_barrier = 0.3;
    par.X_QW = par.X_barrier * (1 - par.QW_indicator);  

  % QW potential from nominal (smoothed) profile
    par.U_QW = par.dEc * par.X_QW;    
      
  % set electric field
    par.F    = 3.0 * 1E6;
    par.U_F  = -par.const.e0 * par.F * par.z; 

  % profile modification      
    xGe = 0.6;
    z_spike = -4.0*par.units.nm;
    sigma_spike = 2*par.units.ML;
    x  = xGe * exp(-0.5*((par.z-z_spike)/sigma_spike).^2);

    xww = 0.1;
    q   = 2*par.k1;
    x   = xww * 0.5 * (1+ cos(q*par.z)) .* par.QW_indicator;

  %%%%%%%%%%%%%%  
  % compute valley splitting: local
    par.envelope_model = 1; % 1 = LK-type | 2 = BF-type
    out_loc = compute_valley_splitting(x, 0, par);
    psi_loc = out_loc.S(:,out_loc.idx_gnd);    

  %%%%%%%%%%%%%%  
  % compute valley splitting: non-local
    par.envelope_model = 2; % 1 = LK-type | 2 = BF-type
    out_nloc = compute_valley_splitting(x, 0, par);
    psi_nloc = out_nloc.S(:,out_nloc.idx_gnd);
    
  %%%%%%%%%%%%%%  
  % compute valley splitting: filtered local

  % compute filtered slowly varying envelope
    psi_filt = exp(-1i*par.k0*par.z)  .* ifft(par.idx_SBZ_p .* fft(exp(1i*par.k0*par.z) .* psi_loc)); 

  % renormalize
    psi_filt = psi_filt/sqrt(psi_filt'*psi_filt);

%%
  %%%%%%%%  
  % plot wave functions
    fig = figure(1); clf; hold all;

    subplot(2,2,1); hold all;

    plot_func = @(x) par.N/4 * abs(x).^2;    
    
    overlap_threshold = 0.6;

    overlap_loc  = zeros(par.neigs,1) ;
    overlap_nloc = zeros(par.neigs,1) ;

    

    plot(par.z/par.units.nm, (par.U_F + par.U_QW + par.dEc*x)/par.units.meV, 'k-','LineWidth',2,'DisplayName','potential U')
    for i = 1 : par.neigs
      out = out_loc;

      overlap_loc(i) = sum(abs(out.S(:,i) .* par.QW_indicator))/sqrt(par.QW_indicator'*par.QW_indicator);
      if overlap_loc(i) > overlap_threshold
      plot(par.z/par.units.nm, out.E(i)/par.units.meV + plot_func(out.S(:,i)), '-','Color',[1 1 1]*0.9,'LineWidth',2,'DisplayName','local')
      end

      out = out_nloc;  
      overlap_nloc(i) = sum(abs(out.S(:,i) .* par.QW_indicator))/sqrt(par.QW_indicator'*par.QW_indicator);
      if overlap_nloc(i) > overlap_threshold
      plot(par.z/par.units.nm, out.E(i)/par.units.meV + plot_func(out.S(:,i)), '--','Color',[1 1 1]*0.7,'LineWidth',2,'DisplayName','non-local')
      end

    end
    
  % highlight ground state  
    out = out_loc;
    plot(par.z/par.units.nm, out.E(out.idx_gnd)/par.units.meV + plot_func(psi_loc), '-','Color',[0 0 1]*0.9,'LineWidth',2,'DisplayName','local')

    plot(par.z/par.units.nm, out.E(out.idx_gnd)/par.units.meV + plot_func(psi_filt), ':','Color',[0 1 0]*0.9,'LineWidth',2,'DisplayName','local filtered')
          
    out = out_nloc;    
    plot(par.z/par.units.nm, out.E(out.idx_gnd)/par.units.meV + plot_func(psi_nloc), '--','Color',[1 0 0]*0.9,'LineWidth',2,'DisplayName','non-local')
  
    title('solid = loc | dashed = nloc')
    xlabel('z (nm)')
    ylabel('energy (ueV)')
    box on
    legend()
    xlim([-16, 25])
    ylim([-10, 205])


  %%%%%%%%
  % zoom / inset
    subplot(2,2,2); hold all;

    plot(par.z/par.units.nm, (par.U_F + par.U_QW + par.dEc*x)/par.units.meV, 'k-','LineWidth',2,'DisplayName','potential U')
    for i = 1 : par.neigs

      out = out_loc;
      if overlap_loc(i) > overlap_threshold
      plot(par.z/par.units.nm, out.E(i)/par.units.meV + plot_func(out.S(:,i)), '-','Color',[1 1 1]*0.9,'LineWidth',2,'DisplayName','local')
      end

      out = out_nloc;  
      if overlap_nloc(i) > overlap_threshold
      plot(par.z/par.units.nm, out.E(i)/par.units.meV + plot_func(out.S(:,i)), '--','Color',[1 1 1]*0.7,'LineWidth',2,'DisplayName','non-local')
      end

    end
    
  % highlight ground state  
    out = out_loc;
    plot(par.z/par.units.nm, out.E(out.idx_gnd)/par.units.meV + plot_func(psi_loc), '-','Color',[0 0 1]*0.9,'LineWidth',2,'DisplayName','local')

    plot(par.z/par.units.nm, out.E(out.idx_gnd)/par.units.meV + plot_func(psi_filt), ':','Color',[0 1 0]*0.9,'LineWidth',2,'DisplayName','local filtered')
          
    out = out_nloc;    
    plot(par.z/par.units.nm, out.E(out.idx_gnd)/par.units.meV + plot_func(psi_nloc), '--','Color',[1 0 0]*0.9,'LineWidth',2,'DisplayName','non-local')
  
    title('solid = loc | dashed = nloc')
    xlabel('z (nm)')
    ylabel('energy (ueV)')
    box on
    
    
    XMIN = -7;
    XMAX = 0;

    YMIN = 41;
    YMAX = 60;

    xlim([XMIN, XMAX])
    ylim([YMIN, YMAX])


    subplot(2,2,1); hold all;
    plot([XMIN XMAX XMAX XMIN XMIN], [YMIN YMIN YMAX YMAX YMIN],'k--')

  %%%%%%%%  
  % plot wave functions power spectral density
    subplot(2,2,3); hold all;


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
    xlim([-1.05, 1.05])
    ylim([3E-14, 4])

    drawnow



  %%%%%%%%%%%%%%%%
  % export
    height = 800;
    width  = 1600;
    set(fig,'Position',[100 100 width height]);
    
    exportgraphics(fig, 'Fig_wave_function.pdf','ContentType','vector')

    savefig(fig, 'Fig_wave_function.fig')

end
 