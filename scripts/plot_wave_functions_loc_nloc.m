function [] = plot_wave_functions_loc_nloc(out_loc, out_nloc, x, par)



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
      
    % local + nonlocal ground state
      psi_loc = out_loc.S(:,out_loc.idx_gnd);
      psi_nloc = out_nloc.S(:,out_nloc.idx_gnd);

    % compute filtered slowly varying envelope
      psi_filt = exp(-1i*par.k0*par.z)  .* ifft(par.idx_SBZ_p .* fft(exp(1i*par.k0*par.z) .* psi_loc)); 

    % renormalize
      psi_filt = psi_filt/sqrt(psi_filt'*psi_filt);
      
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