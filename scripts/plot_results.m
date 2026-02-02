function [] = plot_results(Delta_loc, Delta_nloc, Delta_filt, R_loc, R_nloc, R_filt, Uoff_range, param_range, param_scale, xaxis_label, par)

    cmap = 1-lajolla(length(Uoff_range));
    %cmap = brighten(cmap,0.3);
    cmap = flipud(cmap);

    %cmap = crest(length(Uoff_range));
    %cmap = brighten(cmap,0.5);


    figure(2); clf; hold all;

    for i_U = 1 : length(Uoff_range)

    %%%%%%%%%%%%%%%%%%%%%%
    % valley splitting
      subplot(2,2,1); hold all;

      
      plot(param_range/param_scale, 2*abs(Delta_loc(:,i_U))/par.units.ueV ,'-','Color',cmap(i_U,:),'LineWidth',2,'DisplayName',['local (U_0 = ',num2str(Uoff_range(i_U)/par.units.eV),' eV)'])

      if i_U == length(Uoff_range)
      plot(param_range/param_scale, 2*abs(Delta_nloc(:,i_U))/par.units.ueV,'-','Color',[0.8 0.1 0.0],'LineWidth',2,'DisplayName','non-local')
      plot(param_range/param_scale, 2*abs(Delta_filt(:,i_U))/par.units.ueV,'g-.','Color',[0.6 0 0.6],'LineWidth',2,'DisplayName','local filtered')      
      end

      legend('Location','best')
      box on
      xlabel(xaxis_label)
      ylabel('valley splitting 2 |\Delta| (ueV)')    

    %%%%%%%%%%%%%%%%%%%%%%
    % absolute error
    %{
    subplot(2,2,2); hold all;

      if i_U == 1
      plot(param_range/param_scale, 2*abs(abs(Delta_filt(:,i_U))-abs(Delta_nloc(:,i_U)))/par.units.ueV,'k-','LineWidth',2,'DisplayName','abs error')
      end

      legend()
      box on
      xlabel(xaxis_label)
      ylabel('abs difference |E_{VS}^{filt}- E_{VS}^{nloc}| (ueV)')
    %}
      
    %%%%%%%%%%%%%%%%%%%%%%
    % relative error
    subplot(2,2,2); hold all;

      if i_U == 1
      plot(param_range/param_scale, abs(abs(Delta_filt(:,i_U))-abs(Delta_nloc(:,i_U)))./abs(Delta_nloc(:,i_U)),'k-','LineWidth',2,'DisplayName','rel error')
      end

      legend()
      box on
      xlabel(xaxis_label)
      ylabel('relative error |E_{VS}^{filt}- E_{VS}^{nloc}|/E_{VS}^{nloc}')      
      set(gca,'YScale','log')

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
      xlabel(xaxis_label)
      ylabel('arbitrariness 2|R|')
    

    %%%%%%%%%%%%%%%%%%%%%%
    % R in complex plane
    subplot(2,2,4); hold all;

      if i_U == 2
      plot(real(R_loc(:,i_U)), imag(R_loc(:,i_U)),'ko-','LineWidth',2,'DisplayName','R')
      end

      %legend()
      box on
      xlabel('Re(R)')
      ylabel('Im(R)')
      xline(0,'k--')
      yline(0,'k--')
      axis equal
    end


end