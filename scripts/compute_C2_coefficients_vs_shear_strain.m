function [] = compute_C2_coefficients_vs_shear_strain(par)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % set range
    eps_xy_range = linspace(-1,1,21) * 0.1 * par.units.percent;
  
    nmax = 4;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % sweep
    C2 = zeros(length(par.n_range), length(eps_xy_range));
    B  = zeros(length(par.n_range), length(eps_xy_range));
  
    for ieps = 1 : length(eps_xy_range)
    
    % set new strain tensor  
      par.eps    = par.eps_QW;   
      par.eps(1,2) = eps_xy_range(ieps);
      par.eps(2,1) = par.eps(1,2);
    
    % compute conduction band parameters
      par = compute_conduction_band_parameters(par);
    
    % extract coefficients
      C2(:,ieps) = par.C2;
      B(:,ieps)  = par.B;
    
    end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plot
    figure(100001);clf;hold all;
    cmap = lines(length(par.n_range));
    
    subplot(1,2,1); hold all;
    for in = 1 : length(par.n_range)  
      n = par.n_range(in);
    
      if abs(n)<=nmax
      plot(eps_xy_range/par.units.percent, C2(in,:),'o','MarkerSize',10,'Color',cmap(in,:),'DisplayName',num2str(n))
      end
    
    end
    box on;
    legend()
    xlabel('\varepsilon_{x,y} (%)')
    ylabel('C^{(2)}_n')
    
    subplot(1,2,2); hold all;
    for in = 1 : length(par.n_range)  
      n = par.n_range(in);
    
      if abs(n)<=nmax
      plot(eps_xy_range/par.units.percent, B(in,:),'o','MarkerSize',10,'Color',cmap(in,:),'DisplayName',num2str(n))
      end
    
    end
    box on;
    legend()
    xlabel('\varepsilon_{x,y} (%)')
    ylabel('B_n')
    
  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % least square fit (parabolic)
  
     C2_coeff0 = zeros(length(par.n_range),1);
     C2_coeff1 = zeros(length(par.n_range),1);
     C2_coeff2 = zeros(length(par.n_range),1);
  
     B_coeff0 = zeros(length(par.n_range),1);
     B_coeff1 = zeros(length(par.n_range),1);
     B_coeff2 = zeros(length(par.n_range),1);
  
    for in = 1 : length(par.n_range)  
      n = par.n_range(in);
    
  
      if abs(n)<=nmax
    %%%%%%%
    % fit C2  
      x = eps_xy_range;
      f = C2(in, :);
  
      mean_x  = mean(x);
      mean_x2 = mean(x.^2);
      mean_x3 = mean(x.^3);
      mean_x4 = mean(x.^4);

      mean_f   = mean(f);      
      mean_fx  = mean(f.*x);
      mean_fx2 = mean(f.*x.*x);
  
      mat = [1 mean_x mean_x2; mean_x, mean_x2, mean_x3; mean_x2, mean_x3, mean_x4];
      rhs = [mean_f; mean_fx; mean_fx2];
  
      coeffs = mat\rhs;
  
      C2_coeff0(in) = coeffs(1);
      C2_coeff1(in) = coeffs(2);
      C2_coeff2(in) = coeffs(3);
  
      
    %%%%%%%  
    % fit B  
      x = eps_xy_range;
      f = B(in, :);
  
      mean_x  = mean(x);
      mean_x2 = mean(x.^2);
      mean_x3 = mean(x.^3);
      mean_x4 = mean(x.^4);

      mean_f   = mean(f);      
      mean_fx  = mean(f.*x);
      mean_fx2 = mean(f.*x.^2);
  
      mat = [1 mean_x mean_x2; mean_x, mean_x2, mean_x3; mean_x2, mean_x3, mean_x4];
      rhs = [mean_f; mean_fx; mean_fx2];
  
      coeffs = mat\rhs;
  
      B_coeff0(in) = coeffs(1);
      B_coeff1(in) = coeffs(2);
      B_coeff2(in) = coeffs(3);

      end
    
    end
  
    disp('coefficients for C^(2)_n')
    disp([par.n_range, C2_coeff0, C2_coeff1, C2_coeff2])
    disp('coefficients for B_n')
    disp([par.n_range, B_coeff0, B_coeff1, B_coeff2])
  
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot fit
    
    subplot(1,2,1); hold all;
    for in = 1 : length(par.n_range)  
      n = par.n_range(in);
    
      if abs(n)<=nmax
      plot(eps_xy_range/par.units.percent, C2_coeff0(in) + C2_coeff1(in) * eps_xy_range + C2_coeff2(in) * eps_xy_range.^2,'-','Linewidth',2,'Color',cmap(in,:),'DisplayName',num2str(n))
      end
    
    end
    box on;
    legend()
    xlabel('\varepsilon_{x,y} (%)')
    ylabel('C^{(2)}_n')
    
    subplot(1,2,2); hold all;
    for in = 1 : length(par.n_range)  
      n = par.n_range(in);
    
      if abs(n)<=nmax
      plot(eps_xy_range/par.units.percent, B_coeff0(in) + B_coeff1(in) * eps_xy_range + B_coeff2(in) * eps_xy_range.^2,'-','Linewidth',2,'Color',cmap(in,:),'DisplayName',num2str(n))
      end
    
    end
    box on;
    legend()
    xlabel('\varepsilon_{x,y} (%)')
    ylabel('B_n')


end