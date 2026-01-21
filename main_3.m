clear all
clc

 %% path
    addpath('colormaps/')
    addpath('src/')
    addpath('src/pseudopotential/')
    addpath('scripts/')
    addpath('scripts/silicon_lattice/')
    addpath('scripts/silicon_lattice/arrow3D_pub/')
    addpath('scripts/silicon_lattice/mArrow3')

 %% parameters
    physical_constants
    
  % constants  
    par.const.hbar = reducedPlanckConstant;
    par.const.e0   = elementaryCharge;
    par.const.m0   = electronMass;
    par.const.c0   = vacuumSpeedOfLight;
    par.const.aB   = bohrRadius;
    
  % units  
    par.units.eV   = elementaryCharge;
    par.units.meV  = 1E-3 * par.units.eV;
    par.units.ueV  = 1E-6 * par.units.eV;
    par.units.Ry   = 0.5*hartree;
    par.units.GPa  = 1E9;
    par.units.nm   = 1E-9;   
    par.units.Angstrom = 0.1 * par.units.nm;    
    par.units.percent  = 1E-2;

  % lattice constant
    par.a0_Si = 5.43 * par.units.Angstrom; % Fischetti % Laux (1996)
    par.a0_Ge = 5.65 * par.units.Angstrom; % Fischetti % Laux (1996)
    
    par.a0_SiGe = @(x) par.a0_Si + 0.200326 * par.units.Angstrom * x .* (1-x) + (par.a0_Ge-par.a0_Si)*x.^2; % Rieger % Vogl (1993)
    
    par.a0      = par.a0_Si;
    par.a0_subs = par.a0_SiGe(0.3);

  % add ML as a unit  
    par.units.ML = par.a0/4;
   
  % elasticity 
    par.C11 = 167.5 * par.units.GPa; % Rieger % Vogl (1993)
    par.C12 =  65.0 * par.units.GPa; % Rieger % Vogl (1993)
    par.C44 =  79.6 * par.units.GPa; % Fischetti % Laux (1996)
    
  % pseudopotential
    par.pp.nonlocal = 1; % 0 = off | 1 = on
    par.pp.SOI      = 0; % 0 = off | 1 = on (without SOI, the Ungersboeck pseudopotential model is identical to the model by Rieger & Vogl)

  % Fermi vector and energy
    par.pp.kF = 1.66 * 2*pi/par.a0; % Ungersboeck (2007)
    par.pp.EF = par.const.hbar^2 * par.pp.kF^2 /(2*par.const.m0);

  % local pseudopotential
    par.pp.V_0      = -2/3 * par.pp.EF;       % Rieger % Vogl (1993)
    par.pp.V_sqrt3  = -0.2241 * par.units.Ry; % Rieger % Vogl (1993)
    par.pp.V_sqrt8  = +0.0520 * par.units.Ry; % Rieger % Vogl (1993)
    par.pp.V_sqrt11 = +0.0724 * par.units.Ry; % Rieger % Vogl (1993)

  % compute cubic spline interpolation parameters 
    par.pp.BC_mode = 1; % different BCs for cubic spline interpolation. 1: V'(0)=0 | 2: V''(0) = 0
    [par.pp.w] = compute_cubic_spline_weights(par);

  % plot atomic pseudopotential
    %plot_atomic_pseudopotential(par)
    
  % nonlocal pseudopotential
    par.pp.A0 = 0.03 * par.units.Ry;       % s-well depth
    par.pp.R0 = 1.06 * par.units.Angstrom; % s-well radius

  % spin-orbit interaction
    par.pp.mu   = 0.00023 * par.units.Ry;
    par.pp.zeta = 7.5589  * 1/par.units.Angstrom;

  % set scale for internal computations  
    par.energy_scale = par.units.Ry;

  % cut off energy  
    par.pp.E_cutoff = 12 * par.units.Ry; % this gives 181 plane waves
  
  % strain tensor  
    par.eps_QW = strain_quantum_well(0.3, par);
    par.eps    = par.eps_QW;
   
    par.eps(1,2) = 0.1 * par.units.percent;
    par.eps(2,1) = par.eps(1,2);

  % Kleinman's internal strain parameter
    par.pp.internal_strain = 0.53;
    
  % compute band structure parameters
    [par] = compute_conduction_band_parameters(par);

  % plot band structure
    %plot_bandstructure(par)  

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% heterostructure part

  %%%%%%%%%%%%%%%%%%
  % quantum dot parameters
    par.omega_x = 3.0 * par.units.meV /par.const.hbar;
    par.omega_y = 3.0 * par.units.meV /par.const.hbar;
    par.l_x     = sqrt(par.const.hbar/(par.mass_t*par.omega_x));
    par.l_y     = sqrt(par.const.hbar/(par.mass_t*par.omega_y));

  % Si/Ge conduction band offset
    par.dEc     = 0.5 * par.units.eV;

  %%%%%%%%%%%%%%%%%  
  % quantum well parameters
    par.h_QW      = 75  * par.units.ML; % QW height
    par.X_barrier = 0.3;                % Ge content in barrier
    par.sigma_u   = 0.5 * par.units.nm; % width of upper interface
    par.sigma_l   = 0.5 * par.units.nm; % width of lower interface

  % vertical electric field (in V/m)
    par.F       = 5.0 * 1E6;

  %%%%%%%%%%%%%%%%%  
  % envelope model 
    par.envelope_model = 1; % 1 = LK-type | 2 = BF-type

  %%%%%%%%%%%%%%%%
  % grid
    par.N_BZ  = 2^4; % number of Brillouin zones (determines spatial resolution dz)
    par.N_FBZ = 2^9; % number of points in 1st Brillouin zone (determines domain size L and how well k0 is matched on the grid)
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

     
  %%%%%%%%%%%%%%%%
  % quantum well indicator function
    par.QW_indicator = 0.5 * tanh((par.z+0.5*par.h_QW)/par.sigma_l) + 0.5*tanh((-par.z+0.5*par.h_QW)/par.sigma_u);  % QW at [-h/2 h/2]
    %par.QW_indicator = 0.5 * tanh((par.z+par.h_QW)/par.sigma_l) + 0.5 * tanh((-par.z)/par.sigma_u);                 % QW at [-h 0]
  
    par.z0 = - par.h_QW/2;
    par.QW_indicator = 0.5 * tanh((par.z-par.z0+0.5*par.h_QW)/par.sigma_l) + 0.5*tanh((-par.z+par.z0+0.5*par.h_QW)/par.sigma_u);  % QW at [a-h/2, a+h/2]

  % nominal alloy profile
    par.X_QW = par.X_barrier * (1 - par.QW_indicator);
    
  % QW potential from nominal (smoothed) profile
    par.U_QW = par.dEc * par.X_QW;

  % electric field potential
    par.U_F  = -par.const.e0 * par.F * par.z;    
    
  %%%%%%%%%%%%%%%%
  % numerics
  % discretization method (kinetic energy operator, only for LK-type)
    par.discretization = 2; % 1 = finite difference | 2 = spectral

  % number of eigenvalues
    par.neigs = 6; %  number of eigenvalues   

  % dummy value for artificially negative covariance functions
    par.Gamma_failure = 1E99;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %% make figures for the paper 
  % list of potential offsets
    par.U0_sweep_range = [-1.0 : 0.1 : 1.0] * par.units.eV;

  % interface dependency
  % Fig_interface_width(par)

  % field dependency
  % Fig_electric_field(par)

  % QW depth  
    Fig_QW_depth(par)

  % wiggle well
  % Fig_wiggle_well(par)  
      
  % Ge spike
  % Fig_Ge_spike(par)
   