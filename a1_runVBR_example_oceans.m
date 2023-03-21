%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example script that calculates Vs and Attenuation for a half-space
% cooling oceanic model. It loops through all available anelastic methods
% included in VBR. It also provides an example for specifying a desired
% melt distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% put VBR in the path %%
clear
setup_vbr_paths
vbr_init

% List which anelastic methods to use
methods = {'eburgers_psp','andrade_psp','xfit_mxw','xfit_premelt'};
  
t_age = 50; % [Myrs] seafloor age

dg_mm = 5; % [mm] grain size

phi = 0.01; % melt fraction
zkm_phi_range = [75 150]; % depther over which melt is distributed

% Volatile content. Only affects xfit_premelt by shifting the solidus
H2O_ppm = 0; %;
CO2_ppm = 0; %;

is_correct_Gu_0_ol = 1; % apply correction to the default reference modulus of Isaak? (JBR)

clrs = lines(length(methods));
for imth = 1:length(methods)
    %% Set Thermodynamic State Using Half Space Cooling Model %%
    %%   analytical solution:
    %%     T(z,t)=T_surf + (T_asth - T_surf) * erf(z / (2sqrt(Kappa * t)))
    %%   variables defined below
    
    %Calculate HSC or PLATE cooling geotherm
      % HF settings
      HF.modeltype = 'HSC'; % 'plate' or 'HSC'
      HF.t_Myr = t_age;
      HF.z_plate_km = 100; % plate thickness [km]
      HF.Tp_C = 1380; %[1380 1380]; %[1350 1350]; % Mantle potential temperature [C]
    
    for HFi_t = 1:numel(HF.t_Myr)
        HF.z_m(:,HFi_t) = [(5000:2000:197000),(200000:5000:400000)]; 
        HF.spr_rate_cmyr = 5; % cm/yr
        if strcmpi(HF.modeltype,'hsc')
            [ HF.T_K(:,HFi_t),HF.P_GPa(:,HFi_t),HF.rho_kgm3(:,HFi_t) ] = calc_HSC( HF.Tp_C(HFi_t)+273,HF.t_Myr(HFi_t),HF.spr_rate_cmyr,HF.z_m(:,HFi_t) );
        elseif strcmpi(HF.modeltype,'plate')
            [ HF.T_K(:,HFi_t),HF.P_GPa(:,HFi_t),HF.rho_kgm3(:,HFi_t) ] = calc_platecooling( HF.Tp_C(HFi_t)+273,HF.t_Myr(HFi_t),HF.z_plate_km,HF.spr_rate_cmyr,HF.z_m(:,HFi_t) );
        end
    end
    HF.T_C = HF.T_K - 273;
    HF.z_km = HF.z_m/1000;

    %% Load and set VBR parameters %%
    VBR.in.elastic.methods_list={'anharmonic','anh_poro'};
    VBR.in.viscous.methods_list={'HK2003'};

    method = methods{imth};
    VBR.in.anelastic.methods_list={method}; %{'andrade_psp';'xfit_mxw'};
    VBR.in.anelastic.(method)=Params_Anelastic(method);

    VBR.in.SV.f = logspace(-2.2,-1.3,10); % VBR sweep default

    % copy halfspace model into VBR state variables, adjust units as needed
    VBR.in.SV.T_K = HF.T_K; % set HF temperature, convert to K
    VBR.in.SV.P_GPa = HF.P_GPa; % pressure [GPa]

    % set the other state variables as matrices of same size
    sz=size(HF.T_K);
    VBR.in.SV.rho = 3300 * ones(sz); % density [kg m^-3]
    VBR.in.SV.sig_MPa = 10 * ones(sz); % differential stress [MPa]
    VBR.in.SV.phi = 0.0 * ones(sz); % melt fraction
    Izmelt = HF.z_km(:,1)>=zkm_phi_range(1) & HF.z_km(:,1)<=zkm_phi_range(2);
    VBR.in.SV.phi(Izmelt,1) = phi; 
    
    VBR.in.SV.dg_um = dg_mm * 1e3; % grain size [um]

    % Operations specific to each method
    if strcmp(method,'eburgers_psp')
        VBR.in.anelastic.eburgers_psp.eBurgerFit='bg_only'; % 'bg_only' or 'bg_peak' or 's6585_bg_only'
        GU_lab=VBR.in.anelastic.eburgers_psp.(VBR.in.anelastic.eburgers_psp.eBurgerFit).G_UR; % [GPa] Reference modulus
        PR_lab=VBR.in.anelastic.eburgers_psp.(VBR.in.anelastic.eburgers_psp.eBurgerFit).PR; % [GPa] reference pressure
        TR_lab=VBR.in.anelastic.eburgers_psp.(VBR.in.anelastic.eburgers_psp.eBurgerFit).TR; % [K] reference temperature
    elseif strcmp(method,'andrade_psp')
        GU_lab=VBR.in.anelastic.andrade_psp.G_UR; % [GPa] Reference modulus
        PR_lab=VBR.in.anelastic.andrade_psp.PR; % [GPa] reference pressure
        TR_lab=VBR.in.anelastic.andrade_psp.TR; % [K] reference temperature
    elseif strcmp(method,'xfit_mxw')
        GU_lab=72.66; % [GPa] Reference modulus
        PR_lab=0; % [GPa] reference pressure
        TR_lab=273; % [K] reference temperature
    elseif strcmp(method,'xfit_premelt')
        GU_lab=72.45; % [GPa] Reference modulus
        PR_lab=0; % [GPa] reference pressure
        TR_lab=273; % [K] reference temperature
        
        % Define Solidus
        H2O_solid_ppm=H2O_ppm; % 100 200 ];
        CO2_solid_ppm=CO2_ppm; % 0 0];
        F = [0]; %0.0001; % thermodynamic melt fraction (depletion)
        Cs_H2O=H2O_solid_ppm*1e-4; % ppm to wt%
        Cs_CO2=CO2_solid_ppm*1e-4; % ppm to wt%
        % calculate volatile fractions in melt
        kd_H2O = 1e-2; % equillibrium partition coefficient for H2O
        kd_CO2 = 1e-4; % equillibrium partition coefficient for CO2
        Cf_H2O = Cs_H2O ./ (kd_H2O + F * (1-kd_H2O));
        Cf_CO2 = Cs_CO2 ./ (kd_CO2+ F * (1-kd_CO2));
        % P_Pa = VBR.in.SV.P_GPa(1,1,1)*1e9; % [Pa]
        P_Pa = VBR.in.SV.P_GPa*1e9; % [Pa]
        % [Solidus] = SoLiquidus(P_Pa,Cf_H2O,Cf_CO2,'hirschmann');
        % for ii = 1:length(H2O_solid_ppm)
        [Solidus] = SoLiquidus(P_Pa,Cf_H2O,Cf_CO2,'katz');
        % end
        VBR.in.SV.Tsolidus_K = Solidus.Tsol(:) + 273 ; % [K]
    end

    if is_correct_Gu_0_ol % JBR
        % This corrections accounts for the fact that the laboratory
        % experiments have different reference moduli and P-T conditions.
        % We can "correct" the default moduli of Isaak to agree with the
        % anelasticity laboratory data.
        % Load elastic parameters
        VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); % unrelaxed elasticity
        dGdT=VBR.in.elastic.anharmonic.dG_dT; % [Pa/K]
        dGdP=VBR.in.elastic.anharmonic.dG_dP; % [unitless]
        Tref=VBR.in.elastic.anharmonic.T_K_ref; % [K]
        Pref=VBR.in.elastic.anharmonic.P_Pa_ref/1e9; % [Pa]
        VBR.in.elastic.anharmonic.Gu_0_ol = GU_lab - (TR_lab-Tref) * dGdT/1e9 - (PR_lab-Pref)*dGdP; % olivine reference shear modulus [GPa]
    end

    %% CALL THE VBR CALCULATOR %%
      [VBR] = VBR_spine(VBR) ;

    %% Build figures %%
      % contour T(z,t)
      if imth == 1
        figure(1); clf;
      end
      
      ax(1) = subplot(1,3,1);
      plot(HF.T_C(:,1),HF.z_km,'-','color',clrs(imth,:),'linewidth',2); hold on;
      if isfield(VBR.in.SV,'Tsolidus_K')
        plot(VBR.in.SV.Tsolidus_K-273,HF.z_km,'--k','linewidth',1);
      end
      xlabel('Temperature (C)');
      ylabel('Depth (km)');
      set(gca,'ydir','reverse','fontsize',15,'linewidth',1.5);

      ax(2) = subplot(1,3,2);
      Vsu = VBR.out.elastic.anh_poro.Vsu(:,1); % unrelaxed shear velocity
      Vs = VBR.out.anelastic.(method).Vave(:,1); % relaxed shear velocity
      plot(Vsu/1000,HF.z_km,'--','color',clrs(imth,:),'linewidth',2); hold on;
      plot(Vs/1000,HF.z_km,'-','color',clrs(imth,:),'linewidth',2);
    %   plot(squeeze(VBR.out.anelastic.(VBR.in.anelastic.methods_list{1}).V(:,1,:)),HF.z_km(:,1),'-r');
    %   plot(squeeze(VBR.out.anelastic.(VBR.in.anelastic.methods_list{1}).V(:,2,:)),HF.z_km(:,2),'-b');
      xlabel('Vs (km/s)');
      ylabel('Depth (km)');
      set(gca,'ydir','reverse','fontsize',15,'linewidth',1.5);

      ax(3) = subplot(1,3,3);
      Qinv43 = mean(VBR.out.anelastic.(method).Qinv,2);
      plot(Qinv43,HF.z_km,'-','color',clrs(imth,:),'linewidth',2); hold on;
      xlabel('Q^{-1}');
      ylabel('Depth (km)');
      set(gca,'ydir','reverse','fontsize',15,'linewidth',1.5);

      %%
%       matname = ['./mats_vbr/','VBRout_',VBR.in.anelastic.methods_list{1},'_',VBR.in.anelastic.eburgers_psp.eBurgerFit,'_',HF.modeltype,'_',num2str(HF.t_Myr),'Myr','_T',num2str(HF.Tp_C(1)),'C.mat'];
%       save(matname,'VBR','HF');
end
legend(ax(1),methods,'location','southwest','interpreter','none')

