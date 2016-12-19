

z = [-75:0.5:75] * 1e-6; %[metres]

% betas = [0.05,0.1,0.025,0.25,0.125];
% normalisation_z_positions = [-57.5,-23.5,-126,-6.5,-17.5] * 1e-6; % this is beta dependent, must have same dimensions as "betas" [metres]
% betas = [0.25,0.125,0.1,0.05,0.025];
% normalisation_z_positions = [-6.5,-17.5,-23.5,-57.5,-126] * 1e-6; % this is beta dependent, must have same dimensions as "betas" [metres]
% 
% sigmas = [0.0,0.2,0.4,0.6,0.8,1.0]; % first entry must be zero for comparison
% 
% % C_abss = [25] * 1e2; % [metres^-1]
% C_abss = [10,30,50,70,90] * 1e2; % [metres^-1]

betas = [0.05];
normalisation_z_positions = [-57.5] * 1e-6; % this is beta dependent, must have same dimensions as "betas" [metres]
sigmas = [0.0,0.11,0.22]; % first entry must be zero for comparison
C_abss = [65.02] * 1e2; % [metres^-1]

beamProfile_beforeAttenuation = zeros([length(z),5001 ...
    ,length(betas),length(sigmas),length(C_abss)]);
beamProfile_afterAttenuation = zeros(size(beamProfile_beforeAttenuation));
lightSheetProfile_beforeAttenuation = zeros(size(beamProfile_beforeAttenuation));
lightSheetProfile_afterAttenuation = zeros(size(beamProfile_beforeAttenuation));
beam_attenuationProfile = zeros([length(z)...
    ,length(betas),length(sigmas),length(C_abss)]);
lightSheet_attenuationProfile = zeros([length(z)...
    ,length(betas),length(sigmas),length(C_abss)]);

exponentialDecayFunction = @(params,xData) params(1) .* exp(-1 .* params(2) .* xData) + params(3); % exponential y(x) = Aexp(-Bx)+C
exponentialDecayFunctionMinimisation = @(params,xData,yData) sum((yData - exponentialDecayFunction(params,xData)).^2);


beam_params0 = [0,0,0];
lightSheet_params0 = [0,0,0];
beam_attenuationParams = zeros([length(beam_params0)...
    ,length(betas),length(sigmas),length(C_abss)]);
lightSheet_attenuationParams = zeros(size(beam_attenuationParams));
beam_compensationFactor = zeros([length(betas),length(sigmas),length(C_abss)]);
lightSheet_compensationFactor = zeros(size(beam_compensationFactor));

for betaIdx = 1:length(betas)
    for sigmaIdx = 1:length(sigmas)
        for C_absIdx = 1:length(C_abss)
            
            timer = tic;
            
            beta = betas(betaIdx);
            sigma = sigmas(sigmaIdx);
            C_abs = C_abss(C_absIdx);
            normalisation_z_position = max(normalisation_z_positions(betaIdx),z(1)); % keep within simulation z range.
            fprintf('Beta = %f, sigma = %f, C_abs = % f [cm^-1].\n',beta,sigma,C_abs / 100);
            
            % determine pixel coordinate closest to "normalisation_z_position"
            [~,znorm_Idx] = find(z >= normalisation_z_position - ((z(2) - z(1)) / 2),1,'first');
            zmin_fit = znorm_Idx + 1;
            zmax_fit = find(z >= -1 * normalisation_z_position - (3 * (z(2) - z(1)) / 2),1,'first');
            
            % generate compensated Bessel beam pupil function
            disp('Generating compensated Bessel pupil function.');
            [cartesianPupilFunction,cartesianPupilFunction_amplitude,cartesianPupilFunction_phase,X_cart,Y_cart] = CompensatedBesselLightSheetTheory(beta,sigma);
            disp('Generating compensated Bessel pupil function complete.');
            % simulate beam cross-section and light-sheet cross-section from pupil funciton
            disp('Propagating compensated Bessel pupil function.');
            [beam_CrossSection,lightSheet_CrossSection,transverse_Range,longitudinal_Range] = besselPsfSimulation_ZYang(cartesianPupilFunction,z);
            disp('Propagating compensated Bessel pupil function complete.');
            
            % determine exact transverse localtion (pixel) of beam profile peak if not at transverse_Range = 0;
            [~,maxIdx] = max(beam_CrossSection(floor(length(longitudinal_Range)/2) + 1,:));
            
            % store beam and light-sheet profiles before attenuation
            beamProfile_beforeAttenuation(:,:,betaIdx,sigmaIdx,C_absIdx) = beam_CrossSection / beam_CrossSection(znorm_Idx,maxIdx);
            lightSheetProfile_beforeAttenuation(:,:,betaIdx,sigmaIdx,C_absIdx) = lightSheet_CrossSection / lightSheet_CrossSection(znorm_Idx,maxIdx);
            
            % store beam and light-sheet profiles after attenuation
            absorption_decay = repmat(exp(-C_abs * z).',[1 length(transverse_Range)]);
            beam_CrossSection = beam_CrossSection .* absorption_decay;
            lightSheet_CrossSection = lightSheet_CrossSection .* absorption_decay;
            beam_CrossSection = beam_CrossSection / beam_CrossSection(znorm_Idx,maxIdx);
            lightSheet_CrossSection = lightSheet_CrossSection / lightSheet_CrossSection(znorm_Idx,maxIdx);
            beamProfile_afterAttenuation(:,:,betaIdx,sigmaIdx,C_absIdx) = beam_CrossSection;
            lightSheetProfile_afterAttenuation(:,:,betaIdx,sigmaIdx,C_absIdx) = lightSheet_CrossSection;
            
            % plot decay profiles and infer attenuation compensation factor(sigma)
            beam_attenuationProfile(:,betaIdx,sigmaIdx,C_absIdx) = squeeze(beamProfile_afterAttenuation(:,maxIdx,betaIdx,sigmaIdx,C_absIdx)...
                ./ beamProfile_beforeAttenuation(:,maxIdx,betaIdx,1,C_absIdx));
            lightSheet_attenuationProfile(:,betaIdx,sigmaIdx,C_absIdx) = squeeze(lightSheetProfile_afterAttenuation(:,maxIdx,betaIdx,sigmaIdx,C_absIdx)...
                ./ lightSheetProfile_beforeAttenuation(:,maxIdx,betaIdx,1,C_absIdx));
            beam_params0 = [max(beam_attenuationProfile(:,betaIdx,sigmaIdx,C_absIdx)),C_abs,0];
            lightSheet_params0 = [max(lightSheet_attenuationProfile(:,betaIdx,sigmaIdx,C_absIdx)),C_abs,0];
            params = fminsearch(@(params) exponentialDecayFunctionMinimisation(params,longitudinal_Range(zmin_fit:zmax_fit),beam_attenuationProfile((zmin_fit:zmax_fit),betaIdx,sigmaIdx,C_absIdx).'),beam_params0);
            beam_attenuationParams(:,betaIdx,sigmaIdx,C_absIdx) = params;
            params = fminsearch(@(params) exponentialDecayFunctionMinimisation(params,longitudinal_Range(zmin_fit:zmax_fit),lightSheet_attenuationProfile((zmin_fit:zmax_fit),betaIdx,sigmaIdx,C_absIdx).'),lightSheet_params0);
            lightSheet_attenuationParams(:,betaIdx,sigmaIdx,C_absIdx) = params;
            
            beam_compensationFactor(betaIdx,sigmaIdx,C_absIdx) = beam_attenuationParams(2,betaIdx,sigmaIdx,C_absIdx) - C_abs;
            lightSheet_compensationFactor(betaIdx,sigmaIdx,C_absIdx) = lightSheet_attenuationParams(2,betaIdx,sigmaIdx,C_absIdx) - C_abs;
            
            fprintf('Beam compensation is estimated at %f cm^-1.\nLight-sheet compensation is estimated at %f cm^-1.\n',beam_compensationFactor(betaIdx,sigmaIdx,C_absIdx) / 100,lightSheet_compensationFactor(betaIdx,sigmaIdx,C_absIdx) / 100);
            
            % plot beam and light-sheet axial profiles before and after attenuation
            figure('Name',sprintf('Beta = %f, sigma = %f, C_abs = % f [cm^-1].',beta,sigma,C_abs / 100));
            subplot(3,3,1);
            imagesc(X_cart(1,:),Y_cart(:,1),cartesianPupilFunction_amplitude);
            axis image;
            xlabel('u');
            ylabel('v');
            title('pupil function - amplitude');
            subplot(3,3,4);
            imagesc(X_cart(1,:),Y_cart(:,1),cartesianPupilFunction_phase);
            axis image;
            xlabel('u');
            ylabel('v');
            title('pupil function - phase');
            subplot(3,3,2);
            imagesc(longitudinal_Range * 1e6,transverse_Range * 1e6,beamProfile_beforeAttenuation(:,:,betaIdx,sigmaIdx,C_absIdx).');
            axis image;
            xlabel('x [um]');
            ylabel('z [um]');
            ylim([-30 30]);
            title('beam profile before attenuation');
            subplot(3,3,3);
            imagesc(longitudinal_Range * 1e6,transverse_Range * 1e6,lightSheetProfile_beforeAttenuation(:,:,betaIdx,sigmaIdx,C_absIdx).');
            axis image;
            xlabel('x [um]');
            ylabel('z [um]');
            ylim([-30 30]);
            title('light-sheet profile before attenuation');
            subplot(3,3,5);
            imagesc(longitudinal_Range * 1e6,transverse_Range * 1e6,beamProfile_afterAttenuation(:,:,betaIdx,sigmaIdx,C_absIdx).');
            axis image;
            xlabel('x [um]');
            ylabel('z [um]');
            ylim([-30 30]);
            title('beam profile after attenuation');
            subplot(3,3,6);
            imagesc(longitudinal_Range * 1e6,transverse_Range * 1e6,lightSheetProfile_afterAttenuation(:,:,betaIdx,sigmaIdx,C_absIdx).');
            axis image;
            xlabel('x [um]');
            ylabel('z [um]');
            ylim([-30 30]);
            title('light-sheet profile after attenuation');
            subplot(3,3,8);
            plot(longitudinal_Range * 1e6,beamProfile_beforeAttenuation(:,maxIdx,betaIdx,sigmaIdx,C_absIdx)...
                ,longitudinal_Range * 1e6,beamProfile_afterAttenuation(:,maxIdx,betaIdx,sigmaIdx,C_absIdx));
            xlabel('x [um]');
            ylabel('Intensity [a.u.]');
            title('beam axial profile before/after attenuation');
            subplot(3,3,9);
            plot(longitudinal_Range * 1e6,lightSheetProfile_beforeAttenuation(:,maxIdx,betaIdx,sigmaIdx,C_absIdx)...
                ,longitudinal_Range * 1e6,lightSheetProfile_afterAttenuation(:,maxIdx,betaIdx,sigmaIdx,C_absIdx));
            xlabel('x [um]');
            ylabel('Intensity [a.u.]');
            title('light-sheet axial profile before/after attenuation');
            subplot(3,3,7);
            plot(longitudinal_Range(zmin_fit:zmax_fit) * 1e6,beam_attenuationProfile((zmin_fit:zmax_fit),betaIdx,sigmaIdx,C_absIdx)...
                ,longitudinal_Range(zmin_fit:zmax_fit) * 1e6,exponentialDecayFunction(beam_attenuationParams(:,betaIdx,sigmaIdx,C_absIdx),longitudinal_Range(zmin_fit:zmax_fit))...
                ,longitudinal_Range(zmin_fit:zmax_fit) * 1e6,lightSheet_attenuationProfile((zmin_fit:zmax_fit),betaIdx,sigmaIdx,C_absIdx)...
                ,longitudinal_Range(zmin_fit:zmax_fit) * 1e6,exponentialDecayFunction(lightSheet_attenuationParams(:,betaIdx,sigmaIdx,C_absIdx),longitudinal_Range(zmin_fit:zmax_fit)));
            xlabel('x [um]');
            ylabel('intensity [a.u.]');
            title('beam/light-sheet intensity decay');
            drawnow;shg;
            
            timeTaken = toc(timer);
            fprintf('Time taken is %f seconds (%f minutes).\n',timeTaken,timeTaken / 60);
        end
    end
end




