%%% name:           multiRun_FastLightSheetProfileSimulation_attenComp
%%% author:         Jonathan Nylk
%%% date created:   22/07/2016
%%% description:    This function runs the function
%%%                 "FastLightSheetProfileSimulation_attenuationCompensation"
%%%                 multiple times with various different input parameter
%%%                 sets.
%%%
%%% updates (latest first):
%%%                 
%%%
%%%
%%% END %%%


        zRange = [-50:0.05:50] * 1e-6;  % transverse beam axis [metres]
        xRange = [-175:1:175] * 1e-6; % propagation axis [metres]
        lambda = 532e-9;    % [metres]
        NA = 0.42;
        alpha = 7;
        compensation_type = 'Exponential';
%         compensation_type = 'Linear';
%         compensation_type = 'Peak';
        sigma_exp = 0.7;
        sigma_lin = 0.7;
        sigma_peak = 0.7;
        x_peak = 0;
        amp_peak = 0.5;
        C_abs = 8.42e3;   % [metres^-1]
        outputFolder = 'C:\Users\Jonathan Nylk\Documents\GitHub\gitJonny_OMGsoftware\LightSheet\CompensatedLightsheet\Theory and Simulations\2018-08-10';
        Flag_imageSimulation = 1;
%         deconvolution_lightSheet = 'APriori';
        deconvolution_lightSheet = 'Normal';
        
%         for C_abs = [0, 64.95 * 100] % [metres^-1]
%             for alpha = [7]
%                 for sigma_exp = [0.27, 0.54]
%                     FastLightSheetProfileSimulation_attenuationCompensation(zRange,xRange,lambda,NA,alpha,compensation_type...
%                         ,sigma_exp,sigma_lin,sigma_peak,x_peak,C_abs,outputFolder,Flag_imageSimulation,deconvolution_lightSheet);
%                 end
%             end
%         end

    
    Flag_imageSimulation = 0;
    
    zRange = [-50:0.1:50] * 1e-6;  % transverse beam axis [metres]
    xRange = [-125:1:125] * 1e-6; % propagation axis [metres]
    lambda = 532e-9;    % [metres]
    NA = 0.42;
    alpha = 7;

        for C_abs = [0] % [metres^-1]
            for alpha = [7]
                for sigma_exp = [-0.54:0.02:0.54]
                    sigma_exp
                    FastLightSheetProfileSimulation_attenuationCompensation(zRange,xRange,lambda,NA,alpha,compensation_type...
                        ,sigma_exp,sigma_lin,sigma_peak,x_peak,C_abs,outputFolder,Flag_imageSimulation,deconvolution_lightSheet);
                    close all
                end
            end
        end


%         Flag_imageSimulation = 0;
%         C_abs = 0   % [metres^-1]
%             for alpha = [7]
%                 for sigma_peak = [0.1:0.2:0.9]
%                     for x_peak = [-0.75:0.25:0.75]
%                         for amp_peak = [0.1:0.2:0.9]
%                             FastLightSheetProfileSimulation_attenuationCompensation(zRange,xRange,lambda,NA,alpha,compensation_type...
%                                 ,sigma_exp,sigma_lin,sigma_peak,x_peak,C_abs,outputFolder,Flag_imageSimulation,deconvolution_lightSheet,amp_peak);
%                         end
%                     end
%                 end
%             end
        
