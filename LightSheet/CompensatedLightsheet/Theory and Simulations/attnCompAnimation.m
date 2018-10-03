function attnCompAnimation()

    % load figures from file
    sigmas = [-0.54:0.02:0.54];
    
    
    for sigma_Idx = 1:length(sigmas)
        
        sigma = sigmas(sigma_Idx);
        
        
%         %beam
%         fig = openfig(strcat('C:\Users\Jonathan Nylk\Documents\GitHub\gitJonny_OMGsoftware'...
%             ,'\LightSheet\CompensatedLightsheet\Theory and Simulations\2018-08-10'...
%             ,'\Alpha_7.000000_Abs_0.000000_ExponentialCompensation_sigma_',sprintf('%0.6f',sigma)...
%             ,'\lightSheetProfile.fig'));
%         
%         getAxes = get(gcf,'Children');
% 
%         %attn comp image
%         CData = get(get(getAxes(2),'Children'),'CData');
%         XData = get(get(getAxes(2),'Children'),'XData');
%         YData = get(get(getAxes(2),'Children'),'YData');
% 
%         refImage = get(get(getAxes(1),'Children'),'CData');
% 
%         Xidx = find(XData==0);
%         Yidx = find(YData==1);
%         CData = CData / CData(Yidx,Xidx);
%         CData = CData / 1.6269; % from max compensation
% 
%         CData = uint16(round(CData * (2^16 - 1))); % 16-bit conversion
% 
%         imwrite(CData,strcat('C:\Users\Jonathan Nylk\Documents\GitHub\gitJonny_OMGsoftware'...
%             ,'\LightSheet\CompensatedLightsheet\Theory and Simulations\2018-08-10\Animation\frame_'...
%             ,num2str(1000+sigma_Idx),'.tif'));
%         
%         close(fig)

        fig = openfig(strcat('C:\Users\Jonathan Nylk\Documents\GitHub\gitJonny_OMGsoftware'...
            ,'\LightSheet\CompensatedLightsheet\Theory and Simulations\2018-08-10'...
            ,'\Alpha_7.000000_Abs_0.000000_ExponentialCompensation_sigma_',sprintf('%0.6f',sigma)...
            ,'\pupilAmplitude.fig'));
        
        getAxes = get(gcf,'Children');
% 
        %attn comp pupil
        XData = get(get(getAxes(2),'Children'),'XData');
        YData = get(get(getAxes(2),'Children'),'YData');

        refPupil = get(get(getAxes(1),'Children'),'YData');
        
        Xidx = find(XData==0);
        YData = YData * 0.4768 / refPupil(Xidx);
        fig2 = figure; plot(XData,YData,'w','LineWidth',3);
        ylim([0 1]);
        set(gca,'Color','k')
        axis square;
        frame = getframe(gca);
        
        image_frame = zeros(500,500,3);
        image_frame(71:71+342,71:71+343,:) = frame.cdata;
        imwrite(image_frame,strcat('C:\Users\Jonathan Nylk\Documents\GitHub\gitJonny_OMGsoftware'...
            ,'\LightSheet\CompensatedLightsheet\Theory and Simulations\2018-08-10\Animation\Pupil\frame_'...
            ,num2str(1000+sigma_Idx),'.tif'));
        
        close(fig)
        close(fig2)
       
    end

end