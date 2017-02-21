function [values, R2, RMSE]=GsnFit(x,y)
%'parameters' was not part of the original file generated by cftool, nor
%was the output cf_, which is hopefully the list of fitted parameter
%values.

%CREATEFIT    Create plot of datasets and fits
%   CREATEFIT(X,Y)
%   Creates a plot, similar to the plot in the main curve fitting
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with cftool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1


% Data from dataset "Refolding Curve":
%    X = x:
%    Y = y:
%    Unweighted
%
% This function was automatically generated on 08-Apr-2015 15:49:46

% % Set up figure to receive datasets and fits
% f_ = clf;
% figure(f_);
% set(f_,'Units','Pixels','Position',[118 364 688 488]);
% legh_ = []; legt_ = {};   % handles and text for legend
% xlim_ = [Inf -Inf];       % limits of x axis
% ax_ = axes;
% set(ax_,'Units','normalized','OuterPosition',[0 0 1 1]);
% set(ax_,'Box','on');
% axes(ax_); hold on;
% 
% 
% % --- Plot data originally in dataset "Refolding Curve"
% x = x(:);
% y = y(:);
% h_ = line(x,y,'Parent',ax_,'Color',[0.333333 0 0.666667],...
%     'LineStyle','none', 'LineWidth',1,...
%     'Marker','.', 'MarkerSize',12);
% xlim_(1) = min(xlim_(1),min(x));
% xlim_(2) = max(xlim_(2),max(x));
% legh_(end+1) = h_;
% legt_{end+1} = 'Refolding Curve';
% 
% % Nudge axis limits beyond data limits
% if all(isfinite(xlim_))
%     xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
%     set(ax_,'XLim',xlim_)
% else
%     set(ax_, 'XLim',[0.27999660000000004, 5.7200033999999995]);
% end


% --- Create fit
fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[0   0   0],'Upper',[Inf Inf Inf]);
%fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[0   0],'Upper',[Inf Inf]);

ok_ = isfinite(x) & isfinite(y);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data' );
end
st_ = [1 1 1 1];
%st_ = [1 1 1];
%st_ = [1 1];

%initial values of [A sigma x0], in that order.
set(fo_,'Startpoint',st_);
ft_ = fittype('A/sqrt(2*pi*sigma^2)*exp(-(x-x0)^2/(2*sigma^2))+y0',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'A', 'sigma', 'x0', 'y0'});
% ft_ = fittype('A/sqrt(2*pi*sigma^2)*exp(-(x)^2/(2*sigma^2))',...
%     'dependent',{'y'},'independent',{'x'},...
%     'coefficients',{'A', 'sigma'});

% Fit this model using new data
[cf_, g]= fit(x(ok_),y(ok_),ft_,fo_);
values=coeffvalues(cf_);
R2=g.adjrsquare;
RMSE=g.rmse;


%I added g to report the goodness parameters the otherse to extract values from cf_ and g.

% Or use coefficients from the original fit:
%if 0
%    cv_ = { 0.12872913312004888, 0.26342502437724874, 0.50352676951848441, 0.52324577154471485};
%    cf_ = cfit(ft_,cv_{:});
%end
%Lines 74-77 were originally active, not commented

% % Plot this fit
% h_ = plot(cf_,'fit',0.95);
% legend off;  % turn off legend from plot method call
% set(h_(1),'Color',[1 0 0],...
%     'LineStyle','-', 'LineWidth',2,...
%     'Marker','none', 'MarkerSize',6);
% legh_(end+1) = h_(1);
% legt_{end+1} = 'Dudko 2013 Curve';
% 
% % Done plotting data and fits.  Now finish up loose ends.
% hold off;
% leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'};
% h_ = legend(ax_,legh_,legt_,leginfo_{:});  % create legend
% set(h_,'Interpreter','none');
% xlabel(ax_,'');               % remove x label
% ylabel(ax_,'');               % remove y label
