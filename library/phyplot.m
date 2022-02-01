function hphyplot=phyplot(varargin)
%PHYPLOT takes in the parameters as plot, along with the following:
%
% phyplot(...,'xticks',A),     where A are the locations of the major tick marks on the X axis [VECTOR].
% phyplot(...,'yticks',A),     where A are the locations of the major tick marks on the Y axis [VECTOR].
% phyplot(...,'xwidth',A),     where A is thickness of the X axis [SCALAR].
% phyplot(...,'ywidth',A),     where A is the thickness of the Y axis [SCALAR].
% phyplot(...,'width',A),      where A is the thickness of the X & Y axes [SCALAR].
% phyplot(...,'fontsize',A),   where A is the font size of the axis tickmark labels [SCALAR].
% phyplot(...,'fontname',A),   where A is the font name of the axis tickmark labels [STRING].
% phyplot(...,'fontangle',A),  where A is the font angle of the axis tickmark labels [STRING].
% phyplot(...,'fontweight',A), where A is the font weight of the axis tickmark labels [STRING].
%
% if an output argument is given, phyplot returns handles to the plotted data
%
% h=phyplot(...),              where h are the handles to the plotted data.
%
% Ex:
%   x = 0:15:360;
%   y = 50.*exp(cosd(x-180));
%   subplot(2,2,1);
%   phyplot(x,y,'k');
%   subplot(2,2,2);
%   phyplot(x,y,'k','xticks',0:90:360,'yticks',0:20:140);
%   subplot(2,2,3);
%   phyplot(x,y,'k','xticks',0:90:360,'yticks',0:50:150,'xwidth',4,'ywidth',4);
%   subplot(2,2,4);
%   phyplot(x,y,'k','xticks',0:90:360,'yticks',0:50:150,'width',4,'fontsize',22,'fontname','Times','fontangle','normal','fontweight','bold');
%
% see also phystring.
%
% version 1.0 - romesh.kumbhani@nyu.edu - 2010-10-08
%
%% defaults

xwidth     = 1;
ywidth     = 1;
fontsize   = 10;
fontname   = 'Helvetica';
fontangle  = 'Oblique';         % normal, italic, oblique
fontweight = 'normal';          % normal, bold
h=-1;

%% Parse optional arguments

% find number of optional arguments
nvarargin = size(varargin,2);

if nvarargin < 1
    error('ExpoMatlab:PlottingFunctions:phyplot','You have to give me something to plot.');
else
    params = 'h=plot(varargin{1}';
    i=2;
    while i<=nvarargin
        if ischar(varargin{i})
            switch varargin{i}
                case 'xticks'
                    xticks = varargin{i+1};
                    i=i+2;
                case 'xticklabels'
                    Xticklabels = varargin{i+1};
                    i=i+2;
                case 'xwidth'
                    xwidth = varargin{i+1};
                    i=i+2;
                case 'yticks'
                    yticks = varargin{i+1};
                    i=i+2;
                case 'yticklabels'
                    Yticklabels = varargin{i+1};
                    i=i+2;
                case 'ywidth'
                    ywidth = varargin{i+1};
                    i=i+2;
                case 'width'
                    xwidth = varargin{i+1};
                    ywidth = varargin{i+1};
                    i=i+2;
                case 'fontsize'
                    fontsize = varargin{i+1};
                    i=i+2;
                case 'fontname'
                    fontname = varargin{i+1};
                    i=i+2;
                case 'fontangle'
                    fontangle = varargin{i+1};
                    i=i+2;
                case 'fontweight'
                    fontweight = varargin{i+1};
                    i=i+2;
                otherwise
                    params = sprintf('%s,varargin{%d}',params,i);
                    i=i+1;
            end
        else
            params = sprintf('%s,varargin{%d}',params,i);
            i=i+1;
        end
    end
    params = [params ');'];
end

try
    eval(params);
catch me
    error('ExpoMatlab:PlottingFunctions:phyplot','Invalid Parameters were given. Please check your syntax.');
end

axis square;

if strcmp(get(gca,'NextPlot'),'replace')
    
    
    % alldata = get(gca,'Children');
    % Xdata = get(alldata,'Xdata');
    % Ydata = get(alldata,'Ydata');
    % minimumX = min([Xdata{:}]);
    % maximumX = max([Xdata{:}]);
    % minimumY = min([Ydata{:}]);
    % maximumY = max([Ydata{:}]);
    %axis([minimumX maximumX minimumY maximumY]);
    
    if exist('xticks','var')
        xlim([xticks(1) xticks(end)]);
    end
    
    if exist('yticks','var')
        ylim([yticks(1) yticks(end)]);
    end
    
    % make the axis range 10% larger than "tight"
    axislimits = axis;
    axisranges = [abs(axislimits(2) - axislimits(1)), abs(axislimits(4) - axislimits(3))];
    axis([axislimits(1)-.1*axisranges(1) axislimits(2)+.1*axisranges(1) axislimits(3)-.1*axisranges(2) axislimits(4)+.1*axisranges(2)]);
    axis off;
    
    if exist('xticks','var')
        xticklocations = xticks;
        xticklabels = num2str(xticks');
    else
        xticklocations = get(gca,'XTick');
        xticklabels = get(gca,'XTickLabel');
    end
    
    if exist('Xticklabels','var')
        xticklabels = Xticklabels;
    end
    
    if exist('yticks','var')
        yticklocations = yticks;
        yticklabels = num2str(yticks');
    else
        yticklocations = get(gca,'YTick');
        yticklabels = get(gca,'YTickLabel');
    end
    
    if exist('Yticklabels','var')
        yticklabels = Yticklabels;
    end
    
    %% DRAW X AXIS
    
    if (xwidth > 0)
        line([xticklocations(1) xticklocations(end)],[axislimits(3) axislimits(3)]-.05*axisranges(2),...
            'color','k','linewidth',xwidth);
        for i=1:length(xticklocations)
            line([xticklocations(i) xticklocations(i)],[axislimits(3)-.05*axisranges(2) axislimits(3)-.08*axisranges(2)],...
                'color','k','linewidth',xwidth)
            if (fontsize>0)
                if ~exist('Xticklabels','var')
                    if str2double(xticklabels(i,:)) == round(str2double(xticklabels(i,:)))
                        text(xticklocations(i)*1.005,axislimits(3) - 0.12*axisranges(2),sprintf('%.0f',str2double(xticklabels(i,:))),...
                            'horizontalalignment','center','fontsize',fontsize,'fontname',fontname,'fontangle',fontangle,...
                            'fontweight',fontweight);
                    else
                        text(xticklocations(i)*1.005,axislimits(3) - 0.12*axisranges(2),sprintf('%.1f',str2double(xticklabels(i,:))),...
                            'horizontalalignment','center','fontsize',fontsize,'fontname',fontname,'fontangle',fontangle,...
                            'fontweight',fontweight);
                    end
                else
                    if iscell(xticklabels)
                        text(xticklocations(i)*1.005,axislimits(3) - 0.12*axisranges(2),sprintf('%s',xticklabels{i}),...
                            'horizontalalignment','center','fontsize',fontsize,'fontname',fontname,'fontangle',fontangle,...
                            'fontweight',fontweight);
                    else
                        text(xticklocations(i)*1.005,axislimits(3) - 0.12*axisranges(2),sprintf('%.0f',xticklabels(i)),...
                            'horizontalalignment','center','fontsize',fontsize,'fontname',fontname,'fontangle',fontangle,...
                            'fontweight',fontweight);
                    end
                end
            end
        end
    end
    %% DRAW Y AXIS
    
    if (ywidth > 0)
        line([axislimits(1) axislimits(1)]-.05*axisranges(1),[yticklocations(1) yticklocations(end)],...
            'color','k','linewidth',ywidth);
        for i=1:length(yticklocations)
            line([axislimits(1)-.05*axisranges(1) axislimits(1)-.08*axisranges(1)],[yticklocations(i) yticklocations(i)],...
                'color','k','linewidth',ywidth)
            if (fontsize>0)
                if ~exist('Yticklabels','var')
                    if str2double(yticklabels(i,:)) == round(str2double(yticklabels(i,:)))
                        text(axislimits(1) - 0.11*axisranges(1),yticklocations(i),sprintf('%.0f',str2double(yticklabels(i,:))),...
                            'horizontalalignment','right','fontsize',fontsize,'fontname',fontname,'fontangle',fontangle,...
                            'fontweight',fontweight);
                    else
                        text(axislimits(1) - 0.11*axisranges(1),yticklocations(i),sprintf('%.1f',str2double(yticklabels(i,:))),...
                            'horizontalalignment','right','fontsize',fontsize,'fontname',fontname,'fontangle',fontangle,...
                            'fontweight',fontweight);
                    end
                else
                    text(axislimits(1) - 0.11*axisranges(1),yticklocations(i),sprintf('%s',yticklabels{i}),...
                        'horizontalalignment','right','fontsize',fontsize,'fontname',fontname,'fontangle',fontangle,...
                        'fontweight',fontweight);
                end
                
            end
        end
    end
    
    %% TAKE CARE OF LABELS
    
    Xpos = get(get(gca,'Xlabel'),'Position');
    Ypos = get(get(gca,'Ylabel'),'Position');
    Tpos = get(get(gca,'Title'),'Position');
    
    Xpos(2) = Xpos(2) +.025*axisranges(2);
    Ypos(1) = Ypos(1) -.03*axisranges(1);
    Tpos(2) = Tpos(2) -.025*axisranges(2);
    
    
    set(get(gca,'Xlabel'),'Visible','on','Position',Xpos,'FontAngle',fontangle,'FontName',fontname,'FontSize',fontsize*4/3,'FontWeight',fontweight);
    set(get(gca,'Ylabel'),'Visible','on','Position',Ypos,'FontAngle',fontangle,'FontName',fontname,'FontSize',fontsize*4/3,'FontWeight',fontweight);
    set(get(gca,'Title'), 'Visible','on','Position',Tpos,'FontAngle',fontangle,'FontName',fontname,'FontSize',fontsize*4/3,'FontWeight',fontweight);
    
    axis square;
    
end
%% TAKE CARE OF OUTPUTS

if nargout == 1
    hphyplot = h;
end

