classdef myFigure < handle
    properties
        fig
        axess
        
        colorbars
        legends
        titles = []
        annotations = []
        texts
        
        size
        pos
        hWidth = 1;
        hSetoff = 0;
        vWidth = 1;
        vSetoff = 0;
        
        tickfontsize = 10
        labelfontsize = 12
        legendfontsize = 10
        textfontsize = 12
        titlefontsize = 12
    end
    methods
%         function myFig = myFigure(fig,size,tickfontsize,labelfontsize,legendfontsize,textfontsize,titlefontsize)
        function myFig = myFigure(fig,varargin)
            drawnow
            set(fig,'units','centimeters','PaperUnits','centimeters');
            myFig.fig = fig;
            myFig.update(varargin{:})
            drawnow
        end
        
        function update(self,varargin)
            drawnow
            
            for n = 1:2:length(varargin)
                self.(varargin{n}) = varargin{n+1};
            end
            
            self.resize()
            
            handles = findobj(self.fig,'type','axes');
            filter = false(length(handles),1);
            for i = 1:length(handles) % remove invisble handles
                if strcmp(handles(i).Visible,'off')
                    filter(i) = true;
                end
            end
            handles(filter) = [];
            set(handles,'fontname','times','fontsize',self.tickfontsize,...
                'LabelFontSizeMultiplier',self.labelfontsize/self.tickfontsize,...
                'TickLabelInterpreter','latex')
            if ~isempty(handles)
                set([handles.XLabel],'fontname','times','fontsize',self.labelfontsize,'interpreter','latex')
                set([handles.YLabel],'fontname','times','fontsize',self.labelfontsize,'interpreter','latex')
                set([handles.Title],'fontname','times','fontsize',self.titlefontsize,'interpreter','latex')
            end
            self.axess = handles;

            handles = findobj(self.fig,'type','colorbar');
            set(handles,'fontname','times','fontsize',self.tickfontsize,'TickLabelInterpreter','latex')
            if ~isempty(handles)
                set([handles.Label],'fontname','times','fontsize',self.labelfontsize,'interpreter','latex')
            end
            self.colorbars = handles;
            
            handles = findobj(self.fig,'type','text');
            set(handles,'fontname','times','fontsize',self.textfontsize,'interpreter','latex')
            self.texts = handles;
            
            handles = findobj(self.fig,'type','legend');
            set(handles,'fontname','times','fontsize',self.legendfontsize,'interpreter','latex')
            self.legends = handles;
            
            drawnow
        end
        
        function resize(self)
            if isempty(self.size)
                self.size = self.fig.Position(3:4);
            end
            if isempty(self.pos)
                self.pos = self.fig.Position(1:2);
            end
            pause(0)
            self.fig.Position(1:2) = self.pos;
            self.fig.Position(3:4) = self.size;
            self.fig.PaperSize = self.size;
            self.fig.PaperPosition = [0 0 self.size];
        end
        
        
        
        function title(myFig,titleString,fontSize,alignment)
            figure(myFig.fig);
            handles = myFig.axess;
            
            if nargin<3 || numel(fontSize)==0, fontSize = myFig.titlefontsize; end
            if nargin<4, alignment = []; end
            
            if length(titleString)>1
                for i = 1:length(handles)
                    axes(handles(end-i+1));
                    myFig.titles{end+1} = title(titleString{i},'interpreter','latex','FontSize',fontSize);
                end
            elseif length(titleString)==1 && length(handles)>1
                h = annotation('textbox',[myFig.hSetoff 0.9 myFig.hWidth 0.1],'String',titleString{1,:}, ...
                    'EdgeColor','none','FontSize',fontSize,'interpreter','latex', ...
                    'HorizontalAlignment','center','VerticalAlignment','top');
                if numel(alignment)>0
                    h.Position = [alignment 0.9 0 0.1];
                end
                myFig.titles{end+1} = h;
            else
                myFig.titles{end+1} = title(titleString{1},'interpreter','latex','FontSize',fontSize);
            end
        end
        
        
        
        function print(myFig,title,key)
            figure(myFig.fig);
            savefig(title)
            if nargin < 3
                key = '-dpng';
            end
            
            if strcmp(key,'tikz')
                cleanfigure;
                matlab2tikz([title,'.tikz'],'showInfo',false,'height','\figureheight','width','\figurewidth');
            else
                print(title,key,'-r1000')
            end
        end
        
        
        
        function subPlots(self,arrangement,varargin)
            figure(self.fig);
            
            % default options
            % all measures in cm
            options.x1 = 0; % left space
            options.x2 = 0; % spacing
            options.x3 = 0; % right space
            options.y1 = 0; % bottom space
            options.y2 = 0; % spacing
            options.y3 = 0; % top space
            options.order = 1:prod(arrangement);
            options.aspect_ratio = []; % resizes figure height to accomodate
            
            % input options
            for n = 1:2:length(varargin)
                options.(varargin{n}) = varargin{n+1};
            end
            
            x1 = options.x1;
            x2 = options.x2;
            x3 = options.x3;
            y1 = options.y1;
            y2 = options.y2;
            y3 = options.y3;
            order = flip(options.order);
            
            panel_width = (self.size(1) - x1 - x3 - (arrangement(2)-1)*x2)/arrangement(2);
            if isempty(options.aspect_ratio)
                panel_height = (self.size(2) - y1 - y3 - (arrangement(1)-1)*y2)/arrangement(1);
            else
                panel_height = options.aspect_ratio * panel_width;
                new_figure_height = y1 + y3 + (arrangement(1)-1)*y2 + arrangement(1)*panel_height;
                self.size(2) = new_figure_height;
                self.resize();
            end
            
            count = 1;
            y = self.size(2) - y3 - panel_height;
            for m = 1:arrangement(1)
                x = x1;
                for n = 1:arrangement(2)
                    if order(count)<=length(self.axess)
                        self.axess(order(count)).Units = 'centimeters';
                        self.axess(order(count)).Position = [x y panel_width panel_height];
                    end
                    count = count+1;
                    if count==length(order)+1
                        break
                    end
                    x = x + x2 + panel_width;
                end
                if count==length(order)+1
                    break
                end
                y = y - y2 - panel_height;
            end
        end    
        
        
        function annotation(self,pos,text,fontsize)
            if nargin<4
                fontsize = self.textfontsize;
            end
            width = min([pos(1) abs(1-pos(1))]);
            height = min([pos(2) abs(1-pos(2))]);
            self.annotations{end+1} = annotation('textbox',[pos(1)-width pos(2)-height 2*width 2*height],'HorizontalAlignment','center','VerticalAlignment','middle','String',text,'EdgeColor','none','FontSize',fontsize,'interpreter','latex');
        end
    end
end
