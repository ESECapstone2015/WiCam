<<<<<<< HEAD
classdef gfx
    properties(Constant)       
        
        stem_line_width = 1;
        stem_top_diameter = 7;     
        
        plot_line_width = 1;
        
    end
    
    
    methods(Static)
    
        
        % Plot a discrete-time signal using the stem/lollipop format
        % label: the plot title
        % S: the array of samples to be plotted
        % range [min max]: the range of possible signal values
        % color [R G B]: the color of the plotted signal
        function[] = plotDiscrete(label, S, span, color)
            
            offs = (span(2) + span(1)) / 2;            
            range = (span - offs)*1.2 + offs;
            
            % Draw horizontal grid lines
            GL = zeros(numel(S),1);
            for i = linspace(span(1), span(2), 5)
                GL(:) = i;
                if i == span(1) || i == span(2)
                    plot(GL, 'Color', [.6 .6 .6]);
                else
                    plot(GL, 'Color', [.9 .9 .9]);                    
                end
                hold on;
            end           
            
            label = sprintf('%s (%d samples)', label, numel(S));
            
            stem(S(:,1), 'Color', color, ...
                    'MarkerFaceColor', color, ...
                    'LineWidth', gfx.stem_line_width, ...
                    'MarkerSize', gfx.stem_top_diameter);
             
            axis([0,+inf,range(1),range(2)]);
            %title(['[Discrete] ' label]);
            title(label);
            xlabel('Sample');
            ylabel('Value');       
            hold off;  
        end
        
        
        % Plot a time-domain signal on an unbroken line chart
        % label: the plot title
        % S: an array of [time value] to be plotted (see f.to_time_domain)
        % range [min max]: the range of possible signal values
        % color [R G B]: the color of the plotted signal        
        function[] = plotTD(label, S_t, sample_rate, span, color)
            
            offs = (span(2) + span(1)) / 2;            
            range = (span - offs)*1.2 + offs;
            
            % Draw horizontal grid lines
            GL = S_t;
            for i = linspace(span(1), span(2), 5)
                GL(:,2) = i;
                if i == span(1) || i == span(2)
                    plot(GL(:,1), GL(:,2), 'Color', [.6 .6 .6]);
                else
                    plot(GL(:,1), GL(:,2), 'Color', [.9 .9 .9]);                    
                end
                hold on;
            end     
            
            label = sprintf('%s (%2.2f MSPS)', label, sample_rate/1e6);
            
            plot(S_t(:,1), S_t(:,2), ...
                    'Color',color, ...
                    'LineWidth',gfx.plot_line_width);            
            axis([0,+inf,range(1),range(2)]);
            %title(['[Time Domain] ' label]);
            title(label);
            xlabel('Time (s)');
            ylabel('Amplitude (V)');   
            hold off;  
        end
        
        
        
        function[] = figDiscreteAndTD(label, S, S_t, bits_per_sample, ...
                sample_rate, d_span, td_span, color)
            
            
            persistent figure_number
            
            if isempty(figure_number)
                figure_number = 1;
            end
            
            figure(figure_number);
            figure_number = figure_number + 1;
            
            data_rate = bits_per_sample * sample_rate;
            title_label = sprintf('%s (%2.2f Mbps)',label,data_rate/1e6);
            set(gcf,'numbertitle','off','name',title_label);
            
            subplot(2, 1, 1);
            gfx.plotDiscrete(label, S, d_span, color);
            
            subplot(2, 1, 2);
            gfx.plotTD(label, S_t, sample_rate, td_span, color);

        end
        
        
        
        
        
        
    end
end

=======
classdef gfx
    properties(Constant)       
        
        stem_line_width = 1;
        stem_top_diameter = 7;     
        
        plot_line_width = 1;
        
    end
    
    
    methods(Static)
    
        
        % Plot a discrete-time signal using the stem/lollipop format
        % label: the plot title
        % S: the array of samples to be plotted
        % range [min max]: the range of possible signal values
        % color [R G B]: the color of the plotted signal
        function[] = plotDiscrete(label, S, span, color)
            
            offs = (span(2) + span(1)) / 2;            
            range = (span - offs)*1.2 + offs;
            
            % Draw horizontal grid lines
            GL = zeros(numel(S),1);
            for i = linspace(span(1), span(2), 5)
                GL(:) = i;
                if i == span(1) || i == span(2)
                    plot(GL, 'Color', [.6 .6 .6]);
                else
                    plot(GL, 'Color', [.9 .9 .9]);                    
                end
                hold on;
            end           
            
            label = sprintf('%s (%d samples)', label, numel(S));
            
            stem(S(:,1), 'Color', color, ...
                    'MarkerFaceColor', color, ...
                    'LineWidth', gfx.stem_line_width, ...
                    'MarkerSize', gfx.stem_top_diameter);
             
            axis([0,+inf,range(1),range(2)]);
            %title(['[Discrete] ' label]);
            title(label);
            xlabel('Sample');
            ylabel('Value');       
            hold off;  
        end
        
        
        % Plot a time-domain signal on an unbroken line chart
        % label: the plot title
        % S: an array of [time value] to be plotted (see f.to_time_domain)
        % range [min max]: the range of possible signal values
        % color [R G B]: the color of the plotted signal        
        function[] = plotTD(label, S_t, sample_rate, span, color)
            
            offs = (span(2) + span(1)) / 2;            
            range = (span - offs)*1.2 + offs;
            
            % Draw horizontal grid lines
            GL = S_t;
            for i = linspace(span(1), span(2), 5)
                GL(:,2) = i;
                if i == span(1) || i == span(2)
                    plot(GL(:,1), GL(:,2), 'Color', [.6 .6 .6]);
                else
                    plot(GL(:,1), GL(:,2), 'Color', [.9 .9 .9]);                    
                end
                hold on;
            end     
            
            label = sprintf('%s (%2.2f MSPS)', label, sample_rate/1e6);
            
            plot(S_t(:,1), S_t(:,2), ...
                    'Color',color, ...
                    'LineWidth',gfx.plot_line_width);            
            axis([0,+inf,range(1),range(2)]);
            %title(['[Time Domain] ' label]);
            title(label);
            xlabel('Time (s)');
            ylabel('Amplitude (V)');   
            hold off;  
        end
        
        
        
        function[] = figDiscreteAndTD(label, S, S_t, bits_per_sample, ...
                sample_rate, d_span, td_span, color)
            
            
            persistent figure_number
            
            if isempty(figure_number)
                figure_number = 1;
            end
            
            figure(figure_number);
            figure_number = figure_number + 1;
            
            data_rate = bits_per_sample * sample_rate;
            title_label = sprintf('%s (%2.2f Mbps)',label,data_rate/1e6);
            set(gcf,'numbertitle','off','name',title_label);
            
            subplot(2, 1, 1);
            gfx.plotDiscrete(label, S, d_span, color);
            
            subplot(2, 1, 2);
            gfx.plotTD(label, S_t, sample_rate, td_span, color);

        end
        
        
        
        
        
        
    end
end

>>>>>>> origin/master
