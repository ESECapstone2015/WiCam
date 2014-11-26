classdef f
    properties
    end
    
    methods(Static)
    
        
        
% Function-specific arguments
%
% rect: n_periods, phase, duty_cycle
% tri:  n_periods, phase, duty_cycle
% sin:  n_periods, phase
% rand: <none>


        function [X] = SigGen(type, n_samples, span, varargin)
            X = zeros(n_samples,1);
            
            min_value = span(1);
            max_value = span(2);               
            
            has_periods = 0;
            has_phase = 0;
            has_duty = 0;            
            
            if strcmp(type, 'rect') 
                has_periods = 1;
                has_phase = 1;
                has_duty = 1;
            elseif strcmp(type, 'tri') 
                has_periods = 1;
                has_phase = 1;
                has_duty = 1;
            elseif strcmp(type, 'sin') 
                has_periods = 1;
                has_phase = 1;
            end

            if has_periods == 1        
                n_periods = cell2mat(varargin(1));
            else
                n_periods = 1;
            end
            if has_phase == 1   
                phase = cell2mat(varargin(2));
            else
                phase = 0;
            end;
            if has_duty == 1
                duty_cycle = cell2mat(varargin(3));
            else
                duty_cycle = 0.5;
            end

            n_samples_per_period = floor(n_samples / n_periods);        
            period = zeros(n_samples_per_period+1,1);

            if strcmp(type, 'rect') 
                mid = round(duty_cycle*n_samples_per_period);
                period(1:mid) = max_value;
                period(mid+1:n_samples_per_period) = min_value;
            elseif strcmp(type, 'tri') 
                mid = round(duty_cycle*n_samples_per_period);
                period(1:mid) = linspace(min_value, max_value, mid);
                period(mid:n_samples_per_period+1) = linspace(max_value, min_value, (2+n_samples_per_period)-mid);
            elseif strcmp(type, 'sin') 
                ampl = (max_value - min_value) / 2;
                offset = min_value + ampl;
                period(1:n_samples_per_period+1) = offset + ampl * sin(linspace(0, 2*pi, n_samples_per_period+1));
            elseif strcmp(type, 'rand') 
                ampl = (max_value - min_value);
                period(1:n_samples_per_period) = ampl * rand(n_samples_per_period,1);
            end
            phase_offset = floor(phase * n_samples_per_period);
            phase_offset = mod(phase_offset, n_samples_per_period);
            period = [ period(phase_offset+1:n_samples_per_period); period(1:phase_offset)];
            for i = 1 : n_samples
               X(i) = period( 1 + mod((i-1), n_samples_per_period) ); 
            end

        end


        
        function [bursts] = FindBursts(af)
            
            % Crude large-signal peak detection
            win = abs(sign(af));

            % Peak detector low-pass filter to obtain signal burst windows
            filterSize = 500;
            filterKernel = ones(1, filterSize) / filterSize;
            win = filter(filterKernel, 1, win);
            win = sign(win);

            % Burst window differentiation to identify window start/stop
            win = diff(win);
            [start,c,v] = find(win > 0);
            [stop,c,v] = find(win < 0);
            nBursts = size(start);
            nBursts = nBursts(1);
            ids = [linspace(1, nBursts, nBursts)].';
            bursts = [ [ids], [start], [stop], [stop-start] ]
        end
            
            
        function [fcn] = LinCombo(varargin)
            a = cell2mat(varargin(1 : 2 : end));
            f = cell2mat(varargin(2 : 2 : end));
            npairs = nargin / 2;
            nsmps = size(f) / npairs;
            f = reshape(f, nsmps(2), npairs);
            for i=1:npairs
                f(:,i) = a(i) * f(:,i);
            end   
            f = transpose(f);
            f = sum(f);            
            fcn = f;
        end
        
        function [t] = fromFFT(mag)
            nsamples = numel(mag)-1;
            mag = 10.^(mag/10);
            nPts = 2^nextpow2(nsamples);
            t = ifft(mag, nPts, 'symmetric');
        end
        
        
        function [frq, Y] = FFT2(fcn, smprate)
            nsamples = numel(fcn)-1;
            nPts = 2^nextpow2(nsamples);
            frq = (smprate/2)*linspace(0,1,nPts/2+1);
            %frq = frq.';
            Y = fft(fcn,nPts)/nsamples;
            Y = Y(1:nPts/2+1) + 1e-20;  
            Y = 10*log10(Y);
        end
        
        function [frq, mag] = toFFT(fcn, smprate)
            nsamples = numel(fcn)-1
            nPts = 2^nextpow2(nsamples);
            frq = (smprate/2)*linspace(0,1,nPts/2+1);
            frq = frq.';
            mg = fft(fcn,nPts)/nsamples;
            %mag = 2*abs(mag(1:nPts/2+1));  
            mag = mg(1:nPts/2+1);  
            mag = 10*log(mag);
        end
        
        function [] = plotFFT(fcn, smprate, name)
            [frq, mag] = f.toFFT(fcn, smprate);
            plot(frq.',mag);
            title(['FFT: ', name]);
            xlabel('Frequency');
            ylabel('Magnitude');
        end
        
        function [] = plotFFT2(frq, mag, name)
            plot(frq.',-mag);
            title(['FFT: ', name]);
            xlabel('Frequency');
            ylabel('Magnitude');
        end
        
        function[] = plotSignal(s, name)
            plot(s);
            title(['Signal: ', name]);
            xlabel('Time');
            ylabel('Amplitude');
        end
        
        function [fcn] = Sinusoid(ampl, freq, phase, smprate, nsamples)
            t = 0 : 1/smprate : nsamples/smprate;
            fcn = ampl*sin(2*pi*freq*t + phase);
        end
        
        function [w] = RectWindow(pre, on, post)
            w = [zeros(1,pre), ones(1,on+1), zeros(1,post)];
        end
        
        function [w] = TriWindow(pre, on, post)
            w = [zeros(1,pre), linspace(0, 1, on/2+1), linspace(1, 0, on/2), zeros(1,post)];            
        end    
        
        function [w] = HanningWindow(pre, on, post)
            x = (1:on+1);
            w_on = 0.5 * (1 - cos(2*pi*x/on));
            w = [zeros(1,pre), w_on, zeros(1,post)];    
        end 
        
        function [w] = HammingWindow(pre, on, post)
            x = (1:on+1);
            w_on = 0.53836 - 0.46164 * cos(2*pi*x/on);
            w = [zeros(1,pre), w_on, zeros(1,post)];           
        end 
        
    end    
end

