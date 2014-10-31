classdef f
    properties
    end
    
    methods(Static)
    
        
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
        
        function [T] = to_time_domain(D, sample_rate, pts_per_sample, full_scale_digital, full_scale_analog)
            nsamples = size(D, 1)
            T = zeros(nsamples, 2);
            for i = 1 : nsamples*pts_per_sample
                T(i,1) = i/(sample_rate*pts_per_sample);
                T(i,2) = full_scale_analog*D(ceil(i/pts_per_sample))/full_scale_digital;
            end
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
            frq = frq.';
            Y = fft(fcn,nPts)/nsamples;
            Y = Y(1:nPts/2+1);  
        end
        
        function [frq, mag] = toFFT(fcn, smprate)
            nsamples = numel(fcn)-1;
            nPts = 2^nextpow2(nsamples);
            frq = (smprate/2)*linspace(0,1,nPts/2+1);
            frq = frq.';
            mag = fft(fcn,nPts)/nsamples;
            %mag = 2*abs(mag(1:nPts/2+1));  
            mag = mag(1:nPts/2+1);  
            mag = 10*log(mag);
        end
        
        function [] = plotFFT(fcn, smprate, name)
            [frq, mag] = f.toFFT(fcn, smprate);
            plot(frq.',mag);
            title(['FFT: ', name]);
            xlabel('Frequency');
            ylabel('Magnitude');
        end
        
        function [] = plotFFT2(frq, mag, smprate, name)
            plot(frq.',mag);
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

