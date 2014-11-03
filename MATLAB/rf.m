classdef rf
    properties(Constant)
         
    end
    
    methods(Static)
        
        
        function [n_channels, bits_per_symbol] = get_modulation_format(type)
            
            n_channels = 0;
            bits_per_symbol = 0;
            
            if strncmp('256qam', type, 6)
                n_channels = 2;
                bits_per_symbol = 4;
            elseif strncmp('64qam', type, 5)
                n_channels = 2;
                bits_per_symbol = 3;
            elseif strncmp('16qam', type, 5)
                n_channels = 2;
                bits_per_symbol = 2;
            elseif strncmp('qpsk', type, 4)
                n_channels = 2;
                bits_per_symbol = 1;
            elseif strncmp('bpsk', type, 4)
                n_channels = 1;
                bits_per_symbol = 1;
            end
        end
    
        
        
        function [] = compare_modulation_types(min_bitrate, types, bws, cpss)
            
            fprintf(' %10s %10s %10s %10s %10s %10s %10s\n', ...
                'ModType', 'BW(MHz)', 'Cyc/Sym', 'Spec.Eff', ...
                'BR(Mbps)', 'SNR(dB)', 'FOM');
            
            fprintf(' %10s %10s %10s %10s %10s %10s %10s\n', ...
                '-------', '--------', '-------', '--------', ...
                '--------', '-------', '-------');
                                
            for type = types'
                type = type';
                for bw = bws' 
                    for cps = cpss'
                        
                        [n_channels, bits_per_symbol] = rf.get_modulation_format(type);

                        % QAM requires whole integer cycles per sample
                        if n_channels == 2 && cps ~= round(cps)
                            n_channels = 0;
                            bits_per_symbol = 0;
                        end
                            
                        spec_eff = n_channels * bits_per_symbol / cps;                      
                        bitrate = spec_eff * bw / 2;
                        
                        oversmp_ratio = 10 * log10(cps) / n_channels;
                        
                        snr = 10 * bits_per_symbol * log10(2) - oversmp_ratio;
                        
                        % Set between 1.5 and 2.5
                        % Set lower to favour performance and efficiency
                        % Set higher to favour lower SNR requirements
                        fom_prefer_eff = 2.0;
                        
                        fom = 1000 * log10(bitrate) * log10(spec_eff+1)^fom_prefer_eff / (log10(bw/2)^2 * snr);
                        
                        if bitrate > 0 && bitrate >= min_bitrate                      
                            fprintf(' %10s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n',... 
                                type, bw/1e6, cps, spec_eff, bitrate/1e6, snr, fom);
                        end
                    end
                end
            end     
        end
        
    end
end

