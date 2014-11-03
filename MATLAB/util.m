
classdef util
    properties
    end
    
    methods(Static)
            
        function [BIT_ARRAY] = to_bit_array(WORDS, nbits)
            nwords = numel(WORDS);
            BIT_ARRAY = zeros(nwords, nbits);            
            for i = 1 : nbits
                BIT_ARRAY(:, i) = bitget(WORDS, i);
            end 
        end 
        
        function [WORDS] = to_words(BIT_ARRAY, nbits)
            nwords = size(BIT_ARRAY,1);
            WORDS = zeros(nwords,1);      
            for i = 1 : nbits
                WORDS = bitset(WORDS, i, BIT_ARRAY(:,i));
            end 
        end 
        
        
        
        
%         function [Y] = bp_to_se(X, in_span, out_span)
%             n_samples = numel(X);
%             Y = zeros(n_samples,1);
%             
%             in_ampl = in_span(2) - in_span(1);
%             in_offs = (in_span(2) + in_span(1)) / 2;
%             
%             Y = X 
%             
%             
%             out_ampl = out_span(2) - out_span(1);
%             out_offs = out_span(1);
%             
%             
%         end
        
        
        function [Y] = scale(X, in_span, out_span)
            
            in_ampl = in_span(2) - in_span(1);
            in_offs = in_span(1);
            
            out_ampl = out_span(2) - out_span(1);
            out_offs = out_span(1);
            
            Y = (X - in_offs) * (out_ampl / in_ampl) + out_offs;            
            
        end
        
        
        
        
        % resolution: number of time-domain data points per discrete sample
        function [T] = to_time_domain(D, sample_rate, resolution, in_span, out_span)
            D = util.scale(D, in_span, out_span);
            pts_per_sample = 1 / (sample_rate * resolution);
            nsamples = size(D, 1);
            T = zeros(nsamples, 2);
            for i = 1 : nsamples*pts_per_sample
                T(i,1) = i/(sample_rate*pts_per_sample);
                T(i,2) = D(ceil(i/pts_per_sample));
            end
        end

            
            
        
    end    
end



