<<<<<<< HEAD

classdef fpga
    properties
    end
    
    methods(Static)
    
        % X = [ [ABCD] ]
        %
        % X = fold(X)
        %   => [ [AB]; [CD] ]
        %
        % X = fold(X)
        %   => [ [A]; [B]; [C]; [D] ]
        %        
        function [Y] = fold(X)            
            n = size(X,1)*2;
            b = size(X,2)/2;
            Y = zeros(n,b);
            Y(1:2:n-1) = X(:,1:b);
            Y(2:2:n) = X(:,b+1:2*b);
        end
        
        
        function [I, Q] = map_IQ_8bit_256qam(S)
            
            n_source_samples = size(S,1);            
            n_iq_samples = n_source_samples;
            
            I = zeros(n_iq_samples, 4);
            Q = zeros(n_iq_samples, 4);
            
            B = util.to_bit_array(S, 8);
            %B = fold(B);
            
            
            
            %for i = nbits/2 : -1 : 1
            %    I(:,i) = B(:,2*i-1);
            %    Q(:,i) = B(:,2*i-0);
            %end 
            
            I(:,4:-1:1) = B(:,7:-2:1);
            Q(:,4:-1:1) = B(:,8:-2:2);
            
            
            I = util.to_words(I, 4);
            Q = util.to_words(Q, 4);
        end
        
        
%         
%         
%         function [I, Q] = map_IQ(S, n_source_bits, n_iq_bits, n_channels)
%             
%             n_source_samples = size(S,1);
%             
%             data_size = n_source_samples * n_source_bits;
%                         
%             n_iq_samples = data_size / (n_iq_bits * n_channels);
%             
%             
%             I = zeros(n_iq_samples, n_iq_bits);
%             Q = zeros(n_iq_samples, n_iq_bits);
%         
%             
%             BIT_ARRAY = util.to_bit_array(S, n_source_bits);
%                       
%             n_reshaped_source_samples = n_source_samples * n_iq_bits / (n_source_bits * n_channels);
%             
%             X = zeros(n_reshaped_source_samples, n_iq_bits);
%         
%             X = BIT_ARRAY(
%         
            
            
        function [I, Q] = map_IQ(S, n_source_bits, n_iq_bits, n_channels)
            
            n_source_samples = size(S,1);
            
            data_size = n_source_samples * n_source_bits;
                        
            n_iq_samples = data_size / (n_iq_bits * n_channels);
            
            
            I = zeros(n_iq_samples, n_iq_bits);
            Q = zeros(n_iq_samples, n_iq_bits);
            
            BIT_ARRAY = util.to_bit_array(S, n_source_bits);
            
            
            iq_samples_per_source_sample = n_iq_samples / n_source_samples;
            
            
            for s = 1 : n_source_samples
                
                for i = 1 : iq_samples_per_source_sample
                    
                    for b = 1 : n_iq_bits
                        j = i + (s-1)*iq_samples_per_source_sample;
                        n = n_iq_bits * n_channels * (i-1) + b;
                        I(j,b) = BIT_ARRAY(s,n);
                        Q(j,b) = BIT_ARRAY(s,n+(n_channels-1));
                    end
                end
            end
            
            
%             
%             if strcmp(mode, 'interleave')
%                 for i = nbits/2 : -1 : 1
%                     I(:,i) = BIT_ARRAY(:,2*i-1);
%                     Q(:,i) = BIT_ARRAY(:,2*i-0);
%                 end 
%             elseif strcmp(mode, 'split')
%                 for i = nbits/2 : -1 : 1
%                     I(:,i) = BIT_ARRAY(:,i);
%                     Q(:,i) = BIT_ARRAY(:,nbits/2+i);
%                 end 
%             else
%                 fprintf('map_to_IQ: unrecognized mode\n');
%             end         
            
            I = util.to_words(I, n_iq_bits);
            Q = util.to_words(Q, n_iq_bits);
            
        end
        
        
        
        
        function [WORDS] = unmap_IQ(I, Q, nbits, mode)
            
            nwords = size(I,1);            
            I_BITS = util.to_bit_array(I, nbits/2);
            Q_BITS = util.to_bit_array(Q, nbits/2);
            
            WORD_BITS = zeros(nwords, nbits);
            
            if strcmp(mode, 'interleave')
                for i = nbits/2 : -1 : 1
                    WORD_BITS(:,2*i-1) = I_BITS(:,i);
                    WORD_BITS(:,2*i-0) = Q_BITS(:,i);
                end 
            elseif strcmp(mode, 'split')
                for i = nbits/2 : -1 : 1
                    WORD_BITS(:,i) = I_BITS(:,i);
                    WORD_BITS(:,nbits/2+i) = Q_BITS(:,i);
                end 
            else
                fprintf('map_to_IQ: unrecognized mode\n');
            end         
            
            WORDS = util.to_words(WORD_BITS, nbits);
                        
        end
        
        
        
        
        function [ILC] = line_encode(I, span, symbols_per_cycle, always_cross_zero)

            
            n_samples = size(I,1);

            I = util.scale(I, span, [0,1]);
            
            %I = I - span(1);
            %I = I / (span(2) - span(1));
       
            if symbols_per_cycle == 1

                if always_cross_zero == 1
                    ILC = zeros(2*n_samples,1);
                    ILC(1:2:2*n_samples-1) = I;
                    ILC(2:2:2*n_samples) = 1-I;
                else
                    ILC = zeros(2*n_samples,1);
                    ILC(1:2:2*n_samples-1) = I;
                    ILC(2:2:2*n_samples) = 1-I;
                end    

            elseif symbols_per_cycle == 2

                if always_cross_zero == 1
                    ILC = zeros(n_samples,1);
                    ILC(1:2:n_samples-1) = (1+(1+I(1:2:n_samples-1)))/2;
                    ILC(2:2:n_samples) = (1-(1+I(2:2:n_samples)))/2;
                else        
                    ILC = zeros(n_samples,1);
                    ILC(1:2:n_samples-1) = I(1:2:n_samples-1);
                    ILC(2:2:n_samples) = 1-I(2:2:n_samples);
                end
            end
            
            
            ILC = util.scale(ILC, [0,1], span);
            
           % ILC = ILC * (span(2) - span(1)) / 2;
           % ILC = ILC + span(1);
            

        end
    end    
end


=======

classdef fpga
    properties
    end
    
    methods(Static)
    
        % X = [ [ABCD] ]
        %
        % X = fold(X)
        %   => [ [AB]; [CD] ]
        %
        % X = fold(X)
        %   => [ [A]; [B]; [C]; [D] ]
        %        
        function [Y] = fold(X)            
            n = size(X,1)*2;
            b = size(X,2)/2;
            Y = zeros(n,b);
            Y(1:2:n-1) = X(:,1:b);
            Y(2:2:n) = X(:,b+1:2*b);
        end
        
        
        function [I, Q] = map_IQ_8bit_256qam(S)
            
            n_source_samples = size(S,1);            
            n_iq_samples = n_source_samples;
            
            I = zeros(n_iq_samples, 4);
            Q = zeros(n_iq_samples, 4);
            
            B = util.to_bit_array(S, 8);
            %B = fold(B);
            
            
            
            %for i = nbits/2 : -1 : 1
            %    I(:,i) = B(:,2*i-1);
            %    Q(:,i) = B(:,2*i-0);
            %end 
            
            I(:,4:-1:1) = B(:,7:-2:1);
            Q(:,4:-1:1) = B(:,8:-2:2);
            
            
            I = util.to_words(I, 4);
            Q = util.to_words(Q, 4);
        end
        
        
%         
%         
%         function [I, Q] = map_IQ(S, n_source_bits, n_iq_bits, n_channels)
%             
%             n_source_samples = size(S,1);
%             
%             data_size = n_source_samples * n_source_bits;
%                         
%             n_iq_samples = data_size / (n_iq_bits * n_channels);
%             
%             
%             I = zeros(n_iq_samples, n_iq_bits);
%             Q = zeros(n_iq_samples, n_iq_bits);
%         
%             
%             BIT_ARRAY = util.to_bit_array(S, n_source_bits);
%                       
%             n_reshaped_source_samples = n_source_samples * n_iq_bits / (n_source_bits * n_channels);
%             
%             X = zeros(n_reshaped_source_samples, n_iq_bits);
%         
%             X = BIT_ARRAY(
%         
            
            
        function [I, Q] = map_IQ(S, n_source_bits, n_iq_bits, n_channels)
            
            n_source_samples = size(S,1);
            
            data_size = n_source_samples * n_source_bits;
                        
            n_iq_samples = data_size / (n_iq_bits * n_channels);
            
            
            I = zeros(n_iq_samples, n_iq_bits);
            Q = zeros(n_iq_samples, n_iq_bits);
            
            BIT_ARRAY = util.to_bit_array(S, n_source_bits);
            
            
            iq_samples_per_source_sample = n_iq_samples / n_source_samples;
            
            
            for s = 1 : n_source_samples
                
                for i = 1 : iq_samples_per_source_sample
                    
                    for b = 1 : n_iq_bits
                        j = i + (s-1)*iq_samples_per_source_sample;
                        n = n_iq_bits * n_channels * (i-1) + b;
                        I(j,b) = BIT_ARRAY(s,n);
                        Q(j,b) = BIT_ARRAY(s,n+(n_channels-1));
                    end
                end
            end
            
            
%             
%             if strcmp(mode, 'interleave')
%                 for i = nbits/2 : -1 : 1
%                     I(:,i) = BIT_ARRAY(:,2*i-1);
%                     Q(:,i) = BIT_ARRAY(:,2*i-0);
%                 end 
%             elseif strcmp(mode, 'split')
%                 for i = nbits/2 : -1 : 1
%                     I(:,i) = BIT_ARRAY(:,i);
%                     Q(:,i) = BIT_ARRAY(:,nbits/2+i);
%                 end 
%             else
%                 fprintf('map_to_IQ: unrecognized mode\n');
%             end         
            
            I = util.to_words(I, n_iq_bits);
            Q = util.to_words(Q, n_iq_bits);
            
        end
        
        
        
        
        function [WORDS] = unmap_IQ(I, Q, nbits, mode)
            
            nwords = size(I,1);            
            I_BITS = util.to_bit_array(I, nbits/2);
            Q_BITS = util.to_bit_array(Q, nbits/2);
            
            WORD_BITS = zeros(nwords, nbits);
            
            if strcmp(mode, 'interleave')
                for i = nbits/2 : -1 : 1
                    WORD_BITS(:,2*i-1) = I_BITS(:,i);
                    WORD_BITS(:,2*i-0) = Q_BITS(:,i);
                end 
            elseif strcmp(mode, 'split')
                for i = nbits/2 : -1 : 1
                    WORD_BITS(:,i) = I_BITS(:,i);
                    WORD_BITS(:,nbits/2+i) = Q_BITS(:,i);
                end 
            else
                fprintf('map_to_IQ: unrecognized mode\n');
            end         
            
            WORDS = util.to_words(WORD_BITS, nbits);
                        
        end
        
        
        
        
        function [ILC] = line_encode(I, span, symbols_per_cycle, always_cross_zero)

            
            n_samples = size(I,1);

            I = util.scale(I, span, [0,1]);
            
            %I = I - span(1);
            %I = I / (span(2) - span(1));
       
            if symbols_per_cycle == 1

                if always_cross_zero == 1
                    ILC = zeros(2*n_samples,1);
                    ILC(1:2:2*n_samples-1) = I;
                    ILC(2:2:2*n_samples) = 1-I;
                else
                    ILC = zeros(2*n_samples,1);
                    ILC(1:2:2*n_samples-1) = I;
                    ILC(2:2:2*n_samples) = 1-I;
                end    

            elseif symbols_per_cycle == 2

                if always_cross_zero == 1
                    ILC = zeros(n_samples,1);
                    ILC(1:2:n_samples-1) = (1+(1+I(1:2:n_samples-1)))/2;
                    ILC(2:2:n_samples) = (1-(1+I(2:2:n_samples)))/2;
                else        
                    ILC = zeros(n_samples,1);
                    ILC(1:2:n_samples-1) = I(1:2:n_samples-1);
                    ILC(2:2:n_samples) = 1-I(2:2:n_samples);
                end
            end
            
            
            ILC = util.scale(ILC, [0,1], span);
            
           % ILC = ILC * (span(2) - span(1)) / 2;
           % ILC = ILC + span(1);
            

        end
    end    
end


>>>>>>> origin/master
