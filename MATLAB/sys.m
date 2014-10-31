
classdef sys
    properties
    end
    
    methods(Static)
    
        
        function [BIT_ARRAY] = to_bit_array(WORDS, nbits)
            nwords = size(WORDS,1);
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
        
        
        function [I, Q] = encode_IQ(WORDS, nbits, mode)
            
            nwords = size(WORDS,1);            
            BIT_ARRAY = sys.to_bit_array(WORDS, nbits);
            
            I = zeros(nwords, nbits/2);
            Q = zeros(nwords, nbits/2);
            
            if strcmp(mode, 'interleave')
                for i = nbits/2 : -1 : 1
                    I(:,i) = BIT_ARRAY(:,2*i-1);
                    Q(:,i) = BIT_ARRAY(:,2*i-0);
                end 
            elseif strcmp(mode, 'split')
                for i = nbits/2 : -1 : 1
                    I(:,i) = BIT_ARRAY(:,i);
                    Q(:,i) = BIT_ARRAY(:,nbits/2+i);
                end 
            else
                fprintf('map_to_IQ: unrecognized mode\n');
            end         
            
            I = sys.to_words(I, nbits/2);
            Q = sys.to_words(Q, nbits/2);
            
        end
        
        
        
        
        function [WORDS] = decode_IQ(I, Q, nbits, mode)
            
            nwords = size(I,1);            
            I_BITS = sys.to_bit_array(I, nbits/2);
            Q_BITS = sys.to_bit_array(Q, nbits/2);
            
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
            
            WORDS = sys.to_words(WORD_BITS, nbits);
                        
        end
        
        
    end    
end