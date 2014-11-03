
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   WiCam Wireless Video Camera System Simulation
%
%   This program simulates an end-to-end digital transmission testbench.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Definitions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   application     the HW and SW that uses this transmission system
%
%   baseband        a low-frequency modulated carrier. A baseband carrier 
%                       of frequency f can be modulated by up to +/- f, 
%                       giving a channel bandwidth of 2f.
%
%   data stream     a sequence of samples at the application level
%                       (e.g. video feed)
%
%   I/Q             In-phase / Quadrature carriers. Q is phase-offset from 
%                       I by 90 degrees. Because they are orthogonal, each
%                       can be modulated and demodulated independantly, 
%                       which doubles the number of symbols that can be
%                       transmitted per cycle
%
%   line code       a type of lossless encoding that transforms a signal    
%                       into a format suitable for modulating a carrier
%
%   map             the conversion between a data stream sample to/from 
%                       I/Q samples. An N-bit stream sample will map to K 
%                       M-bit I/Q sample pairs, where 2^M is the number of 
%                       possible amplitude levels of that I/Q component
%
%   sample          a data point representing the value of some signal at 
%                       one point in time
%
%   sequence        an array of discrete samples
%
%   source          the origin of application data (e.g. camera)
%
%   symbol          a sample that modulates a carrier. An N-bit symbol 
%                       corresponds to 2^N levels of modulation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Signal path
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [TX-SOURCE]
%    1) Generate source data stream samples and save to file
%    2) Read back data stream samples from file into sample array
% [TX-FPGA]
%    3) Map data stream samples into I/Q samples
%    4) Line-code I/Q samples
%    5) Generate baseband samples from line-coded I/Q
%    6) Filter/condition/pre-emph baseband samples
% [TX-DAC]
%    7) Convert conditioned baseband samples to baseband analog I/Q
% [TX-FILTER]
%    8) Filter baseband analog I/Q
% [TX-RF]
%    9) Low-pass-filter baseband I/Q
%   10) Direct-upconvert filtered baseband I/Q to 2.4 GHz RF
%   11) Amplify 2.4 GHz RF and send to antenna
% [AIR]
%   12) Atentuate 2.4 GHz RF and subject to environmental noise
% [RX-RF]
%   13) Receive 2.4 GHz RF on antenna and amplify
%   14) Direct-downconvert 2.4 GHz RF to baseband I/Q
%   15) Low-pass filter baseband I/Q
% [RX-ADC]
%   16) Convert baseband I/Q into digital line-coded I/Q samples
% [RX-FPGA]
%   17) Decode line-coded I/Q samples
%   18) Recombine I/Q samples into data stream samples
% [RX-SINK]
%   19) Save data stream to file
% [FINALLY]
%   20) Compare source and sink data and report errors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% TODO: define signal structure
% struct Signal {
%   samples
%   sample_rate
%   symbol_rate
%   bits_per_symbol
% }


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Always start from scratch 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;      % clear workspace
clc;            % clear console

rng(0);         % seed random number generator with a constant


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Universal constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% What is this?
pi = 3.1416;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Duration of one time-domain instant
td_resolution = 5e-9;  

% Colours for plotting
black = [0.0 0.0 0.0];
blue = [0.0 0.0 1.0];
green = [0.0 0.7 0.0]; 


% Plot colours
source_plot_color = black;
I_plot_color = blue;
Q_plot_color = green;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datafile = 'data_in.txt';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modulation scheme comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All modulation types supported (or TODO)
mod_types = ['256qam'; '64qam '; '16qam '; 'qpsk  '; 'bpsk  '];

% Possible channel bandwidths (programmable in the RF IC)
bandwidths = [10e6; 15e6; 20e6];

% Number of baseband cycles per I/Q symbol
% Optionally reduce throughput for increased noise tolerance
cycles_per_symbols = [0.5; 1; 2; 3; 4; 5; 6];

% Only report modulation schemes that can achieve this bitrate.
% (For this section only; actual system values are defined later)
min_rf_bitrate = 20e6; % 20Mbps = 320 x 240 x 4 bits x 60/sec

% Generate a report comparing all combinations of modulation, bandwidth, 
% and cycles per symbol that meet the specified minimum bitrate
rf.compare_modulation_types( ...
    min_rf_bitrate, mod_types, bandwidths, cycles_per_symbols );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modulation specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the modulation scheme. See mod_types above for supported schemes.
modulation_type = '256qam';
          
% Define the channel bandwidth. See bandwidths above for allowed values.
bandwidth = 10e6;

% Define number of baseband cycles per I/Q symbol. See symbols_per_cycles
% above for allowed values. Be aware that the chosen modulation type may 
% not support all values. e.g. QAM requires a whole integer value
%baseband_cycles_per_iq_symbol = 1;


% Get the channel parameters for the chosen modulation type
[n_channels, bits_per_iq_symbol] = rf.get_modulation_format(modulation_type);

baseband_carrier_frequency = bandwidth / 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RF I/O specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nominal I/Q input common-mode voltage
rf_iq_vcm_in = 1.1; % 1.1Vcm

% Nominal I/Q input voltage
rf_iq_vrms_in = 100e-3; % 100mVrms



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source data sample format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bits per sample from the "camera"
bits_per_source_sample = 8;    

% Range of possible source values
source_sample_span = [0 2^bits_per_source_sample-1];  

% Number of source samples per second
source_sample_rate = 5e6;  % 2.5 Msps   

% Source bits per second
source_data_rate = bits_per_source_sample * source_sample_rate;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I/Q sample format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I/Q bit mapping method
%IQ_mode = 'interleave';   % Even bits to I, odd bits to Q
%IQ_mode = 'split';        % High bits to I, low bits to Q 
 

% iq_samples_per_symbol =  iq_samples_per_cycle / iq_symbols_per_cycle;


% bits_per_iq_sample = bits_per_iq_symbol / iq_samples_per_symbol;

bits_per_iq_sample = bits_per_iq_symbol;

% Range of possible I/Q component values
iq_sample_span = [0 2^bits_per_iq_sample-1];        


% Range of possible I/Q component values
iq_samples_per_source_sample = bits_per_source_sample / (bits_per_iq_sample * n_channels);

% Number of I/Q samples per second
iq_sample_rate = source_sample_rate * iq_samples_per_source_sample;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DAC input format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of DAC input bits per channel
n_dac_bits = bits_per_iq_sample;

% Range of possible DAC input values
dac_input_span = [0 2^n_dac_bits-1];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DAC output format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Desired DAC output voltage swing (based on RF I/Q input spec)
dac_full_scale_voltage = rf_iq_vrms_in * 2 * sqrt(2);

% TODO: implement current-output DAC

% Nominal DAC full-scale output current
%dac_full_scale_current = 8e-3;  % 8mA

% Calculate ideal load resistance
%dac_load_resistance = dac_full_scale_voltage / dac_full_scale_current;

% Range of possible DAC output voltages
dac_output_voltage_span = rf_iq_vcm_in + ... 
    [-dac_full_scale_voltage/2 dac_full_scale_voltage/2];

% Number of DAC samples per second
dac_sample_rate = 80e6;  % 80 Msps   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [TX-SOURCE]
%    1) Generate source data stream samples and save to file
%    2) Read back data stream samples from file into sample array
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    1) Generate source data stream samples and save to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate and combine a variety of signals
X = [];
X = [X; f.SigGen('rect', 32, source_sample_span, 16, 0, 0.5)];    
X = [X; f.SigGen('tri',  64, source_sample_span, 1.5, 0, 0.5)];
X = [X; f.SigGen('sin',  64, source_sample_span, 1.5, -0.25)];
X = [X; f.SigGen('rand', 96, source_sample_span)];

% Clip to span bounds and round to discrete integers
X = min(X, source_sample_span(2));
X = max(X, source_sample_span(1));
X = round(X);

% Write samples to file one per line as readable hex (e.g. 00, 01, FF)
chars_per_source_sample = ceil(bits_per_source_sample / 4);
fmt = sprintf('%%0%dX\\n', chars_per_source_sample);
fileId = fopen(datafile,'w');
fprintf(fileId,fmt, X);
fclose(fileId);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    2) Read back data stream samples from file into sample array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read input data file of readable hex
fmt = sprintf('%%%dc', chars_per_source_sample);
fileId = fopen(datafile,'r');
X = hex2dec(textscan(fileId, fmt));
fclose(fileId);

% Plot input stream
X_t = util.to_time_domain(X, source_sample_rate, td_resolution, ...
    source_sample_span, source_sample_span);
gfx.figDiscreteAndTD('Input Stream', X, X_t, ...
    bits_per_source_sample, source_sample_rate, ...
    source_sample_span, source_sample_span, source_plot_color);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [TX-FPGA]
%    3) Map data stream samples into I/Q samples
%    4) Line-code I/Q samples
%    5) Generate baseband samples from line-coded I/Q
%    6) Filter/condition/pre-emph baseband samples
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    3) Map data stream samples into I/Q samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%[IS, QS] = fpga.map_IQ(X, bits_per_source_sample, bits_per_iq_sample, n_channels); 

[IS, QS] = fpga.map_IQ_8bit_256qam(X); 




IS_t = util.to_time_domain(IS, iq_sample_rate, td_resolution, ...
    iq_sample_span, iq_sample_span);
QS_t = util.to_time_domain(QS, iq_sample_rate, td_resolution, ...
    iq_sample_span, iq_sample_span);

gfx.figDiscreteAndTD('I Samples', IS, IS_t, bits_per_iq_sample, ...
    iq_sample_rate, iq_sample_span, iq_sample_span, I_plot_color);
gfx.figDiscreteAndTD('Q Samples', QS, QS_t, bits_per_iq_sample, ...
    iq_sample_rate, iq_sample_span, iq_sample_span, Q_plot_color);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    4) Line-code I/Q samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%symbols_per_cycle = 1;
%always_cross_zero = 1;

%ILC = fpga.line_encode(I, iq_sample_span, symbols_per_cycle, always_cross_zero);
%QLC = fpga.line_encode(Q, iq_sample_span, symbols_per_cycle, always_cross_zero);

iqlc_sample_span = iq_sample_span;

iqlc_samples_per_iq_sample = 1;
iqlc_sample_rate = iq_sample_rate * iqlc_samples_per_iq_sample;
%iqlc_samples_per_baseband_cycle = ...
%    iqlc_samples_per_iq_sample / baseband_cycles_per_iq_symbol;


bits_per_iqlc_sample = bits_per_iq_sample;

ILC = IS;
QLC = QS;

ILCN = ILC;
QLCN = QLC;

ILC_t = util.to_time_domain(ILC, iqlc_sample_rate, ...
    td_resolution, iq_sample_span, iq_sample_span);
ILCN_t = util.to_time_domain(ILCN, iqlc_sample_rate, ...
    td_resolution, iq_sample_span, iq_sample_span);

QLC_t = util.to_time_domain(QLC, iqlc_sample_rate, ...
td_resolution, iq_sample_span, iq_sample_span);
QLCN_t = util.to_time_domain(QLCN, iqlc_sample_rate, ...
    td_resolution, iq_sample_span, iq_sample_span);

gfx.figDiscreteAndTD('I Line-Coded', ILC, ILC_t, bits_per_iqlc_sample, ...
iqlc_sample_rate, iq_sample_span, iq_sample_span, I_plot_color);
gfx.figDiscreteAndTD('Q Line-Coded', QLC, QLC_t, bits_per_iqlc_sample, ...
    iqlc_sample_rate, iq_sample_span, iq_sample_span, Q_plot_color);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    5) Generate baseband samples from line-coded I/Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

baseband_sample_rate = dac_sample_rate;

baseband_samples_per_cycle = ...
    baseband_sample_rate / baseband_carrier_frequency;


iqlc_samples_per_baseband_cycle = source_data_rate / (bits_per_iqlc_sample * baseband_carrier_frequency * n_channels);


baseband_samples_per_iqlc_sample = ...
    baseband_samples_per_cycle / iqlc_samples_per_baseband_cycle;

n_baseband_samples = numel(ILC) * baseband_samples_per_iqlc_sample;

n_baseband_cycles = n_baseband_samples / baseband_samples_per_cycle;

IC = f.SigGen('sin',n_baseband_samples,[-1,1],n_baseband_cycles,0,0.5);
QC = f.SigGen('rect',n_baseband_samples,[-1,1],n_baseband_cycles,0.25,0.5);


IC_t = util.to_time_domain(IC, baseband_sample_rate, ...
td_resolution, [-1,1], [-1,1]);
QC_t = util.to_time_domain(QC, baseband_sample_rate, ...
td_resolution, [-1,1], [-1,1]);

gfx.figDiscreteAndTD('Baseband I Carrier', IC, IC_t, bits_per_iqlc_sample / iqlc_samples_per_baseband_cycle, ...
    baseband_sample_rate, [-1,1], [-1,1], I_plot_color);
gfx.figDiscreteAndTD('Baseband Q Carrier', QC, QC_t, bits_per_iqlc_sample / iqlc_samples_per_baseband_cycle, ...
    baseband_sample_rate, [-1,1], [-1,1], Q_plot_color);

I = util.to_time_domain(ILC, iqlc_sample_rate, 1/baseband_sample_rate, ...
    iqlc_sample_span, [-1,1]);
Q = util.to_time_domain(QLC, iqlc_sample_rate, 1/baseband_sample_rate, ...
    iqlc_sample_span, [-1,1]);

% Ignore the timescale information generated by to_time_domain
I = I(:,2);
Q = Q(:,2);

n = size(I,1);
q_offset = baseband_samples_per_iqlc_sample/4;

Q = [zeros(q_offset,1); Q(1:n-q_offset)];



I = I .* IC;
Q = Q .* QC;

I_t = util.to_time_domain(I, dac_sample_rate, ...
td_resolution, [-1,1], [-1,1]);
Q_t = util.to_time_domain(Q, dac_sample_rate, ...
td_resolution, [-1,1], [-1,1]);

gfx.figDiscreteAndTD('Baseband I', I, I_t, bits_per_iqlc_sample / iqlc_samples_per_baseband_cycle, ...
    baseband_sample_rate, [-1,1], [-1,1], I_plot_color);
gfx.figDiscreteAndTD('Baseband Q', Q, Q_t, bits_per_iqlc_sample / iqlc_samples_per_baseband_cycle, ...
    baseband_sample_rate, [-1,1], [-1,1], Q_plot_color);






% 
% [freq,ILC_FFT] = f.FFT2(ILC, symbols_per_cycle*iq_sample_rate);
% 
% figure(7)
% title = 'I Line-Coded FFT'; 
% set(gcf,'numbertitle','off','name',title);
% subplot(2, 1, 1);
% f.plotFFT(real(ILC_FFT), symbols_per_cycle*iq_sample_rate, 'Magnitude');
% subplot(2, 1, 2);
% f.plotFFT(imag(ILC_FFT), symbols_per_cycle*iq_sample_rate, 'Phase');
% 
% 
% [freq,QLC_FFT] = f.FFT2(QLC, symbols_per_cycle*iq_sample_rate);
% 
% figure(8)
% title = 'Q Line-Coded FFT'; 
% set(gcf,'numbertitle','off','name',title);
% subplot(2, 1, 1);
% f.plotFFT(real(QLC_FFT), symbols_per_cycle*iq_sample_rate, 'Magnitude');
% subplot(2, 1, 2);
% f.plotFFT(imag(QLC_FFT), symbols_per_cycle*iq_sample_rate, 'Phase');
% 


%ILC_t(:,2) = cm_IQ + ILC_t(:,2);
%ILCN_t(:,2) = cm_IQ + ILCN_t(:,2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    6) Filter/condition/pre-emph line-coded I/Q samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%ILC_t = filter(a, [1 a-1], ILC_t);

%Y_t = conv(ILC_t(:,2), S);
%ILC_t(:,2) = Y_t(1:size(ILC_t,1));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [TX-DAC]
%    7) Convert conditioned I/Q samples to baseband analog I/Q
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [TX-RF]
%    7) Low-pass-filter baseband I/Q
%    8) Direct-upconvert filtered baseband I/Q to 2.4 GHz RF
%    9) Amplify 2.4 GHz RF and send to antenna
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [AIR]
%   10) Atentuate 2.4 GHz RF and subject to environmental noise
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [RX-RF]
%   11) Receive 2.4 GHz RF on antenna and amplify
%   12) Direct-downconvert 2.4 GHz RF to baseband I/Q
%   13) Low-pass filter baseband I/Q
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [RX-ADC]
%   14) Convert baseband I/Q into digital line-coded I/Q samples
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [RX-FPGA]
%   15) Decode line-coded I/Q samples
%   16) Recombine I/Q samples into data stream samples
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Y = fpga.unmap_IQ(I, Q, bits_per_source_sample, IQ_mode)

%stem(Y, 'MarkerFaceColor',[0 0 1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [RX-SINK]
%   17) Save data stream to file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [FINALLY]
%   18) Compare source and sink data and report errors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% if Y ~= X
%     fprintf('RX errors found\n');
% else
%     fprintf('No RX errors found\n');
% end
 
 
 