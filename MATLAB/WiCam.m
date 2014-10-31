
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% WiCam Wireless Video Camera System Simulation
%
% This program simulates a full end-to-end digital transmission testbench.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Signal path:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [TX-SOURCE]
%    1) Source data stream is generated and saved to file
%    2) Data stream is read back from file into a sample array
% [TX-FPGA]
%    3) Samples are seperated into I/Q components
%    4) I/Q samples are line-coded 
%    5) Line-coded I/Q samples are filtered/conditioned/pre-emph'd.
% [TX-DAC]
%    6) Conditioned I/Q samples are converted to baseband analog I/Q
% [TX-RF]
%    7) Baseband is low-pass-filtered
%    8) Filtered baseband is direct-upconverted to 2.4 GHz RF
%    9) RF is amplified and sent to antenna
% [AIR]
%   10) RF is attentuated and subjected to environmental noise
% [RX-RF]
%   11) RF is received on antenna and amplified
%   12) RF is direct-downconverted to baseband I/Q
%   13) Baseband is low-pass filtered 
% [RX-ADC]
%   14) Baseband is converted into digital I/Q samples
% [RX-FPGA]
%   15) I/Q samples are decoded
%   16) Decoded I/Q samples are recombined into the original data stream
% [RX-SINK]
%   17) Data is saved to file
% [FINALLY]
%   18) Compare source and sink data and report errors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Always start from scratch 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;      % clear workspace
clc;            % clear console


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Universal constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pi = 3.1416;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Digital data format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bits_per_sample = 8;    % 

bits_per_IQ = bits_per_sample / 2;

fsd_S = 2^bits_per_sample-1;  % Sample full-scale digital value

fsd_IQ = 2^bits_per_IQ-1;     % I/Q full-scale digital value
fsa_IQ = 1;                 % I/Q differential full-scale voltage
cm_IQ = 1;                 % I/Q common-mode voltage


sample_rate = 20e6;     % Number of I/Q samples per second

td_res = 1e-9;          % Time-domain resolution: duration of one point


% Number of time-domain data points per discrete sample
pts_per_sample = 1/(sample_rate*td_res);   


IQ_mode = 'interleave';
%IQ_mode = 'split';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a test data file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datafile = 'data_in.txt';

xpos = 1;   % Position of the next data sample

% 1) Clock burst
clock_burst_len = 32;
F = [];
F(1:2:clock_burst_len) = 0;
F(2:2:clock_burst_len) = 1;
X(xpos:xpos+clock_burst_len-1) = F;
xpos = xpos + clock_burst_len;

% 2) Ramp
ramp_len = 32;
F = [];
F = linspace(0,1,ramp_len);
X(xpos:xpos+ramp_len-1) = F;
xpos = xpos + ramp_len;

% 3) Quadratic growth
quad_growth_len = 32;
F = [];
F = linspace(0,1,quad_growth_len);
F = F.^2;
X(xpos:xpos+quad_growth_len-1) = F;
xpos = xpos + quad_growth_len;

% 4) Exponential growth
exp_growth_len = 32;
exp_growth_start = log(1/bits_per_sample);
exp_growth_stop = 3;
F = [];
F = linspace(exp_growth_start,exp_growth_stop,exp_growth_len);
F = exp(F) / exp(exp_growth_stop);
X(xpos:xpos+exp_growth_len-1) = F;
xpos = xpos + exp_growth_len;

% 5) Sinusoid
sinusoid_len = 64;
n_periods = 1.5;
F = [];
F = linspace(0,2*pi*n_periods,sinusoid_len);
F = 0.5*(1-0.5*cos(F));
X(xpos:xpos+sinusoid_len-1) = F;
xpos = xpos + sinusoid_len;

% 6) Random
rng(0);
random_len = 100;
F = [];
F = rand(random_len,1);
X(xpos:xpos+random_len-1) = F;
xpos = xpos + random_len;


% Scale to full-scale and round to integers
X = X .* fsd_S;
X = round(X);
X = mod(X, 2^bits_per_sample);

% Write to file
fileId = fopen(datafile,'w');
fprintf(fileId,'%02X\n', X);
fclose(fileId);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process test data file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Read input data file
fileId = fopen(datafile,'r');
X=hex2dec(textscan(fileId,'%2c'))    % readable hex data
fclose(fileId);
 
figure(1);
title = 'Input Data'; 
set(gcf,'numbertitle','off','name',title);
subplot(2, 1, 1);
stem(X, 'Color',[0 0 0], 'MarkerFaceColor',[0 0 0]);
axis([0,+inf,0,255]);

X_t = f.to_time_domain(X, sample_rate, pts_per_sample, fsd_S, fsd_S);



% Generate a windowed sinc
td_filter_len = 300;   
td_filter_speed = 5;
S = sinc(linspace(-td_filter_speed*pi,td_filter_speed*pi,td_filter_len+1));
S = S / (sum(S));
W = f.HammingWindow(0,td_filter_len,0);   
S = S .* W;


figure(2);
title = 'Sinc Window'; 
set(gcf,'numbertitle','off','name',title);
subplot(2, 1, 1);
stem(S, 'Color',[0 0 0], 'MarkerFaceColor',[0 0 0]);
subplot(2, 1, 2);
plot(S);


%Y_t = conv(X_t(:,2), S);
%X_t(:,2) = Y_t(1:size(X_t,1));



fc = 22.5e6;
tau = 1/(2*pi*fc);
a = 1/(sample_rate*pts_per_sample*tau);

%X_t = filter(a, [1 a-1], X_t);

figure(1);
subplot(2, 1, 2);
plot(X_t(:,1),X_t(:,2), 'Color',[0 0 0], 'LineWidth',1);
axis([0,+inf,0,255]);


[I, Q] = sys.encode_IQ(X, bits_per_sample, IQ_mode); 
IQ = [I, Q];


figure(3);
title = 'I Samples'; 
set(gcf,'numbertitle','off','name',title);
subplot(2, 1, 1);
stem(I, 'Color',[0 0 1], 'MarkerFaceColor',[0 0 1], 'LineWidth', 1, 'MarkerSize', 7);
axis([0,+inf,0,fsd_IQ]);

figure(4);
title = 'Q Samples'; 
set(gcf,'numbertitle','off','name',title);
subplot(2, 1, 1);
stem(Q, 'Color',[0 .7 0], 'MarkerFaceColor',[0 .8 0], 'LineWidth', 1, 'MarkerSize', 7);
axis([0,+inf,0,fsd_IQ]);

I_t = f.to_time_domain(I, sample_rate, pts_per_sample, fsd_IQ, fsa_IQ);
%I_t = filter(a, [1 a-1], I_t);

Q_t = f.to_time_domain(Q, sample_rate, pts_per_sample, fsd_IQ, fsa_IQ);
%Q_t = filter(a, [1 a-1], Q_t);

figure(3);
subplot(2, 1, 2);
plot(I_t(:,1), I_t(:,2), 'Color',[0 0 1], 'LineWidth',1);
axis([0,+inf,0,fsa_IQ]);

figure(4);
subplot(2, 1, 2);
plot(Q_t(:,1), Q_t(:,2), 'Color',[0 .7 0], 'LineWidth',1);
axis([0,+inf,0,fsa_IQ]);



n_samples = size(I,1);


symbols_per_cycle = 1;

always_cross_zero = 1;

if symbols_per_cycle == 1
    
    if always_cross_zero == 1
        ILC = zeros(2*n_samples,1);
        ILC(1:2:2*n_samples-1) = (fsd_IQ+(1+I))/2;
        ILC(2:2:2*n_samples) = (fsd_IQ-(1+I))/2;

        QLC = zeros(2*n_samples,1);
        QLC(1:2:2*n_samples-1) = (fsd_IQ+(1+Q))/2;
        QLC(2:2:2*n_samples) = (fsd_IQ-(1+Q))/2;  
    
    else
        ILC = zeros(2*n_samples,1);
        ILC(1:2:2*n_samples-1) = I;
        ILC(2:2:2*n_samples) = fsd_IQ-I;

        QLC = zeros(2*n_samples,1);
        QLC(1:2:2*n_samples-1) = Q;
        QLC(2:2:2*n_samples) = fsd_IQ-Q;
        
    end    
        
elseif symbols_per_cycle == 2
    
    if always_cross_zero == 1
        ILC = zeros(n_samples,1);
        ILC(1:2:n_samples-1) = (fsd_IQ+(1+I(1:2:n_samples-1)))/2;
        ILC(2:2:n_samples) = (fsd_IQ-(1+I(2:2:n_samples)))/2;

        QLC = zeros(n_samples,1);
        QLC(1:2:n_samples-1) = (fsd_IQ+(1+Q(1:2:n_samples-1)))/2;
        QLC(2:2:n_samples) = (fsd_IQ-(1+Q(2:2:n_samples)))/2;
        
    else        
        ILC = zeros(n_samples,1);
        ILC(1:2:n_samples-1) = I(1:2:n_samples-1);
        ILC(2:2:n_samples) = fsd_IQ-I(2:2:n_samples);

        QLC = zeros(n_samples,1);
        QLC(1:2:n_samples-1) = Q(1:2:n_samples-1);
        QLC(2:2:n_samples) = fsd_IQ-Q(2:2:n_samples);
        
    end
end




ILCN = ILC;
QLCN = QLC;


ILC_t = f.to_time_domain(ILC, symbols_per_cycle*sample_rate, pts_per_sample, fsd_IQ, fsa_IQ);
ILCN_t = f.to_time_domain(ILCN, symbols_per_cycle*sample_rate, pts_per_sample, fsd_IQ, fsa_IQ);


QLC_t = f.to_time_domain(QLC, symbols_per_cycle*sample_rate, pts_per_sample, fsd_IQ, fsa_IQ);
QLCN_t = f.to_time_domain(QLCN, symbols_per_cycle*sample_rate, pts_per_sample, fsd_IQ, fsa_IQ);


%ILC_t = filter(a, [1 a-1], ILC_t);



%Y_t = conv(ILC_t(:,2), S);
%ILC_t(:,2) = Y_t(1:size(ILC_t,1));




figure(5);
title = 'I Line-Coded'; 
set(gcf,'numbertitle','off','name',title);
subplot(2, 1, 1);
stem(ILC, 'Color',[0 0 1], 'MarkerFaceColor',[0 0 1], 'LineWidth', 1, 'MarkerSize', 7);
axis([0,+inf,0,fsd_IQ]);
subplot(2, 1, 2);
plot(ILC_t(:,1), ILC_t(:,2), 'Color',[0 0 1], 'LineWidth',1);
axis([0,+inf,0,fsa_IQ]);


figure(6);
title = 'Q Line-Coded'; 
set(gcf,'numbertitle','off','name',title);
subplot(2, 1, 1);
stem(QLC, 'Color',[0 .7 0], 'MarkerFaceColor',[0 .7 0], 'LineWidth', 1, 'MarkerSize', 7);
axis([0,+inf,0,fsd_IQ]);
subplot(2, 1, 2);
plot(QLC_t(:,1), QLC_t(:,2), 'Color',[0 .7 0], 'LineWidth',1);
axis([0,+inf,0,fsa_IQ]);




[freq,ILC_FFT] = f.FFT2(ILC, symbols_per_cycle*sample_rate);

figure(7)
title = 'I Line-Coded FFT'; 
set(gcf,'numbertitle','off','name',title);
subplot(2, 1, 1);
f.plotFFT(real(ILC_FFT), symbols_per_cycle*sample_rate, 'Magnitude');
subplot(2, 1, 2);
f.plotFFT(imag(ILC_FFT), symbols_per_cycle*sample_rate, 'Phase');


[freq,QLC_FFT] = f.FFT2(QLC, symbols_per_cycle*sample_rate);

figure(8)
title = 'Q Line-Coded FFT'; 
set(gcf,'numbertitle','off','name',title);
subplot(2, 1, 1);
f.plotFFT(real(QLC_FFT), symbols_per_cycle*sample_rate, 'Magnitude');
subplot(2, 1, 2);
f.plotFFT(imag(QLC_FFT), symbols_per_cycle*sample_rate, 'Phase');





%ILC_t(:,2) = cm_IQ + ILC_t(:,2);
%ILCN_t(:,2) = cm_IQ + ILCN_t(:,2);



Y = sys.decode_IQ(I, Q, bits_per_sample, IQ_mode)


%stem(Y, 'MarkerFaceColor',[0 0 1]);

if Y ~= X
    fprintf('RX errors found\n');
else
    fprintf('No RX errors found\n');
end
 
 
 