clc 
clear all
close all
%  rand('state',0);       % Set RNG state for repeatability
%s = rng(211); 
 %UFMC initial Parameters
numFFT = 512;        % number of FFT points
subbandSize = 20;    % must be > 1 
numSubbands = 10;    % numSubbands*subbandSize <= numFFT
subbandOffset = 156; % numFFT/2-subbandSize*numSubbands/2 for band center


%Channel Coding Parameters
coderate=1/2;
ConstraintLength=3;
CodeGenerator=[6 7];

% Dolph-Chebyshev window design parameters
filterLen = 43;      % similar to cyclic prefix length
slobeAtten = 40;     % side-lobe attenuation, dB
prototypeFilter = chebwin(filterLen, slobeAtten);

%QAM Order
bitsPerSubCarrier = 8;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM

% QAM Symbol mapper
qamMapper = comm.RectangularQAMModulator('ModulationOrder', ...
    2^bitsPerSubCarrier, 'BitInput', true, ...
    'NormalizationMethod', 'Average power');

% Transmit-end processing

%  Initialize arrays
inpData = zeros(bitsPerSubCarrier*subbandSize, numSubbands);
txSig = complex(zeros(numFFT+filterLen-1, 1));

snrdB = 1:2:32;

%Figure settings
hFig = figure;
axis([-0.5 0.5 -100 20]);
hold on; 
grid on

xlabel('Normalized frequency');
ylabel('PSD (dBW/Hz)')
title(['UFMC, ' num2str(numSubbands) ' Subbands, '  ...
    num2str(subbandSize) ' Subcarriers each'])
save initpara

data=randi([0 1], bitsPerSubCarrier*subbandSize, numSubbands)

%UFMC Modulator Function
txdata=ufmc_mod(data)
save txdata

set(hFig, 'Position', figposition([20 50 25 30]));
hold off;

PAPR = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
[~,~,paprUFMC] = PAPR(txdata);
disp(['Peak-to-Average-Power-Ratio (PAPR) for UFMC = ' num2str(paprUFMC) ' dB']);


% Massive MIMO SETTING THE PARAMETERS FOR THE SIMULATION
fc=1e9;        %Carrier frequency
c=3e8;        %Speed of light
l=c/fc;        %Wavelength
d=l/2;        %Rx array spacing
N=20;         %Receive array size
M=16;         %Transmit array size (users)
theta=2*pi*(rand(1,M));         %Angular separation of users
n=1:N;                          %Rx array number
n=transpose(n);                 %Row vector to column vector
sigma=0.1;

% RECEIVE SIGNAL MODEL (LINEAR)
H=exp(-i*(n-1)*2*pi*d*cos(theta)/l);  %Channel matrix of size NxM
%  H=H/norm(H'*H);
H=(1/sqrt(2))*randn(N,M)+i*(1/sqrt(2))*randn(N,M);
w111=reshape([txdata;zeros(6,1)],M,35);
x=H*w111;                             %Receive vector of length N

% Add WGN
 %  Loop over each snrdB
for loop=1:length(snrdB)
     rx = awgn(x, snrdB(loop), 'measured');
    

     
     % LINEAR ARRAY PROCESSING - METHODS
    
     %Zero forcing detector
     y=pinv(H'*H)*H'*rx;
     
     % 4-Minimum Mean Square Error (MMMSE)
     %y=(H'*H+(2*sigma^2)*eye([M,M]))^(-1)*(H')*rxSig;

     rxdata = ufmc_demod(y);
     
% Demapping and BER computation
qamDemod = comm.RectangularQAMDemodulator('ModulationOrder', ...
    2^bitsPerSubCarrier, 'BitOutput', true, ...
    'NormalizationMethod', 'Average power');
BER = comm.ErrorRate;

% Perform hard decision and measure errors
rxSig = qamDemod(rxdata);
ber = BER(inpData(:), rxSig)
ber1(loop)=ber(1);
disp(['UFMC Reception, BER = ' num2str(ber(1)) ' at SNR = ' ...
    num2str(snrdB(loop)) ' dB']);
end
% Restore RNG state
%rng(s);
figure (100)
hold on
semilogy(snrdB,smooth(ber1),'r')
set(gca,'Yscale','log')
grid on
xlabel('SNR')
ylabel('BER')


