clc;
clear all;
% close all;
s = rng(211);       % Set RNG state for repeatability

%Set parameters
numFFT = 512;           % Number of FFT points
numRBs = 20;             % Number of resource blocks
rbSize = 10;             % Number of subcarriers per resource block
cpLen = 32;              % Cyclic prefix length in samples

bitsPerSubCarrier = 8;   % 2: QPSK, 4: 16QAM, 6: 64QAM, 8: 256QAM
snrdB = 1:1:21;              % SNR in dB

toneOffset = 2.5;        % Tone offset or excess bandwidth (in subcarriers)
L = 513;                 % Filter length (=filterOrder+1), odd

%Filter Design
numDataCarriers = numRBs*rbSize;    % number of data subcarriers in sub-band
halfFilt = floor(L/2);
n = -halfFilt:halfFilt;

% Sinc function prototype filter
pb = sinc((numDataCarriers+2*toneOffset).*n./numFFT);

% Sinc truncation window
w = (0.5*(1+cos(2*pi.*n/(L-1)))).^0.6;

% Normalized lowpass filter coefficients
fnum = (pb.*w)/sum(pb.*w);

% Use dsp filter objects for filtering
filtTx = dsp.FIRFilter('Structure', 'Direct form symmetric', ...
    'Numerator', fnum);
filtRx = clone(filtTx); % Matched filter for the Rx


qamMapper = comm.RectangularQAMModulator('ModulationOrder', ...
    2^bitsPerSubCarrier, 'BitInput', true, ...
    'NormalizationMethod', 'Average power');

% Transmit-end processing
%  Initialize arrays
inpData = zeros(bitsPerSubCarrier*numRBs*rbSize,1);
txSig = complex(zeros(numFFT+L-1, 1));
snrdBarr =1:1:30;              % SNR in dB

save initpara
m=bitsPerSubCarrier*numRBs;
n1=rbSize;
frsz=m*n1
data11=randi([0 1],1320,1)
N3 = 132;
tempdata=[];
for clp=1:10

dataplr=data11((clp-1)*N3+1:clp*N3);
% Find the power of 2 that is equal to or greater than N
n2 = ceil(log2(N3));



% Polar encoding
encodedBits = nrPolarEncode(dataplr, 2^n2);
tempdata=[tempdata ;encodedBits(:)];
end


data1=[tempdata(:) ; zeros(frsz-rem(length(tempdata),frsz),1)];

for loopx=1:length(snrdBarr)


dataagg=[];
for dtl=1:ceil(length(data1)/frsz)
data=reshape(data1((dtl-1)*frsz+1:dtl*frsz),m*n1,1)

txSig=fofdm_mod(data)

% Compute peak-to-average-power ratio (PAPR)
pm = powermeter(Measurement="Peak-to-average power ratio",ComputeCCDF=true);
paprFOFDM = pm(txSig);
disp(['Peak-to-Average-Power-Ratio for F-OFDM = ' num2str(paprFOFDM) ' dB']);

% SETTING THE PARAMETERS FOR THE SIMULATION
fc=1e9;        %Carrier frequency
c=3e8;        %Speed of light
l=c/fc;        %Wavelength
d=l/2;        %Rx array spacing
N=32;         %Receive array size
M=16;         %Transmit array size (users)
theta=2*pi*(rand(1,M));         %Angular separation of users
n=1:N;                          %Rx array number
n=transpose(n);                 %Row vector to column vector
sigma=0.1;

% RECEIVE SIGNAL MODEL (LINEAR)
H=exp(-i*(n-1)*2*pi*d*cos(theta)/l);  %Channel matrix of size NxM
%  H=H/norm(H'*H);
H=(1/sqrt(2))*randn(N,M)+i*(1/sqrt(2))*randn(N,M);
w111=reshape([txSig],M,66);
x=H*w111;                             %Receive vector of length N
% Add WGN
tempx=zeros(size(w111));


rxSig = awgn(x, snrdBarr(loopx),'Measured');

% LINEAR ARRAY PROCESSING - METHODS
% 1-MATCHED FILTER
 y=pinv(H'*H)*H'*rxSig;

rxBits =fofdm_demod(y)

BER = comm.ErrorRate;

ber = BER(data(:), rxBits)

% disp(['UFMC Reception, BER = ' num2str(ber(1)) ' at SNR = ' ...
%     num2str(snrdBarr(loopx)) ' dB']);
ber_lp(dtl)=ber(1)

dataagg=[dataagg ;rxBits(:)];

end
tempdatad=[]
N1=2^n2;
for clp=1:10
% Polar decoding
decodedBits = nrPolarDecode((-100*(2*dataagg((clp-1)*N1+1:N1*clp)-1)), 132,256,8);
tempdatad=[tempdatad; decodedBits(:)];
end
SE(loopx)=(1-mean(ber_lp))*bitsPerSubCarrier;
ber1(loopx)=mean(ber_lp)
[nerr berc(loopx)]=biterr(data11(:),tempdatad(:))
end

 rng(s);
figure(100)

hold on
% semilogy(snrdBarr,smooth(ber1),'r')

semilogy(snrdBarr,smooth(berc),'g')
set(gca,'Yscale','log')
xlabel('SNR')
ylabel('BER')
grid on 

figure(101)
hold on ;
% ylim([0.5 0.5 10]);
plot(snrdBarr,smooth(SE),'k')
 set(gca,'Yscale','log')
grid on
xlabel('SNR')
ylabel('SE')
title('SE vs SNR for PC-F-OFDM Waveform')
