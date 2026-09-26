function rxSig=ufmc_demod(y)
numFFT = 512;        % number of FFT points
subbandSize = 20;    % must be > 1 
numSubbands = 10;    % numSubbands*subbandSize <= numFFT
subbandOffset = 156; % numFFT/2-subbandSize*numSubbands/2 for band center
load txdata


 y1=y(:);
 y1(end-5:end)=[];

  yRxPadded = [y1; zeros(2*numFFT-numel(txdata),1)];

% Perform FFT and downsample by 2
RxSymbols2x = fftshift(fft(yRxPadded));
RxSymbols = RxSymbols2x(1:2:end);

% Select data subcarriers
dataRxSymbols = RxSymbols(subbandOffset+(1:numSubbands*subbandSize));

% Plot received symbols constellation
 constDiagRx = comm.ConstellationDiagram('ShowReferenceConstellation', ...
     false, 'Position', figposition([20 15 25 30]), ...
     'Title', 'UFMC Pre-Equalization Symbols', ...
    'Name', 'UFMC Reception', ...
    'XLimits', [-150 150], 'YLimits', [-150 150]);
 constDiagRx(dataRxSymbols);

% Use zero-forcing equalizer after OFDM demodulation
rxf = [prototypeFilter.*exp(1i*2*pi*0.5*(0:filterLen-1)'/numFFT); ...
       zeros(numFFT-filterLen,1)];
prototypeFilterFreq = fftshift(fft(rxf));
prototypeFilterInv = 1./prototypeFilterFreq(numFFT/2-subbandSize/2+(1:subbandSize));

% Equalize per subband - undo the filter distortion
dataRxSymbolsMat = reshape(dataRxSymbols,subbandSize,numSubbands);
EqualizedRxSymbolsMat = bsxfun(@times,dataRxSymbolsMat,prototypeFilterInv);
EqualizedRxSymbols = EqualizedRxSymbolsMat(:);

% Plot equalized symbols constellation
 constDiagEq = comm.ConstellationDiagram('ShowReferenceConstellation', ...
     false, 'Position', figposition([46 15 25 30]), ...
     'Title', 'UFMC Equalized Symbols', ...
     'Name', 'UFMC Equalization');
 constDiagEq(EqualizedRxSymbols);


rxSig=EqualizedRxSymbols;
end
