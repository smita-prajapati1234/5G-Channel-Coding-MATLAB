function rxSig=ufmc_demodcc(y)
load initpara



 y1=y(:);
 y1(end-5:end)=[];

  yRxPadded = [y1; zeros(2*numFFT-numel(txSig),1)];

% Perform FFT and downsample by 2
RxSymbols2x = fftshift(fft(yRxPadded));
RxSymbols = RxSymbols2x(1:2:end);

% Select data subcarriers
dataRxSymbols = RxSymbols(subbandOffset+(1:numSubbands*subbandSize));



% Use zero-forcing equalizer after OFDM demodulation
rxf = [prototypeFilter.*exp(1i*2*pi*0.5*(0:filterLen-1)'/numFFT); ...
       zeros(numFFT-filterLen,1)];
prototypeFilterFreq = fftshift(fft(rxf));
prototypeFilterInv = 1./prototypeFilterFreq(numFFT/2-subbandSize/2+(1:subbandSize));

% Equalize per subband - undo the filter distortion
dataRxSymbolsMat = reshape(dataRxSymbols,subbandSize,numSubbands);
EqualizedRxSymbolsMat = bsxfun(@times,dataRxSymbolsMat,prototypeFilterInv);
EqualizedRxSymbols = EqualizedRxSymbolsMat(:);


% Demapping and BER computation
qamDemod = comm.RectangularQAMDemodulator('ModulationOrder', ...
    2^bitsPerSubCarrier, 'BitOutput', true, ...
    'NormalizationMethod', 'Average power');


% Perform hard decision and measure errors
rxSig1 = qamDemod(EqualizedRxSymbols);

temp2 = zeros(bitsPerSubCarrier*subbandSize*coderate, numSubbands);
temp=reshape(rxSig1,size(inpData1));
for lp=1:size(inpData1,2)
    temp11 = randdeintrlv(temp(:,lp),1008);
temp2(:,lp)=vitdec(temp11,poly2trellis(ConstraintLength,CodeGenerator),32,'trunc','hard');
end
rxSig=temp2(:);
end

