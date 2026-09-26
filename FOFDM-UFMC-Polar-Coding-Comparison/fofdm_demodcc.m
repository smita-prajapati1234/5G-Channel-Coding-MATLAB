function rxSig=fofdm_demodcc(y)
load initpara
offset = (numFFT-numDataCarriers)/2; % for band center
y1=y(:);
%  y1(end-7:end)=[];

% Receive matched filter
rxSigFilt = filtRx(y1);

% Account for filter delay
rxSigFiltSync = rxSigFilt(filterLen:end);

% Remove cyclic prefix
rxSymbol = rxSigFiltSync(cpLen+1:end);

% Perform FFT
RxSymbols = fftshift(fft(rxSymbol));

% Select data subcarriers
dataRxSymbols = RxSymbols(offset+(1:numDataCarriers));
rxBits = qamdemod(dataRxSymbols, 2^bitsPerSubCarrier, 'OutputType', 'bit', ...
    'UnitAveragePower', true);

temp2 = zeros(bitsPerSubCarrier*numRBs*coderate, rbSize);
temp=reshape(rxBits,size(inpData1));
for lp=1:size(inpData1,2)
    temp11 = randdeintrlv(temp(:,lp),1008);
temp2(:,lp)=vitdec(temp11,poly2trellis(ConstraintLength,CodeGenerator),32,'trunc','hard');
end
rxSig=temp2(:);
end