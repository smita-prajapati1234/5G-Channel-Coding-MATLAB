function rxBits=fofdm_demod(y)
load initpara
offset = (numFFT-numDataCarriers)/2; % for band center
% y1=y(:);
% y1(end-5:end)=[];
y1=y(:);
% Receive matched filter
rxSigFilt = filtRx(y1);

% Account for filter delay
rxSigFiltSync = rxSigFilt(L:end);

% Remove cyclic prefix
rxSymbol = rxSigFiltSync(cpLen+1:end);

% Perform FFT
RxSymbols = fftshift(fft(rxSymbol));

% Select data subcarriers
dataRxSymbols = RxSymbols(offset+(1:numDataCarriers));
rxBits = qamdemod(dataRxSymbols, 2^bitsPerSubCarrier, 'OutputType', 'bit', ...
    'UnitAveragePower', true);
end