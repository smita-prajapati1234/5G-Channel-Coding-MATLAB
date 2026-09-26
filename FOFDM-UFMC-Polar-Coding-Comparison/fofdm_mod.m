function txSig=fofdm_mod(data)
load initpara

% QAM Symbol mapper
symbolsIn = qammod(data, 2^bitsPerSubCarrier, 'InputType', 'bit', ...
    'UnitAveragePower', true);

% Pack data into an OFDM symbol
offset = (numFFT-numDataCarriers)/2; % for band center
symbolsInOFDM = [zeros(offset,1); symbolsIn; ...
    zeros(numFFT-offset-numDataCarriers,1)];
ifftOut = ifft(ifftshift(symbolsInOFDM));

% Prepend cyclic prefix
txSigOFDM = [ifftOut(end-cpLen+1:end); ifftOut];

% Filter, with zero-padding to flush tail. Get the transmit signal
txSigFOFDM = filtTx([txSigOFDM; zeros(L-1,1)]);
txSig=txSigFOFDM;
end