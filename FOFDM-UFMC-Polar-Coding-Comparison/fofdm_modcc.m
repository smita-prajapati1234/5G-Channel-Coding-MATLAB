function txSig=fofdm_modcc(data)
load initpara

data12=data(:);
cdata=convenc(data12,poly2trellis(ConstraintLength,CodeGenerator));
cdata1 = randintrlv(cdata,1008);

% QAM Symbol mapper
symbolsIn = qammod(cdata1, 2^bitsPerSubCarrier, 'InputType', 'bit', ...
    'UnitAveragePower', true);

% Pack data into an OFDM symbol
offset = (numFFT-numDataCarriers)/2; % for band center
symbolsInOFDM = [zeros(offset,1); symbolsIn; ...
    zeros(numFFT-offset-numDataCarriers,1)];
ifftOut = ifft(ifftshift(symbolsInOFDM));

% Prepend cyclic prefix
txSigOFDM = [ifftOut(end-cpLen+1:end); ifftOut];

% Filter, with zero-padding to flush tail. Get the transmit signal
txSigFOFDM = filtTx([txSigOFDM; zeros(filterLen-1,1)]);
txSig=txSigFOFDM;
end