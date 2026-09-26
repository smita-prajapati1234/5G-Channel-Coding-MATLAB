function txSig=ufmc_modcc(data)
load initpara

%  Loop over each subband
for bandIdx = 1:numSubbands

    bitsIn = data(:,bandIdx);
%      bitsIn = randi([0 1], bitsPerSubCarrier*subbandSize*coderate, 1);
    cdata=convenc(bitsIn,poly2trellis(ConstraintLength,CodeGenerator));
    cdata1 = randintrlv(cdata,1008);

    symbolsIn = qamMapper( cdata1);
    inpData1(:,bandIdx) = cdata1; % log bits for comparison
        
    % Pack subband data into an OFDM symbol
    offset = subbandOffset+(bandIdx-1)*subbandSize; 
    symbolsInOFDM = [zeros(offset,1); symbolsIn; ...
                     zeros(numFFT-offset-subbandSize, 1)];
    ifftOut = ifft(ifftshift(symbolsInOFDM));
    
    % Filter for each subband is shifted in frequency
    bandFilter = prototypeFilter.*exp( 1i*2*pi*(0:filterLen-1)'/numFFT* ...
                 ((bandIdx-1/2)*subbandSize+0.5+subbandOffset+numFFT/2) );    
    filterOut = conv(bandFilter,ifftOut);
    
     
    % Sum the filtered subband responses to form the aggregate transmit signal
    txSig = txSig + filterOut;     
end
PAPR = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
[~,~,paprUFMC] = PAPR(txSig);
disp(['Peak-to-Average-Power-Ratio (PAPR) for UFMC = ' num2str(paprUFMC) ' dB']);


