function txSig=ufmc_mod(data)
load initpara

%  Loop over each subband
for bandIdx = 1:numSubbands

    bitsIn = data(:,bandIdx);
    symbolsIn = qamMapper( bitsIn);
    inpData(:,bandIdx) = bitsIn; % log bits for comparison
        
    % Pack subband data into an OFDM symbol
    offset = subbandOffset+(bandIdx-1)*subbandSize; 
    symbolsInOFDM = [zeros(offset,1); symbolsIn; ...
                     zeros(numFFT-offset-subbandSize, 1)];
    ifftOut = ifft(ifftshift(symbolsInOFDM));
    
    % Filter for each subband is shifted in frequency
    bandFilter = prototypeFilter.*exp( 1i*2*pi*(0:filterLen-1)'/numFFT* ...
                 ((bandIdx-1/2)*subbandSize+0.5+subbandOffset+numFFT/2) );    
    filterOut = conv(bandFilter,ifftOut);
    
    % Plot power spectral density (PSD) per subband
    [psd,f] = periodogram(filterOut, rectwin(length(filterOut)), ...
                          numFFT*2, 1, 'centered'); 
    plot(f,10*log10(psd)); 
    
    % Sum the filtered subband responses to form the aggregate transmit signal
    txSig = txSig + filterOut;     
end

