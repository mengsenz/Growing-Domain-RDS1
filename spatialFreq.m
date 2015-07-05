function [ freq ] = spatialFreq( concentration, L)
%SPATIALFREQ compute the spatial frequency of the strongest wave mode
%(cycle/unit distance) for 1-D data (column-wise)
%   usage: [ freq ] = spatialFreq( concentration, L)
%   input:
%           concentration:  u or v
%           L:              length of the domain.

% --- check input
tResolution = size(concentration,2);
concentration = [flipud(concentration);concentration]; 
nx = size(concentration,1);
nfft=nx;
% --- Fourier Transform for reduced time points
fIdx = 1:nfft/2+1;
fuNew=fft(concentration,nfft,1);
fuNew(1,:)=0;
fuNew=abs(fuNew(fIdx,:)); 
% --- Find highest frequency peak
freq = zeros(1,tResolution);
for i=1:tResolution
    [peak, locs] = findpeaks(fuNew(:,i));
    if isempty(peak), 
        freq(i)=0;
    else
        [~, maxpeakIdx] = max(peak);
        freq(i) = (locs(maxpeakIdx)-1);
    end
end
freq=freq/2/nfft*nx./L;
end

