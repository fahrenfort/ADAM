function outIm = phaseScramble(im, perc, newSize)
% function outIm = phaseScramble(im, perc, newSize)
%
% Scramble phase of image im, leaving color space intact
% im is phase scrambled by percentage perc(can be an array, in which case output is cell array)
% if perc = 0, no phase information is retained
% the output image is resized to newSize*newSize
%
% By J.J.Fahrenfort, UvA 2010, based on a script by Sennay Ghebreab

if nargin<2
    perc=1;
end

% determine plot dimensions
dim1 = 3;
dim2 = 4;
%figure;

% Generate random phase structure for this image
RandomPhase = angle(fft2(rand(size(im,1), size(im,2)))); 

% for each color layer
for layer = 1:size(im,3)
    % Obtain fast-Fourier transform
    ImFourier(:,:,layer) = fft2(im(:,:,layer));
    % Obtain amplitude spectrum
    Amp(:,:,layer) = abs(ImFourier(:,:,layer));
    % Obtain phase spectrum
    Phase(:,:,layer) = angle(ImFourier(:,:,layer));
end

% Phase scramble images
count = 0;
for cPerc = 1:numel(perc)
    clear ImScrambled;
    % for each color layer
    for layer = 1:size(im,3)
        % Scramble phase
        %if perc(cPerc) == 1
        %    NewPhase(:,:,layer) = RandomPhase; % total annihilation
        %else
            NewPhase(:,:,layer) = Phase(:,:,layer) + RandomPhase.*perc(cPerc); % add random phase to original phase
        %end
        % Combine Amp and Phase then perform inverse Fourier
        ImScrambled(:,:,layer) = ifft2(Amp(:,:,layer).*exp(sqrt(-1)*(NewPhase(:,:,layer))));         
    end
    % Get rid of imaginery
    ImScrambled = abs(real(ImScrambled));
    % Scale from 0 to 1
    outIm{cPerc} = ImScrambled./max(ImScrambled(:));
    
    % Resize and display images
    [ oldSize1 oldSize2 dummy ] = size(outIm{cPerc});
    if all([newSize newSize] ~= [oldSize1 oldSize2])
        outIm{cPerc} = imresize(outIm{cPerc},[newSize newSize]);
        disp('resized');
    end
    if mod(cPerc,2)==0 || cPerc==numel(perc)
        count = count + 1;
        %subplot(dim1,dim2,count);
        %imshow(onlyCircle(outIm{cPerc}));
    end
end

% If only one image, remove from cell
if cPerc == 1
    outIm = outIm{1};
%else
%    subplot(dim1,dim2,count+1);
end

% make fullscreen
%set(gcf, 'Position', get(0,'Screensize'));