function bmp_contrast_yuv(filename, contrast, brightness, linear, prefix, showFigures);
%Adjusts global contrast of image, saving output as bitmap with prefix
%  filename: name(s) of bitmap [optional]
%  contrast: intensity contrast 0..1: 0.1 = lower contrast, 0.5= unchanged, 0.9 higher contrast
%  brightness: intensity brightness 0..1: 0.1 = darker, 0.5= unchanged, 0.9  brighter
%        Special value: if brightness is <0, then adaptive so image does not get darker/brighter
%  linear: true/false - select between linear transform or nonlinear
%  prefix: string appended to output name, e.g. if 'e' then cat.jpg ->  ecat.jpg
%  showFigures: Determine whether histogram and function are graphed
%Example
% bmp_contrast; %dialogs will ask the user to select images and settings
% bmp_contrast('dog.png');
% bmp_contrast_yuv('cat.jpg',0.75,-0.5,false,'N95'); %non-linear, auto-brightness
% bmp_contrast_yuv('cat.jpg',0.75,-0.5,true,'L95'); %linear, auto-brightness

%initialize any values that were not explicitly specified
if (nargin < 1)  
   [files,pth] = uigetfile({'*.bmp;*.jpg;*.png;*.tiff;';'*.*'},'Select the Image[s]', 'MultiSelect', 'on'); 
   files = cellstr(files); %make cellstr regardless of whether user selects single or multiple images
else
    [pth,nam, ext] = fileparts(filename);
    files = cellstr([nam, ext]); 
end;
if (nargin < 4)
    prompt = {'Enter contrast (0..1; <0.5 lower, >0.5 higher):','Enter brightness (0..1; <0.5 darker, >0.5 brighter, negative for automatic):','Linear (1) or Nonlinear (0) transform'};
    dlg_title = 'Values for adjusting the image(s)';
    num_lines = 1;
    def = {'0.8','-0.5','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    contrast = str2num(answer{1});
    brightness = str2num(answer{2});
    linear = str2num(answer{3});
end;
if (nargin < 5)  
    prefix = 'c';
end;
if (nargin < 6)
    showFigures = true;
end;

%apply to each of the image(s)
for i=1:size(files,2)
    nam = strvcat(deblank(files(:,i)));
    inName = fullfile(pth, [nam]);
    outName = fullfile(pth, [prefix nam ]);
    Im = imread(inName);
    ImSize = size(Im);
    %scale to range 0..1
    if isa(Im,'uint16')
        scale = 1/65535;
    elseif isa(Im,'uint8')
        scale = 1/255;
    else 
        fprintf('Unsupported data format\n');
        return;
    end;
    Im = double(Im) .* scale;
    %convert RGG to YUV - we will only adjust Y (luminance)
    if length(ImSize) == 3 %3D data has three 2D slices: Red,Green,Blue convert to Y,U,V
        [Y,U,V]=YUV_RGB_sub(Im);        
    else % if not RGB, assume grayscale, where image only stores luminance
        Y = Im;
    end;
    %make a histogram - much more rapid for autobalance, and creates useful graph
    Y1d = reshape(Y,size(Y,1)*size(Y,2),1);
    x = 0:1/255:1; 
    h = hist(Y1d,x);
    %compute mean
    delta=1e-6; % a very small number, so zeros do not cause problems...
    %logMeanY= logMean_sub (h); %you can compute from either histogram or based on each voxel
    logMeanY=exp(mean(mean(log(Y+delta))));%overall luminance: log average of image
    %meanY = mean(mean(Y)) ;
    fprintf(' %s: Before transform: log Mean intensity %f\n', nam, logMeanY);
    %compute transform
    if brightness < 0, %autobalance
        %exhaustively test 256 possible brightness levels
        % find brightness level that preserves original luminance
        %for speed, compute effects on 1D histogram, rather than 2D image
        bestFit = Inf;
        hR = h;
        for i = 0:255,
        	out = makeTransform_sub (contrast, i/256, linear);
            for x = 1:256, %initialize - several input intensities may be set to same output
                hR(x) = 0;
            end; 
            for x = 1:256,
                pos = round(out(x)*255) + 1; 
                hR(pos) = hR(pos) + h(x);
            end; %for x
            logMeanYr= logMean_sub (hR);
            fit = abs(logMeanYr-logMeanY);
            if (fit < bestFit) 
                brightness = i/256;
                bestFit = fit;
            end;
        end;
        %fprintf(' Auotbalance set the to brightness %f\n',brightness);
    end; %if brightness < 0: autobalance         
     out = makeTransform_sub (contrast, brightness, linear);
     %apply transform
     for y = 1:ImSize(2),
         for x = 1:ImSize(1),
                  Y(x,y) = out(1+round(255*Y(x,y)));
         end;
     end;
     %report effects
     logMeanY=exp(mean(mean(log(Y+delta))));%overall luminance: log average of image
     fprintf(' %s: After transform: log Mean intensity %f  with contrast=%.2f and brightness=%.2f\n', [prefix nam ], logMeanY, contrast, brightness);
     if showFigures %graphic display of adjustment
        subplot(1,3,2);
        plotTransform_sub (out, ['Transform for ' nam sprintf(' with contrast=%.2f and brightness=%.2f', contrast, brightness )]);        
        hOrig = hist(Y1d,x);
        Y1dT = reshape(Y,size(Y,1)*size(Y,2),1);
        hTrans = hist(Y1dT,x);
        if max(hOrig) > max(hTrans)
            histMax = max(hOrig);
        else
            histMax = max(hTrans);
        end;
        subplot(1,3,1);
        hist(Y1d,x);
        axis([0 1 0 histMax]);
        title('Original histogram');
        subplot(1,3,3);
        hist(Y1dT,x);
        axis([0 1 0 histMax]);
        title('Transformed histogram');
     end; %if showFigures
     %save data
     if length(ImSize) == 3 %RGB
        imAdjusted = RGB_YUV_sub(Y,U,V);
        imwrite(imAdjusted,outName);
     else %not RGB - therefore Grayscale
        imwrite(Y,outName); 
     end;
end; %for each image

function [out]= logMean_sub (histo8bit);
delta=1e-6; % a very small number, so zeros do not cause problems...
%logMeanY=exp(mean(mean(log(Y+delta))));   
n = sum(histo8bit); %number of pixels
logV = 0.0;
for i = 1:256
    logV = logV + (histo8bit(i)*log( ((i-1)/255)+delta)); 
end;
out = exp(logV/n);
%end subfunction logMean_sub

function plotTransform_sub (out, label);
%OPTIONAL: plot lookup table
in = 0:1/255:1;
%figure; %to save this image rather than overwrite
p =plot(in, out);
axis([0 1 0 1]);
xlabel('Input Intensity');
ylabel('Output intensity');
title( label);
set(p,'LineWidth',2);
set(gcf,'Color',[1 1 1]);
%end subfunction plotTransform_sub

function [out]= makeTransform_sub (contrast, brightness, linear);
%Create lookup table for chosen contrast/brightness
if brightness < 0.00001
    brightness = 0.00001;
end;
if contrast < 0.00001
    contrast = 0.00001;
end;
in = 0:1/255:1;
out = 0:1/255:1;
if linear
    midpoint = brightness;
    if contrast == 1.0
        slope = 1;
    else
        deg = contrast * 90; %0= no contrast = horizontal line, 1=no gray = vertical line
        rad = deg*pi/180; %degrees to radians
        slope = tan(rad); %http://en.wikipedia.org/wiki/Slope
    end;
    for i = 1:256,
        v = (((in(i) -midpoint)* slope)+ 0.5);
        
        if v > 1
            v = 1;
        elseif v < 0
                v = 0;
        end;
        out(i) = v;    
        %out(i) = in(i) * 2;
    end; %for i: each possible intensity
else %if not linear than non-linear
    %http://dept-info.labri.fr/~schlick/DOC/gem2.html
    %http://dept-info.labri.fr/~schlick/publi.html
    % "Fast Alternatives to Perlin's Bias and Gain Functions"
    %   Christophe Schlick: Graphics Gems IV, p379-382, April 1994 
    for i = 1:256,
        lT = in(i);
        %apply bias (brightness)
        lT = (lT/((1/brightness-2)*(1-lT)+1)) ;
        %next apply gain (contrast)
        gain = 1.0-contrast;
        if lT <= 0.5 
            lG = (lT/((1/gain-2)*(1-2*lT)+1));
        else
            lG = (( (1/gain-2)*(1-2*lT)-lT ) / ( (1/gain-2)*(1-2*lT)-1 ) );
        end;
        if lT == 0
            lG = 0;
        else
            lG = lG / lT;
        end;
        lV = (lT*lG);
        if lV > 1
            lV = 1;
        end;
        if lV < 0
            lV = 0;
        end;
        out(i) = lV;
        
    end; %for i: each possible intensity
end; %nonlinear 
%end subfunction makeTransform_sub

%I originally used RGB<->YUV by Mohammed Mustafa Siddeq:
% http://www.mathworks.com/matlabcentral/fileexchange/33352
%However, I now include a more accurate conversion 

function [Y,U,V]=YUV_RGB_sub(Im)
%http://en.wikipedia.org/wiki/YUV
    Im=double(Im);
    R=Im(:,:,1); G=Im(:,:,2); B=Im(:,:,3);
    % transfom layers to YUV
    Y = 0.299*R + 0.587*G + 0.114*B;
    U=(B-Y)*0.492;
    V=(R-Y)*0.877;
% end YUV_RGB_sub

function Im=RGB_YUV_sub(Y,U,V)
% http://en.wikipedia.org/wiki/YUV
R = Y + 1.13983*V;  %1.1403*V = V/0.877
G = Y - 0.39465*U - 0.58060*V;
B = Y + 2.03211*U; %2.0325*U = U/0.492
Im(:,:,1)=R; Im(:,:,2)=G; Im(:,:,3)=B; 
