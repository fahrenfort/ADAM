kLowContrastSteps = 5;
kHighContrastSteps = 5;
kChangePerStep = 0.05; %e.g. 0.1 will have contrast steps 0.4, 0.5, 0.6, 0.7...
%have user select files
[files,pth] = uigetfile({'*.bmp;*.jpg;*.png;*.tiff;';'*.*'},'Select the Image[s]', 'MultiSelect', 'on'); 
files = cellstr(files); %make cellstr regardless of whether user selects single or multiple images

%for s = -kLowContrastSteps : kHighContrastSteps
    %step = 0.5 + (kChangePerStep * s);
    step = .25;
    label = ['g', int2str(round( step*100) )];
    for f=1:size(files,2)
        nam = fullfile(pth, strvcat(deblank(files(:,f))) );
        fprintf('file %s gain %.2f label %s\n',nam, step, label);
        if step < 0.5 %if reducing contrast, use linear transform
            bmp_contrast(nam,step,-0.5,true,label,false); 
        else %if increasing contrast, use non-linear transform
            bmp_contrast(nam,step,-0.5,false,label, false);
        end;
    end;
%end;