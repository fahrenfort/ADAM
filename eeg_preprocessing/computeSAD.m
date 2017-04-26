
% computeSAD() - Computes Spatial Average Difference feature 
%
% Usage:
%   >> [rapp,var_front,var_back,mean_front,mean_back]=computeSAD(topog,chanlocs,n);
%
% Inputs:
%   topog      - topographies vector
%   chanlocs   - EEG.chanlocs struct
%   n          - number of ICs
%   nchannels  - number of channels
%
% Outputs:
%   rapp       - SAD values
%   var_front  - Frontal Area variance values
%   var_back   - Posterior Area variance values
%   mean_front - Frontal Area average values
%   mean_back  - Posterior Area average values
%
%
% Copyright (C) 2009-2014 Andrea Mognon (1) and Marco Buiatti (2), 
% (1) Center for Mind/Brain Sciences, University of Trento, Italy
% (2) INSERM U992 - Cognitive Neuroimaging Unit, Gif sur Yvette, France
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


function [rapp,var_front,var_back,mean_front,mean_back]=computeSAD(topog,chanlocs,n)

nchannels=length(chanlocs);

%% Define scalp zones

% Find electrodes in Frontal Area (FA)
dimfront=0; %number of FA electrodes
index1=zeros(1,nchannels); %indexes of FA electrodes

for k=1:nchannels
    if (abs(chanlocs(1,k).theta)<60) && (chanlocs(1,k).radius>0.40) %electrodes are in FA
        dimfront=dimfront+1; %count electrodes
        index1(1,dimfront)=k; 
    end
end

if sum(index1) == 0
    disp('WARNING: relaxing criterion for frontal electrodes in ADJUST');
    for k=1:nchannels
        if (abs(chanlocs(1,k).theta)<60) && (chanlocs(1,k).radius>0.30) %electrodes are in FA
            dimfront=dimfront+1; %count electrodes
            index1(1,dimfront)=k;
        end
    end
end

 % Find electrodes in Posterior Area (PA)
    dimback=0;
    index3=zeros(1,nchannels);
    for h=1:nchannels 
        if (abs(chanlocs(1,h).theta)>110) 
            dimback=dimback+1; 
            index3(1,dimback)=h; 
        end
    end
 
    if dimfront*dimback==0
        disp('ERROR: no channels included in some scalp areas.')
        disp('Check channels distribution and/or change scalp areas definitions in computeSAD.m and computeSED_NOnorm.m')
        disp('ADJUST session aborted.')
        return
    end
    
%% Outputs

     rapp=zeros(1,n); % SAD
      mean_front=zeros(1,n); % FA electrodes mean value
      mean_back=zeros(1,n); % PA electrodes mean value
      var_front=zeros(1,n); % FA electrodes variance value
      var_back=zeros(1,n); % PA electrodes variance value

%% Output computation

for i=1:n % for each topography
    
 %create FA electrodes vector
    front=zeros(1,dimfront);
    for h=1:dimfront
        front(1,h)=topog(i,index1(1,h));
    end
    
  %create PA electrodes vector
    back=zeros(1,dimback);
    for h=1:dimback
        back(1,h)=topog(i,index3(1,h));
    end
    
     
    
   %compute features
    
    rapp(1,i)=abs(mean(front))-abs(mean(back)); % SAD
    mean_front(1,i)=mean(front);
    mean_back(1,i)=mean(back);
    var_back(1,i)=var(back);
    var_front(1,i)=var(front);
    
end

