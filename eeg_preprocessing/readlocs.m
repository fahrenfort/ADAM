% readlocs() - read electrode location coordinates and other information from a file. 
%              Several standard file formats are supported. Users may also specify 
%              a custom column format. Defined format examples are given below 
%              (see File Formats).
% Usage:
%   >>  eloc = readlocs( filename );
%   >>  EEG.chanlocs = readlocs( filename, 'key', 'val', ... ); 
%   >>  [eloc, labels, theta, radius, indices] = ...
%                                               readlocs( filename, 'key', 'val', ... );
% Inputs:
%   filename   - Name of the file containing the electrode locations
%                {default: 2-D polar coordinates} (see >> help topoplot )
%
% Optional inputs:
%   'filetype'  - ['loc'|'sph'|'sfp'|'xyz'|'asc'|'polhemus'|'besa'|'chanedit'|'custom'] 
%                 Type of the file to read. By default the file type is determined 
%                 using the file extension (see below under File Formats),
%                  'loc'   an EEGLAB 2-D polar coordinates channel locations file 
%                          Coordinates are theta and radius (see definitions below).
%                  'sph'   Matlab spherical coordinates (Note: spherical
%                          coordinates used by Matlab functions are different 
%                          from spherical coordinates used by BESA - see below).
%                  'sfp'   EGI Cartesian coordinates (NOT Matlab Cartesian - see below).
%                  'xyz'   Matlab/EEGLAB Cartesian coordinates (NOT EGI Cartesian).
%                          z is toward nose; y is toward left ear; z is toward vertex
%                  'asc'   Neuroscan polar coordinates.
%                  'polhemus' or 'polhemusx' - Polhemus electrode location file recorded 
%                          with 'X' on sensor pointing to subject (see below and readelp()).
%                  'polhemusy' - Polhemus electrode location file recorded with 
%                          'Y' on sensor pointing to subject (see below and readelp()).
%                  'besa' BESA-'.elp' spherical coordinates. (Not MATLAB spherical -
%                           see below).
%                  'chanedit' - EEGLAB channel location file created by pop_chanedit().
%                  'custom' - Ascii file with columns in user-defined 'format' (see below).
%   'importmode' - ['eeglab'|'native'] for location files containing 3-D cartesian electrode
%                  coordinates, import either in EEGLAB format (nose pointing toward +X). 
%                  This may not always be possible since EEGLAB might not be able to 
%                  determine the nose direction for scanned electrode files. 'native' import
%                  original carthesian coordinates (user can then specify the position of
%                  the nose when calling the topoplot() function; in EEGLAB the position
%                  of the nose is stored in the EEG.chaninfo structure). {default 'eeglab'}
%   'format'    -  [cell array] Format of a 'custom' channel location file (see above).
%                  {default: if no file type is defined. The cell array contains
%                  labels defining the meaning of each column of the input file.
%                           'channum'   [positive integer] channel number.
%                           'labels'    [string] channel name (no spaces).
%                           'theta'     [real degrees] 2-D angle in polar coordinates.
%                                       positive => rotating from nose (0) toward left ear
%                           'radius'    [real] radius for 2-D polar coords; 0.5 is the head
%                                       disk radius and limit for topoplot() plotting).
%                           'X'         [real] Matlab-Cartesian X coordinate (to nose).
%                           'Y'         [real] Matlab-Cartesian Y coordinate (to left ear).
%                           'Z'         [real] Matlab-Cartesian Z coordinate (to vertex).
%                           '-X','-Y','-Z' Matlab-Cartesian coordinates pointing opposite
%                                       to the above.
%                           'sph_theta' [real degrees] Matlab spherical horizontal angle.
%                                       positive => rotating from nose (0) toward left ear.
%                           'sph_phi'   [real degrees] Matlab spherical elevation angle.
%                                       positive => rotating from horizontal (0) upwards.
%                           'sph_radius' [real] distance from head center (unused).
%                           'sph_phi_besa' [real degrees] BESA phi angle from vertical.
%                                       positive => rotating from vertex (0) towards right ear.
%                           'sph_theta_besa' [real degrees] BESA theta horiz/azimuthal angle.
%                                       positive => rotating from right ear (0) toward nose.
%                           'ignore'    ignore column}.
%     The input file may also contain other channel information fields.
%                           'type'      channel type: 'EEG', 'MEG', 'EMG', 'ECG', others ...
%                           'calib'     [real near 1.0] channel calibration value.
%                           'gain'      [real > 1] channel gain.
%                           'custom1'   custom field #1.
%                           'custom2', 'custom3', 'custom4', etc.    more custom fields
%   'skiplines' - [integer] Number of header lines to skip (in 'custom' file types only).
%                 Note: Characters on a line following '%' will be treated as comments.
%   'readchans' - [integer array] indices of electrodes to read. {default: all}
%   'center'    - [(1,3) real array or 'auto'] center of xyz coordinates for conversion 
%                 to spherical or polar, Specify the center of the sphere here, or 'auto'. 
%                 This uses the center of the sphere that best fits all the electrode 
%                 locations read. {default: [0 0 0]}
% Outputs:
%   eloc        - structure containing the channel names and locations (if present).
%                 It has three fields: 'eloc.labels', 'eloc.theta' and 'eloc.radius' 
%                 identical in meaning to the EEGLAB struct 'EEG.chanlocs'.
%   labels      - cell array of strings giving the names of the electrodes. NOTE: Unlike the
%                 three outputs below, includes labels of channels *without* location info.
%   theta       - vector (in degrees) of polar angles of the electrode locations.
%   radius      - vector of polar-coordinate radii (arc_lengths) of the electrode locations 
%   indices     - indices, k, of channels with non-empty 'locs(k).theta' coordinate
%
% File formats:
%   If 'filetype' is unspecified, the file extension determines its type.
%
%   '.loc' or '.locs' or '.eloc': 
%               polar coordinates. Notes: angles in degrees: 
%               right ear is 90; left ear -90; head disk radius is 0.5. 
%               Fields:   N    angle  radius    label
%               Sample:   1    -18    .511       Fp1   
%                         2     18    .511       Fp2  
%                         3    -90    .256       C3
%                         4     90    .256       C4
%                           ...
%               Note: In previous releases, channel labels had to contain exactly 
%               four characters (spaces replaced by '.'). This format still works, 
%               though dots are no longer required.
%   '.sph':
%               Matlab spherical coordinates. Notes: theta is the azimuthal/horizontal angle
%               in deg.: 0 is toward nose, 90 rotated to left ear. Following this, performs
%               the elevation (phi). Angles in degrees.
%               Fields:   N    theta    phi    label
%               Sample:   1      18     -2      Fp1
%                         2     -18     -2      Fp2
%                         3      90     44      C3
%                         4     -90     44      C4
%                           ...
%   '.elc':
%               Cartesian 3-D electrode coordinates scanned using the EETrak software. 
%               See readeetraklocs().
%   '.elp':     
%               Polhemus-.'elp' Cartesian coordinates. By default, an .elp extension is read
%               as PolhemusX-elp in which 'X' on the Polhemus sensor is pointed toward the 
%               subject. Polhemus files are not in columnar format (see readelp()).
%   '.elp':
%               BESA-'.elp' spherical coordinates: Need to specify 'filetype','besa'.
%               The elevation angle (phi) is measured from the vertical axis. Positive 
%               rotation is toward right ear. Next, perform azimuthal/horizontal rotation 
%               (theta): 0 is toward right ear; 90 is toward nose, -90 toward occiput. 
%               Angles are in degrees.  If labels are absent or weights are given in 
%               a last column, readlocs() adjusts for this. Default labels are E1, E2, ...
%               Fields:   Type  label      phi  theta   
%               Sample:   EEG   Fp1        -92   -72    
%                         EEG   Fp2         92    72   
%                         EEG   C3         -46    0  
%                         EEG   C4          46    0 
%                           ...
%   '.xyz': 
%               Matlab/EEGLAB Cartesian coordinates. Here. x is towards the nose, 
%               y is towards the left ear, and z towards the vertex. Note that the first
%               column (x) is -Y in a Matlab 3-D plot, the second column (y) is X in a 
%               matlab 3-D plot, and the third column (z) is Z.
%               Fields:   channum   x           y         z     label
%               Sample:   1       .950        .308     -.035     Fp1
%                         2       .950       -.308     -.035     Fp2
%                         3        0           .719      .695    C3
%                         4        0          -.719      .695    C4
%                           ...
%   '.asc', '.dat':     
%               Neuroscan-.'asc' or '.dat' Cartesian polar coordinates text file.
%   '.sfp': 
%               BESA/EGI-xyz Cartesian coordinates. Notes: For EGI, x is toward right ear, 
%               y is toward the nose, z is toward the vertex. EEGLAB converts EGI 
%               Cartesian coordinates to Matlab/EEGLAB xyz coordinates. 
%               Fields:   label   x           y          z
%               Sample:   Fp1    -.308        .950      -.035    
%                         Fp2     .308        .950      -.035  
%                         C3     -.719        0          .695  
%                         C4      .719        0          .695  
%                           ...
%   '.ced':   
%               ASCII file saved by pop_chanedit(). Contains multiple MATLAB/EEGLAB formats.
%               Cartesian coordinates are as in the 'xyz' format (above).
%               Fields:   channum  label  theta  radius   x      y      z    sph_theta   sph_phi  ...
%               Sample:   1        Fp1     -18    .511   .950   .308  -.035   18         -2       ...
%                         2        Fp2      18    .511   .950  -.308  -.035  -18         -2       ...
%                         3        C3      -90    .256   0      .719   .695   90         44       ...
%                         4        C4       90    .256   0     -.719   .695  -90         44       ...
%                           ...
%               The last columns of the file may contain any other defined fields (gain,
%               calib, type, custom).
%
%    Fieldtrip structure: 
%               If a Fieltrip structure is given as input, an EEGLAB
%               chanlocs structure is returned
%
% Author: Arnaud Delorme, Salk Institute, 8 Dec 2002
%
% See also: readelp(), writelocs(), topo2sph(), sph2topo(), sph2cart()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 28 Feb 2002
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


function [eloc, labels, theta, radius, indices] = readlocs( filename, varargin ); 

if nargin < 1
	help readlocs;
	return;
end;

% NOTE: To add a new channel format:
% ----------------------------------
% 1) Add a new element to the structure 'chanformat' (see 'ADD NEW FORMATS HERE' below):
% 2)  Enter a format 'type' for the new file format, 
% 3)  Enter a (short) 'typestring' description of the format
% 4)  Enter a longer format 'description' (possibly multiline, see ex. (1) below)
% 5)  Enter format file column labels in the 'importformat' field (see ex. (2) below)
% 6)  Enter the number of header lines to skip (if any) in the 'skipline' field
% 7)  Document the new channel format in the help message above.
% 8)  After testing, please send the new version of readloca.m to us
%       at eeglab@sccn.ucsd.edu with a sample locs file.
% The 'chanformat' structure is also used (automatically) by the writelocs() 
% and pop_readlocs() functions. You do not need to edit these functions.

chanformat(1).type         = 'polhemus';
chanformat(1).typestring   = 'Polhemus native .elp file';
chanformat(1).description  = [ 'Polhemus native coordinate file containing scanned electrode positions. ' ...
                               'User must select the direction ' ...
                               'for the nose after importing the data file.' ];
chanformat(1).importformat = 'readelp() function';
% ---------------------------------------------------------------------------------------------------
chanformat(2).type         = 'besa';
chanformat(2).typestring   = 'BESA spherical .elp file';
chanformat(2).description  = [ 'BESA spherical coordinate file. Note that BESA spherical coordinates ' ...
                               'are different from Matlab spherical coordinates' ];
chanformat(2).skipline     = 0; % some BESA files do not have headers
chanformat(2).importformat = { 'type' 'labels' 'sph_theta_besa' 'sph_phi_besa' 'sph_radius' };
% ---------------------------------------------------------------------------------------------------
chanformat(3).type         = 'xyz';
chanformat(3).typestring   = 'Matlab .xyz file';
chanformat(3).description  = [ 'Standard 3-D cartesian coordinate files with electrode labels in ' ...
                               'the first column and X, Y, and Z coordinates in columns 2, 3, and 4' ];
chanformat(3).importformat = { 'channum' '-Y' 'X' 'Z' 'labels'};
% ---------------------------------------------------------------------------------------------------
chanformat(4).type         = 'sfp';
chanformat(4).typestring   = 'BESA or EGI 3-D cartesian .sfp file';
chanformat(4).description  = [ 'Standard BESA 3-D cartesian coordinate files with electrode labels in ' ...
                               'the first column and X, Y, and Z coordinates in columns 2, 3, and 4.' ...
                               'Coordinates are re-oriented to fit the EEGLAB standard of having the ' ...
                               'nose along the +X axis.' ];
chanformat(4).importformat = { 'labels' '-Y' 'X' 'Z' };
chanformat(4).skipline     = 0;
% ---------------------------------------------------------------------------------------------------
chanformat(5).type         = 'loc';
chanformat(5).typestring   = 'EEGLAB polar .loc file';
chanformat(5).description  = [ 'EEGLAB polar .loc file' ];
chanformat(5).importformat = { 'channum' 'theta' 'radius' 'labels' };
% ---------------------------------------------------------------------------------------------------
chanformat(6).type         = 'sph';
chanformat(6).typestring   = 'Matlab .sph spherical file';
chanformat(6).description  = [ 'Standard 3-D spherical coordinate files in Matlab format' ];
chanformat(6).importformat = { 'channum' 'sph_theta' 'sph_phi' 'labels' };
% ---------------------------------------------------------------------------------------------------
chanformat(7).type         = 'asc';
chanformat(7).typestring   = 'Neuroscan polar .asc file';
chanformat(7).description  = [ 'Neuroscan polar .asc file, automatically recentered to fit EEGLAB standard' ...
                               'of having ''Cz'' at (0,0).' ];
chanformat(7).importformat = 'readneurolocs';
% ---------------------------------------------------------------------------------------------------
chanformat(8).type         = 'dat';
chanformat(8).typestring   = 'Neuroscan 3-D .dat file';
chanformat(8).description  = [ 'Neuroscan 3-D cartesian .dat file. Coordinates are re-oriented to fit ' ...
                               'the EEGLAB standard of having the nose along the +X axis.' ];
chanformat(8).importformat = 'readneurolocs';
% ---------------------------------------------------------------------------------------------------
chanformat(9).type         = 'elc';
chanformat(9).typestring   = 'ASA .elc 3-D file';
chanformat(9).description  = [ 'ASA .elc 3-D coordinate file containing scanned electrode positions. ' ...
                               'User must select the direction ' ...
                               'for the nose after importing the data file.' ];
chanformat(9).importformat = 'readeetraklocs';
% ---------------------------------------------------------------------------------------------------
chanformat(10).type         = 'chanedit';
chanformat(10).typestring   = 'EEGLAB complete 3-D file';
chanformat(10).description  = [ 'EEGLAB file containing polar, cartesian 3-D, and spherical 3-D ' ...
                               'electrode locations.' ];
chanformat(10).importformat = { 'channum' 'labels'  'theta' 'radius' 'X' 'Y' 'Z' 'sph_theta' 'sph_phi' ...
                               'sph_radius' 'type' };
chanformat(10).skipline     = 1;
% ---------------------------------------------------------------------------------------------------
chanformat(11).type         = 'custom';
chanformat(11).typestring   = 'Custom file format';
chanformat(11).description  = 'Custom ASCII file format where user can define content for each file columns.';
chanformat(11).importformat = '';
% ---------------------------------------------------------------------------------------------------
% ----- ADD MORE FORMATS HERE -----------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------

listcolformat = { 'labels' 'channum' 'theta' 'radius' 'sph_theta' 'sph_phi' ...
      'sph_radius' 'sph_theta_besa' 'sph_phi_besa' 'gain' 'calib' 'type' ...
      'X' 'Y' 'Z' '-X' '-Y' '-Z' 'custom1' 'custom2' 'custom3' 'custom4' 'ignore' 'not def' };

% ----------------------------------
% special mode for getting the info
% ----------------------------------
if isstr(filename) & strcmp(filename, 'getinfos')
   eloc = chanformat;
   labels = listcolformat;
   return;
end;

g = finputcheck( varargin, ...
   { 'filetype'	   'string'  {}                 '';
     'importmode'  'string'  { 'eeglab','native' } 'eeglab';
     'defaultelp'  'string'  { 'besa','polhemus' } 'polhemus';
     'skiplines'   'integer' [0 Inf] 			[];
     'elecind'     'integer' [1 Inf]	    	[];
     'format'	   'cell'	 []					{} }, 'readlocs');
if isstr(g), error(g); end;  
if ~isempty(g.format), g.filetype = 'custom'; end;

if isstr(filename)
   
   % format auto detection
	% --------------------
   if strcmpi(g.filetype, 'autodetect'), g.filetype = ''; end;
   g.filetype = strtok(g.filetype);
   periods = find(filename == '.');
   fileextension = filename(periods(end)+1:end);
   g.filetype = lower(g.filetype);
   if isempty(g.filetype)
       switch lower(fileextension),
        case {'loc' 'locs' 'eloc'}, g.filetype = 'loc'; % 5/27/2014 Ramon: 'eloc' option introduced.
        case 'xyz', g.filetype = 'xyz'; 
          fprintf( [ 'WARNING: Matlab Cartesian coord. file extension (".xyz") detected.\n' ... 
                  'If importing EGI Cartesian coords, force type "sfp" instead.\n'] );
        case 'sph', g.filetype = 'sph';
        case 'ced', g.filetype = 'chanedit';
        case 'elp', g.filetype = g.defaultelp;
        case 'asc', g.filetype = 'asc';
        case 'dat', g.filetype = 'dat';
        case 'elc', g.filetype = 'elc';
        case 'eps', g.filetype = 'besa';
        case 'sfp', g.filetype = 'sfp';
        otherwise, g.filetype =  ''; 
       end;
       fprintf('readlocs(): ''%s'' format assumed from file extension\n', g.filetype); 
   else 
       if strcmpi(g.filetype, 'locs'),  g.filetype = 'loc'; end
       if strcmpi(g.filetype, 'eloc'),  g.filetype = 'loc'; end
   end;
   
   % assign format from filetype
   % ---------------------------
   if ~isempty(g.filetype) & ~strcmpi(g.filetype, 'custom') ...
           & ~strcmpi(g.filetype, 'asc') & ~strcmpi(g.filetype, 'elc') & ~strcmpi(g.filetype, 'dat')
      indexformat = strmatch(lower(g.filetype), { chanformat.type }, 'exact');
      g.format = chanformat(indexformat).importformat;
      if isempty(g.skiplines)
         g.skiplines = chanformat(indexformat).skipline;
      end;
      if isempty(g.filetype) 
         error( ['readlocs() error: The filetype cannot be detected from the \n' ...
                 '                  file extension, and custom format not specified']);
      end;
   end;
   
   % import file
   % -----------
   if strcmp(g.filetype, 'asc') | strcmp(g.filetype, 'dat')
       eloc = readneurolocs( filename );
       eloc = rmfield(eloc, 'sph_theta'); % for the conversion below
       eloc = rmfield(eloc, 'sph_theta_besa'); % for the conversion below
       if isfield(eloc, 'type')
           for index = 1:length(eloc)
               type = eloc(index).type;
               if type == 69,     eloc(index).type = 'EEG';
               elseif type == 88, eloc(index).type = 'REF';
               elseif type >= 76 & type <= 82, eloc(index).type = 'FID';
               else eloc(index).type = num2str(eloc(index).type);
               end;
           end;
       end;
   elseif strcmp(g.filetype, 'elc')
       eloc = readeetraklocs( filename );
       %eloc = read_asa_elc( filename ); % from fieldtrip
       %eloc = struct('labels', eloc.label, 'X', mattocell(eloc.pnt(:,1)'), 'Y', ...
       %                        mattocell(eloc.pnt(:,2)'), 'Z', mattocell(eloc.pnt(:,3)'));
       eloc = convertlocs(eloc, 'cart2all');
       eloc = rmfield(eloc, 'sph_theta'); % for the conversion below
       eloc = rmfield(eloc, 'sph_theta_besa'); % for the conversion below
   elseif strcmp(lower(g.filetype(1:end-1)), 'polhemus') | ...
           strcmp(g.filetype, 'polhemus')
       try, 
           [eloc labels X Y Z]= readelp( filename );
           if strcmp(g.filetype, 'polhemusy')
               tmp = X; X = Y; Y = tmp;
           end;
           for index = 1:length( eloc )
               eloc(index).X = X(index);
               eloc(index).Y = Y(index);	
               eloc(index).Z = Z(index);	
           end;
       catch, 
           disp('readlocs(): Could not read Polhemus coords. Trying to read BESA .elp file.');
           [eloc, labels, theta, radius, indices] = readlocs( filename, 'defaultelp', 'besa', varargin{:} );
       end;
   else      
       % importing file
       % --------------
       if isempty(g.skiplines), g.skiplines = 0; end;
       if strcmpi(g.filetype, 'chanedit')
           array = loadtxt( filename, 'delim', 9, 'skipline', g.skiplines, 'blankcell', 'off');
       else
           array = load_file_or_array( filename, g.skiplines);
       end;
       if size(array,2) < length(g.format)
           fprintf(['readlocs() warning: Fewer columns in the input than expected.\n' ...
                    '                    See >> help readlocs\n']);
       elseif size(array,2) > length(g.format)
           fprintf(['readlocs() warning: More columns in the input than expected.\n' ...
                    '                    See >> help readlocs\n']);
       end;
       
       % removing lines BESA
       % -------------------
       if isempty(array{1,2})
           disp('BESA header detected, skipping three lines...');
           array = load_file_or_array( filename, g.skiplines-1);
           if isempty(array{1,2})
               array = load_file_or_array( filename, g.skiplines-1);
           end;
       end;

       % xyz format, is the first col absent
       % -----------------------------------
       if strcmp(g.filetype, 'xyz')
           if size(array, 2) == 4
               array(:, 2:5) = array(:, 1:4);
           end;
       end;
       
       % removing comments and empty lines
       % ---------------------------------
       indexbeg = 1;
       while isempty(array{indexbeg,1}) | ...
               (isstr(array{indexbeg,1}) & array{indexbeg,1}(1) == '%' )
           indexbeg = indexbeg+1;
       end;
       array = array(indexbeg:end,:);
       
       % converting file
       % ---------------
       for indexcol = 1:min(size(array,2), length(g.format))
           [str mult] = checkformat(g.format{indexcol});
           for indexrow = 1:size( array, 1)
               if mult ~= 1
                   eval ( [ 'eloc(indexrow).'  str '= -array{indexrow, indexcol};' ]);
               else
                   eval ( [ 'eloc(indexrow).'  str '= array{indexrow, indexcol};' ]);
               end;
           end;
       end;
   end;
   
   % handling BESA coordinates
   % -------------------------
   if isfield(eloc, 'sph_theta_besa')
       if isfield(eloc, 'type')
           if isnumeric(eloc(1).type)
               disp('BESA format detected ( Theta | Phi )');
               for index = 1:length(eloc)
                   eloc(index).sph_phi_besa   = eloc(index).labels;
                   eloc(index).sph_theta_besa = eloc(index).type;
                   eloc(index).labels         = '';
                   eloc(index).type           = '';
               end;
               eloc = rmfield(eloc, 'labels');
           end;
       end;
       if isfield(eloc, 'labels')       
           if isnumeric(eloc(1).labels)
               disp('BESA format detected ( Elec | Theta | Phi )');
               for index = 1:length(eloc)
                   eloc(index).sph_phi_besa   = eloc(index).sph_theta_besa;
                   eloc(index).sph_theta_besa = eloc(index).labels;
                   eloc(index).labels         = eloc(index).type;
                   eloc(index).type           = '';
                   eloc(index).radius         = 1;
               end;           
           end;
       end;
       
       try
           eloc = convertlocs(eloc, 'sphbesa2all');
           eloc = convertlocs(eloc, 'topo2all'); % problem with some EGI files (not BESA files)
       catch, disp('Warning: coordinate conversion failed'); end;
       fprintf('Readlocs: BESA spherical coords. converted, now deleting BESA fields\n');   
       fprintf('          to avoid confusion (these fields can be exported, though)\n');   
       eloc = rmfield(eloc, 'sph_phi_besa');
       eloc = rmfield(eloc, 'sph_theta_besa');

       % converting XYZ coordinates to polar
       % -----------------------------------
   elseif isfield(eloc, 'sph_theta')
       try
           eloc = convertlocs(eloc, 'sph2all');  
       catch, disp('Warning: coordinate conversion failed'); end;
   elseif isfield(eloc, 'X')
       try
           eloc = convertlocs(eloc, 'cart2all');  
       catch, disp('Warning: coordinate conversion failed'); end;
   else 
       try
           eloc = convertlocs(eloc, 'topo2all');  
       catch, disp('Warning: coordinate conversion failed'); end;
   end;
   
   % inserting labels if no labels
   % -----------------------------
   if ~isfield(eloc, 'labels')
       fprintf('readlocs(): Inserting electrode labels automatically.\n');
       for index = 1:length(eloc)
           eloc(index).labels = [ 'E' int2str(index) ];
       end;
   else 
       % remove trailing '.'
       for index = 1:length(eloc)
           if isstr(eloc(index).labels)
               tmpdots = find( eloc(index).labels == '.' );
               eloc(index).labels(tmpdots) = [];
           end;
       end;
   end;
   
   % resorting electrodes if number not-sorted
   % -----------------------------------------
   if isfield(eloc, 'channum')
       if ~isnumeric(eloc(1).channum)
           error('Channel numbers must be numeric');
       end;
       allchannum = [ eloc.channum ];
       if any( sort(allchannum) ~= allchannum )
           fprintf('readlocs(): Re-sorting channel numbers based on ''channum'' column indices\n');
           [tmp newindices] = sort(allchannum);
           eloc = eloc(newindices);
       end;
       eloc = rmfield(eloc, 'channum');      
   end;
else
    if isstruct(filename)
        % detect Fieldtrip structure and convert it
        % -----------------------------------------
        if isfield(filename, 'pnt')
            neweloc = [];
            for index = 1:length(filename.label)
                neweloc(index).labels = filename.label{index};
                neweloc(index).X      = filename.pnt(index,1);
                neweloc(index).Y      = filename.pnt(index,2);
                neweloc(index).Z      = filename.pnt(index,3);
            end;
            eloc = neweloc;
            eloc = convertlocs(eloc, 'cart2all');
        else
            eloc = filename;
        end;
    else
        disp('readlocs(): input variable must be a string or a structure');
    end;        
end;
if ~isempty(g.elecind)
	eloc = eloc(g.elecind);
end;
if nargout > 2
    if isfield(eloc, 'theta')
         tmptheta = { eloc.theta }; % check which channels have (polar) coordinates set
    else tmptheta = cell(1,length(eloc));
    end;
    if isfield(eloc, 'theta')
         tmpx = { eloc.X }; % check which channels have (polar) coordinates set
    else tmpx = cell(1,length(eloc));
    end;
    
    indices           = find(~cellfun('isempty', tmptheta));
    indices           = intersect_bc(find(~cellfun('isempty', tmpx)), indices);
    indices           = sort(indices);
    
    indbad            = setdiff_bc(1:length(eloc), indices);
    tmptheta(indbad)  = { NaN };
    theta             = [ tmptheta{:} ];
end;
if nargout > 3
    if isfield(eloc, 'theta')
         tmprad = { eloc.radius }; % check which channels have (polar) coordinates set
    else tmprad = cell(1,length(eloc));
    end;
    tmprad(indbad)    = { NaN };
    radius            = [ tmprad{:} ];
end;

%tmpnum = find(~cellfun('isclass', { eloc.labels }, 'char'));
%disp('Converting channel labels to string');
for index = 1:length(eloc)
    if ~isstr(eloc(index).labels)
        eloc(index).labels = int2str(eloc(index).labels);
    end;
end;
labels = { eloc.labels };
if isfield(eloc, 'ignore')
    eloc = rmfield(eloc, 'ignore');
end;

% process fiducials if any
% ------------------------
fidnames = { 'nz' 'lpa' 'rpa' 'nasion' 'left' 'right' 'nazion' 'fidnz' 'fidt9' 'fidt10' 'cms' 'drl' };
for index = 1:length(fidnames)
    ind = strmatch(fidnames{index}, lower(labels), 'exact');
    if ~isempty(ind), eloc(ind).type = 'FID'; end;
end;

return;

% interpret the variable name
% ---------------------------
function array = load_file_or_array( varname, skiplines );
	 if isempty(skiplines),
       skiplines = 0;
    end;
    if exist( varname ) == 2
        array = loadtxt(varname,'verbose','off','skipline',skiplines,'blankcell','off');
    else % variable in the global workspace
         % --------------------------
         try, array = evalin('base', varname);
	     catch, error('readlocs(): cannot find the named file or variable, check syntax');
		 end;
    end;     
return;

% check field format
% ------------------
function [str, mult] = checkformat(str)
	mult = 1;
	if strcmpi(str, 'labels'),         str = lower(str); return; end;
	if strcmpi(str, 'channum'),        str = lower(str); return; end;
	if strcmpi(str, 'theta'),          str = lower(str); return; end;
	if strcmpi(str, 'radius'),         str = lower(str); return; end;
	if strcmpi(str, 'ignore'),         str = lower(str); return; end;
	if strcmpi(str, 'sph_theta'),      str = lower(str); return; end;
	if strcmpi(str, 'sph_phi'),        str = lower(str); return; end;
	if strcmpi(str, 'sph_radius'),     str = lower(str); return; end;
	if strcmpi(str, 'sph_theta_besa'), str = lower(str); return; end;
	if strcmpi(str, 'sph_phi_besa'),   str = lower(str); return; end;
	if strcmpi(str, 'gain'),           str = lower(str); return; end;
	if strcmpi(str, 'calib'),          str = lower(str); return; end;
	if strcmpi(str, 'type') ,          str = lower(str); return; end;
	if strcmpi(str, 'X'),              str = upper(str); return; end;
	if strcmpi(str, 'Y'),              str = upper(str); return; end;
	if strcmpi(str, 'Z'),              str = upper(str); return; end;
	if strcmpi(str, '-X'),             str = upper(str(2:end)); mult = -1; return; end;
	if strcmpi(str, '-Y'),             str = upper(str(2:end)); mult = -1; return; end;
	if strcmpi(str, '-Z'),             str = upper(str(2:end)); mult = -1; return; end;
	if strcmpi(str, 'custom1'), return; end;
	if strcmpi(str, 'custom2'), return; end;
	if strcmpi(str, 'custom3'), return; end;
	if strcmpi(str, 'custom4'), return; end;
    error(['readlocs(): undefined field ''' str '''']);
   
    % finputcheck() - check Matlab function {'key','value'} input argument pairs
%
% Usage: >> result = finputcheck( varargin, fieldlist );
%        >> [result varargin] = finputcheck( varargin, fieldlist, ... 
%                                              callingfunc, mode, verbose );
% Input:
%   varargin  - Cell array 'varargin' argument from a function call using 'key', 
%               'value' argument pairs. See Matlab function 'varargin'.
%               May also be a structure such as struct(varargin{:})
%   fieldlist - A 4-column cell array, one row per 'key'. The first
%               column contains the key string, the second its type(s), 
%               the third the accepted value range, and the fourth the 
%               default value.  Allowed types are 'boolean', 'integer', 
%               'real', 'string', 'cell' or 'struct'.  For example,
%                       {'key1' 'string' { 'string1' 'string2' } 'defaultval_key1'}
%                       {'key2' {'real' 'integer'} { minint maxint } 'defaultval_key2'} 
%  callingfunc - Calling function name for error messages. {default: none}.
%  mode        - ['ignore'|'error'] ignore keywords that are either not specified 
%                in the fieldlist cell array or generate an error. 
%                {default: 'error'}.
%  verbose     - ['verbose', 'quiet'] print information. Default: 'verbose'.
%
% Outputs:
%   result     - If no error, structure with 'key' as fields and 'value' as 
%                content. If error this output contain the string error.
%   varargin   - residual varagin containing unrecognized input arguments.
%                Requires mode 'ignore' above.
%
% Note: In case of error, a string is returned containing the error message
%       instead of a structure.
%
% Example (insert the following at the beginning of your function):
%	result = finputcheck(varargin, ...
%               { 'title'         'string'   []       ''; ...
%                 'percent'       'real'     [0 1]    1 ; ...
%                 'elecamp'       'integer'  [1:10]   [] });
%   if isstr(result)
%       error(result);
%   end
%
% Note: 
%   The 'title' argument should be a string. {no default value}
%   The 'percent' argument should be a real number between 0 and 1. {default: 1}
%   The 'elecamp' argument should be an integer between 1 and 10 (inclusive).
%
%   Now 'g.title' will contain the title arg (if any, else the default ''), etc.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 10 July 2002

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 10 July 2002, arno@salk.edu
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

function [g, varargnew] = finputcheck( vararg, fieldlist, callfunc, mode, verbose )

	if nargin < 2
		help finputcheck;
		return;
	end;
	if nargin < 3
		callfunc = '';
	else 
		callfunc = [callfunc ' ' ];
	end;
    if nargin < 4
        mode = 'do not ignore';
    end;
    if nargin < 5
        verbose = 'verbose';
    end;
	NAME = 1;
	TYPE = 2;
	VALS = 3;
	DEF  = 4;
	SIZE = 5;
	
	varargnew = {};
	% create structure
	% ----------------
	if ~isempty(vararg)
        if isstruct(vararg)
            g = vararg;
        else
            for index=1:length(vararg)
                if iscell(vararg{index})
                    vararg{index} = {vararg{index}};
                end;
            end;
            try
                g = struct(vararg{:});
            catch
                vararg = removedup(vararg, verbose);
                try
                    g = struct(vararg{:});
                catch
                    g = [ callfunc 'error: bad ''key'', ''val'' sequence' ]; return;
                end;
            end;
        end;
	else 
		g = [];
	end;
	
	for index = 1:size(fieldlist,NAME)
		% check if present
		% ----------------
		if ~isfield(g, fieldlist{index, NAME})
			g = setfield( g, fieldlist{index, NAME}, fieldlist{index, DEF});
		end;
		tmpval = getfield( g, {1}, fieldlist{index, NAME});
		
		% check type
		% ----------
        if ~iscell( fieldlist{index, TYPE} )
            res = fieldtest( fieldlist{index, NAME},  fieldlist{index, TYPE}, ...
                           fieldlist{index, VALS}, tmpval, callfunc );
            if isstr(res), g = res; return; end;
        else 
            testres = 0;
            tmplist = fieldlist;
            for it = 1:length( fieldlist{index, TYPE} )
                if ~iscell(fieldlist{index, VALS})
                     res{it} = fieldtest(  fieldlist{index, NAME},  fieldlist{index, TYPE}{it}, ...
                                           fieldlist{index, VALS}, tmpval, callfunc );
                else res{it} = fieldtest(  fieldlist{index, NAME},  fieldlist{index, TYPE}{it}, ...
                                           fieldlist{index, VALS}{it}, tmpval, callfunc );
                end;
                if ~isstr(res{it}), testres = 1; end;
            end;
            if testres == 0,
                g = res{1};
                for tmpi = 2:length(res)
                    g = [ g 10 'or ' res{tmpi} ];
                end;
                return; 
            end;
        end;
	end;
    
    % check if fields are defined
	% ---------------------------
	allfields = fieldnames(g);
	for index=1:length(allfields)
		if isempty(strmatch(allfields{index}, fieldlist(:, 1)', 'exact'))
			if ~strcmpi(mode, 'ignore')
				g = [ callfunc 'error: undefined argument ''' allfields{index} '''']; return;
			end;
			varargnew{end+1} = allfields{index};
			varargnew{end+1} = getfield(g, {1}, allfields{index});
		end;
	end;


function g = fieldtest( fieldname, fieldtype, fieldval, tmpval, callfunc );
	NAME = 1;
	TYPE = 2;
	VALS = 3;
	DEF  = 4;
	SIZE = 5;
    g = [];
    
    switch fieldtype
     case { 'integer' 'real' 'boolean' 'float' }, 
      if ~isnumeric(tmpval) && ~islogical(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be numeric' ]; return;
      end;
      if strcmpi(fieldtype, 'boolean')
          if tmpval ~=0 && tmpval ~= 1
              g = [ callfunc 'error: argument ''' fieldname ''' must be 0 or 1' ]; return;
          end;  
      else 
          if strcmpi(fieldtype, 'integer')
              if ~isempty(fieldval)
                  if (any(isnan(tmpval(:))) && ~any(isnan(fieldval))) ...
                          && (~ismember(tmpval, fieldval))
                      g = [ callfunc 'error: wrong value for argument ''' fieldname '''' ]; return;
                  end;
              end;
          else % real or float
              if ~isempty(fieldval) && ~isempty(tmpval)
                  if any(tmpval < fieldval(1)) || any(tmpval > fieldval(2))
                      g = [ callfunc 'error: value out of range for argument ''' fieldname '''' ]; return;
                  end;
              end;
          end;
      end;  
      
      
     case 'string'
      if ~isstr(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be a string' ]; return;
      end;
      if ~isempty(fieldval)
          if isempty(strmatch(lower(tmpval), lower(fieldval), 'exact'))
              g = [ callfunc 'error: wrong value for argument ''' fieldname '''' ]; return;
          end;
      end;

      
     case 'cell'
      if ~iscell(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be a cell array' ]; return;
      end;
      
      
     case 'struct'
      if ~isstruct(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be a structure' ]; return;
      end;
      
     case 'function_handle'
      if ~isa(tmpval, 'function_handle')
          g = [ callfunc 'error: argument ''' fieldname ''' must be a function handle' ]; return;
      end;
      
     case '';
     otherwise, error([ 'finputcheck error: unrecognized type ''' fieldname '''' ]);
    end;

% remove duplicates in the list of parameters
% -------------------------------------------
function cella = removedup(cella, verbose)
% make sure if all the values passed to unique() are strings, if not, exist
%try
    [tmp indices] = unique_bc(cella(1:2:end));
    if length(tmp) ~= length(cella)/2
        myfprintf(verbose,'Note: duplicate ''key'', ''val'' parameter(s), keeping the last one(s)\n');
    end;
    cella = cella(sort(union(indices*2-1, indices*2)));
%catch
    % some elements of cella were not string
%    error('some ''key'' values are not string.');
%end;    

function myfprintf(verbose, varargin)

if strcmpi(verbose, 'verbose')
    fprintf(varargin{:});
end;

% loadtxt() - load ascii text file into numeric or cell arrays
%
% Usage:
%   >> array = loadtxt( filename, 'key', 'val' ...);
%
% Inputs:
%    filename - name of the input file
%
% Optional inputs
%   'skipline' - number of lines to skip {default:0}. If this number is
%                negative the program will only skip non-empty lines 
%                (can be usefull for files transmitted from one platform
%                to an other, as CR may be inserted at every lines).
%   'convert'  - 'on' standard text conversion, see note 1
%                'off' no conversion, considers text only
%                'force' force conversion, NaN are returned 
%                for non-numeric inputs {default:'on'}
%   'delim'    - ascii character for delimiters. {default:[9 32]
%                i.e space and tab}. It is also possible to enter 
%                strings, Ex: [9 ' ' ','].
%   'blankcell' - ['on'|'off'] extract blank cells {default:'on'}
%   'verbose'  - ['on'|'off'] {default:'on'}
%   'convertmethod' - ['str2double'|'str2num'] default is 'str2double'
%   'nlines'   - [integer] number of lines to read {default: all file}
%
% Outputs:
%    array - cell array. If the option 'force' is given, the function
%            retrun a numeric array.
%
% Notes: 1) Since it uses cell arrays, the function can handle text input.
%        The function reads each token and then try to convert it to a 
%        number. If the conversion is unsucessfull, the string itself
%        is included in the array.
%        2) The function adds empty entries for rows that contains
%        fewer columns than others.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 29 March 2002

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 29 March 2002
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

function array = loadtxt( filename, varargin );

if nargin < 1
	help loadtxt;
	return;
end;	
if ~isempty(varargin)
   try, g = struct(varargin{:});
   catch, disp('Wrong syntax in function arguments'); return; end;
else
    g = [];
end;

g = finputcheck( varargin, { 'convert'   'string'   { 'on';'off';'force' }   'on';
                             'skipline'  'integer'  [0 Inf]          0;
                             'verbose'   'string'   { 'on';'off' }   'on';
                             'uniformdelim' 'string'   { 'on';'off' }   'off';                             
                             'blankcell' 'string'   { 'on';'off' }   'on';
                             'convertmethod' 'string'   { 'str2double';'str2num' }   'str2double';
                             'delim'     { 'integer';'string' } []               [9 32];
                             'nlines'    'integer'  []               Inf });
if isstr(g), error(g); end;
if strcmpi(g.blankcell, 'off'), g.uniformdelim = 'on'; end;
g.convert = lower(g.convert);
g.verbose = lower(g.verbose);
g.delim = char(g.delim);

% open the file
% -------------
if exist(filename) ~=2, error( ['file ' filename ' not found'] ); end;  
fid=fopen(filename,'r','ieee-le');
if fid<0, error( ['file ' filename ' found but error while opening file'] ); end;  

index = 0;
while index < abs(g.skipline)
    tmpline = fgetl(fid); 
    if g.skipline > 0 | ~isempty(tmpline)
        index = index + 1;
    end;    
end; % skip lines ---------

inputline = fgetl(fid);
linenb = 1;
if strcmp(g.verbose, 'on'), fprintf('Reading file (lines): '); end;
while isempty(inputline) | inputline~=-1
     colnb = 1;
     if ~isempty(inputline)
         tabFirstpos = 1;
         
         % convert all delimiter to the first one
         if strcmpi(g.uniformdelim, 'on')
             for index = 2:length(g.delim)
                 inputline(find(inputline == g.delim(index))) = g.delim(1);
             end;
         end;
         
         while ~isempty(deblank(inputline))
             if strcmpi(g.blankcell,'off'), inputline = strtrim(inputline); end;
             if tabFirstpos && length(inputline) > 1 && all(inputline(1) ~= g.delim), tabFirstpos = 0; end;
             [tmp inputline tabFirstpos] = mystrtok(inputline, g.delim, tabFirstpos);
             switch g.convert
                case 'off', array{linenb, colnb} = tmp;
                case 'on',  
                     if strcmpi(g.convertmethod, 'str2double')
                         tmp2 = str2double(tmp);
                         if isnan( tmp2 )  , array{linenb, colnb} = tmp;
                         else                array{linenb, colnb} = tmp2;
                         end;
                     else
                         tmp2 = str2num(tmp);
                         if isempty( tmp2 )  , array{linenb, colnb} = tmp;
                         else                  array{linenb, colnb} = tmp2;
                         end;
                     end;
                case 'force', array{linenb, colnb} = str2double(tmp);
             end;
             colnb = colnb+1;
         end;
	     linenb = linenb +1;
     end;
     inputline = fgetl(fid);
     if linenb > g.nlines
         inputline = -1;
     end;
     if ~mod(linenb,10) & strcmp(g.verbose, 'on'), fprintf('%d ', linenb); end;
end;        
if strcmp(g.verbose, 'on'),  fprintf('%d\n', linenb-1); end;
if strcmp(g.convert, 'force'), array = [ array{:} ]; end;
fclose(fid); 

% problem strtok do not consider tabulation
% -----------------------------------------
function [str, strout, tabFirstpos] = mystrtok(strin, delim, tabFirstpos)
    % remove extra spaces at the beginning
    while any(strin(1) == delim) && strin(1) ~= 9 && strin(1) ~= ','
         strin = strin(2:end);
    end;
    % for tab and coma, consider empty cells
    if length(strin) > 1 && any(strin(1) == delim)
        if tabFirstpos || any(strin(2) == delim)
            str = '';
            strout = strin(2:end);
            if strin(2) ~= 9 && strin(2) ~= ','
                tabFirstpos = 0;
                strout = strtrim(strout);
            end;
        else
            [str, strout] = strtok(strin, delim);
        end;
    else
        [str, strout] = strtok(strin, delim);
    end;

    % convertlocs() - Convert electrode locations between coordinate systems
%                 using the EEG.chanlocs structure.
%
% Usage: >> newchans = convertlocs( EEG, 'command');
%
% Input:
%   chanlocs  - An EEGLAB EEG dataset OR a EEG.chanlocs channel locations structure
%   'command' - ['cart2topo'|'sph2topo'|'sphbesa2topo'| 'sph2cart'|'topo2cart'|'sphbesa2cart'|
%               'cart2sph'|'sphbesa2sph'|'topo2sph'| 'cart2sphbesa'|'sph2sphbesa'|'topo2sphbesa'|
%               'cart2all'|'sph2all'|'sphbesa2all'|'topo2all']
%                These command modes convert between four coordinate frames: 3-D Cartesian 
%                (cart), Matlab spherical (sph), Besa spherical (sphbesa), and 2-D polar (topo)
%               'auto' -- Here, the function finds the most complex coordinate frame 
%                 and constrains all the others to this one. It searches first for Cartesian 
%                 coordinates, then for spherical and finally for polar. Default is 'auto'.
%
% Optional input
%   'verbose' - ['on'|'off'] default is 'off'.
%
% Outputs:
%   newchans - new EEGLAB channel locations structure
%
% Ex:  CHANSTRUCT = convertlocs( CHANSTRUCT, 'cart2topo');
%      % Convert Cartesian coordinates to 2-D polar (topographic). 
%
% Author: Arnaud Delorme, CNL / Salk Institute, 22 Dec 2002
%
% See also: readlocs()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 22 Dec 2002, arno@salk.edu
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

function chans = convertlocs(chans, command, varargin)

if nargin < 1
   help convertlocs;
   return;
end;

if nargin < 2
   command = 'auto';
end;
if nargin == 4 && strcmpi(varargin{2}, 'on')
    verbose = 1;
else
    verbose = 0; % off
end;

% test if value exists for default
% --------------------------------
if strcmp(command, 'auto')
    if isfield(chans, 'X') && ~isempty(chans(1).X)
        command = 'cart2all';
        if verbose
            disp('Make all coordinate frames uniform using Cartesian coords');
        end;
    else
        if isfield(chans, 'sph_theta') && ~isempty(chans(1).sph_theta)
            command = 'sph2all';
            if verbose
                disp('Make all coordinate frames uniform using spherical coords');
            end;
        else
            if isfield(chans, 'sph_theta_besa') && ~isempty(chans(1).sph_theta_besa)
                command = 'sphbesa2all';
                if verbose
                    disp('Make all coordinate frames uniform using BESA spherical coords');
                end;
            else
                command = 'topo2all';
                if verbose
                    disp('Make all coordinate frames uniform using polar coords');
                end;
            end;
        end;
    end;
end;

% convert
% -------         
switch command
 case 'topo2sph',
   theta  = {chans.theta};
   radius = {chans.radius};
   indices = find(~cellfun('isempty', theta));
   [sph_phi sph_theta] = topo2sph( [ [ theta{indices} ]' [ radius{indices}]' ] );
   if verbose
       disp('Warning: electrodes forced to lie on a sphere for polar to 3-D conversion');
   end;
   for index = 1:length(indices)
      chans(indices(index)).sph_theta  = sph_theta(index);
      chans(indices(index)).sph_phi    = sph_phi  (index);
   end;
   if isfield(chans, 'sph_radius'),
       meanrad = mean([ chans(indices).sph_radius ]);
       if isempty(meanrad), meanrad = 1; end;
   else
       meanrad = 1;
   end;
   sph_radius(1:length(indices)) = {meanrad};
case 'topo2sphbesa',
   chans = convertlocs(chans, 'topo2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
case 'topo2cart'
   chans = convertlocs(chans, 'topo2sph', varargin{:}); % search for spherical coords
   if verbose
       disp('Warning: spherical coordinates automatically updated');
   end;
   chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords
case 'topo2all',
   chans = convertlocs(chans, 'topo2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords
case 'sph2cart',
   sph_theta  = {chans.sph_theta};
   sph_phi    = {chans.sph_phi};
   indices = find(~cellfun('isempty', sph_theta));
   if ~isfield(chans, 'sph_radius'), sph_radius(1:length(indices)) = {1};
   else                              sph_radius = {chans.sph_radius};
   end;
   inde = find(cellfun('isempty', sph_radius));
   if ~isempty(inde)
       meanrad = mean( [ sph_radius{:} ]);
       sph_radius(inde) = { meanrad };
   end;
   [x y z] = sph2cart([ sph_theta{indices} ]'/180*pi, [ sph_phi{indices} ]'/180*pi, [ sph_radius{indices} ]');
   for index = 1:length(indices)
      chans(indices(index)).X = x(index);
      chans(indices(index)).Y = y(index);
      chans(indices(index)).Z = z(index);
   end;
case 'sph2topo',
 if verbose
     % disp('Warning: all radii constrained to one for spherical to topo transformation');
 end;
 sph_theta  = {chans.sph_theta};
 sph_phi    = {chans.sph_phi};
 indices = find(~cellfun('isempty', sph_theta));
 [chan_num,angle,radius] = sph2topo([ ones(length(indices),1)  [ sph_phi{indices} ]' [ sph_theta{indices} ]' ], 1, 2); % using method 2
 for index = 1:length(indices)
     chans(indices(index)).theta  = angle(index);
     chans(indices(index)).radius = radius(index);
     if ~isfield(chans, 'sph_radius') || isempty(chans(indices(index)).sph_radius)
         chans(indices(index)).sph_radius = 1;
     end;
 end;
case 'sph2sphbesa',
   % using polar coordinates
   sph_theta  = {chans.sph_theta};
   sph_phi    = {chans.sph_phi};
   indices = find(~cellfun('isempty', sph_theta));
   [chan_num,angle,radius] = sph2topo([ones(length(indices),1)  [ sph_phi{indices} ]' [ sph_theta{indices} ]' ], 1, 2);
   [sph_theta_besa sph_phi_besa] = topo2sph([angle radius], 1, 1);
   for index = 1:length(indices)
      chans(indices(index)).sph_theta_besa  = sph_theta_besa(index);
      chans(indices(index)).sph_phi_besa    = sph_phi_besa(index);
   end;   
case 'sph2all',
   chans = convertlocs(chans, 'sph2topo', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords
case 'sphbesa2sph',
   % using polar coordinates
   sph_theta_besa  = {chans.sph_theta_besa};
   sph_phi_besa    = {chans.sph_phi_besa};
   indices = find(~cellfun('isempty', sph_theta_besa));
   [chan_num,angle,radius] = sph2topo([ones(length(indices),1)  [ sph_theta_besa{indices} ]' [ sph_phi_besa{indices} ]' ], 1, 1);
   %for index = 1:length(chans)
   %   chans(indices(index)).theta  = angle(index);
   %   chans(indices(index)).radius = radius(index);
   %   chans(indices(index)).labels = int2str(index);
   %end;   
   %figure; topoplot([],chans, 'style', 'blank', 'electrodes', 'labelpoint');
   
   [sph_phi sph_theta] = topo2sph([angle radius], 2);
   for index = 1:length(indices)
      chans(indices(index)).sph_theta  = sph_theta(index);
      chans(indices(index)).sph_phi    = sph_phi  (index);      
   end;
case 'sphbesa2topo',
   chans = convertlocs(chans, 'sphbesa2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2topo', varargin{:}); % search for spherical coords
case 'sphbesa2cart',
   chans = convertlocs(chans, 'sphbesa2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords   
case 'sphbesa2all',
   chans = convertlocs(chans, 'sphbesa2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2all', varargin{:}); % search for spherical coords
case 'cart2topo',
   chans = convertlocs(chans, 'cart2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2topo', varargin{:}); % search for spherical coords
case 'cart2sphbesa',
   chans = convertlocs(chans, 'cart2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
case 'cart2sph',
    if verbose
        disp('WARNING: If XYZ center has not been optimized, optimize it using Edit > Channel Locations');
	end;
    X  = {chans.X};
    Y  = {chans.Y};
    Z  = {chans.Z};
    indices = find(~cellfun('isempty', X));
    [th phi radius] = cart2sph( [ X{indices} ], [ Y{indices} ], [ Z{indices} ]);
	for index = 1:length(indices)
		 chans(indices(index)).sph_theta     = th(index)/pi*180;
		 chans(indices(index)).sph_phi       = phi(index)/pi*180;
		 chans(indices(index)).sph_radius    = radius(index);
	end;
case 'cart2all',
   chans = convertlocs(chans, 'cart2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2all', varargin{:}); % search for spherical coords
end;
