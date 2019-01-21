function [p,data]=nt_dataview(data,p)
%[p,data]=nt_dataview(data,p) - view data sets
%
% 
%  DATA: matrix, struct, file or directory to view
%  P: parameter structure
%  
% 
VERBOSE=1;

if nargin < 1;
    % look for file in entire machine
    [p,data]=nt_dataview('/');
    return
end

if nargin < 2; p=[]; end
if isempty(data); return; end

% % update string to recreate this view
% if ~isfield(p,'recreate');
%     % first call, call string depends on number of arguments
%     if nargout==2; 
%         s1='[p,data]='; 
%     elseif nargout==1; 
%         s1='[p]='; 
%     else
%         s1='';
%     end
%     if isa(data,'char')
%         s2=data;
%     else
%         s2=inputname(1);
%     end
%     if nargout==2;
%         s3=',p);';
%     else
%         s3=');';
%     end
%     p.recreate=[s1,'nt_dataview(',s2,s3];
% else
%     % append new call
%     p.recreate=([p.recreate,'; [p,data]=nt_data_view(data,p);']);
% end


% name to display on dialog or window
if ~isfield(p,'data_name')
    if isa(data,'char');
        p.data_name=data;
    elseif isa(data,'struct');
        p.data_name='structure';
    elseif isnumeric(data) %isa(data,'double');
        p.data_name='matrix';
    else
        error('argument should be string, struct, or numeric');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%   MATRIX  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isnumeric(data) 
    
    if ~isa(data, 'double')
        warning('converting data to double');
        data=double(data);
    end
    
    if ~isempty(p) && isfield(p, 'matrix_modify'); 
        eval([p.matrix_modify,';']); % modify data
    end
    
    nt_whoss;
    if VERBOSE; disp('matrix'); end
    nDims=ndims(data);
    if nDims==2 && (size(data,1)==1 || size(data,2)==1); nDims=1; end
    if nDims>4; nDims=4; end
    
    % positions
    posFig=[0 100, 1000, 400];
    posButtonReturn=[50, 20, 100, 25];
    posButtonAssign=[200, 20, 300, 25];
    posEdit=[100 70 800 30];
    
    % put up window
    dialogH=figure('Visible','on', 'WindowStyle','normal', 'Name',p.data_name, 'NumberTitle','off',...
        'Pointer','arrow', 'Position',posFig,'color', [1 1 1], 'KeyPressFcn',@doFigureKeyPress,...
        'IntegerHandle','off', 'CloseRequestFcn' ,@doDelete);
    editH=[]; buttonReturnH=[];
    set(gcf,'userdata',p);
    done=0;
    while ~done
    
        % plot summary statistics for all dimensions
        plot_params.bottom=0.35; plot_params.height=0.5;
        
        %%%%%%%%%%%%%%%%%%  HERE's where we plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nt_statmatrix(data,plot_params);
        
        if isfield(p,'return') && p.return; break; end

        %set(dialogH,'HandleVisibility','callback');
        whichButton=[];
        returnPressed=[];
        escapePressed=[];
        otherKeyPressed=[];

        % create return button
        buttonReturnH=uicontrol('style', 'pushbutton', 'string', 'Return','position', posButtonReturn,...
            'KeyPressFcn',@doControlKeyPress, 'Callback',@doButtonPress, 'HorizontalAlignment','center', 'parent', dialogH,...
            'fontSize', 14);

        % create return button
        buttonAssignH=uicontrol('style', 'pushbutton', 'string', 'Assign to p,data in workspace','position', posButtonAssign,...
            'KeyPressFcn',@doControlKeyPress, 'Callback',@doButtonPress, 'HorizontalAlignment','center', 'parent', dialogH,...
            'fontSize', 14);

        editString=['data=data( 1 : ',num2str(size(data,1))];
        for iDim=2:nDims
            editString=[editString,', 1 : ',num2str(size(data,iDim))]; % full index for that dimension
        end
        editString=[editString,' );'];
        editH=uicontrol('Style','edit','String',editString,'position', posEdit, 'parent', dialogH,...
            'callback', @editCallback, 'foregroundcolor',[1 1 1]*0 );
        % wait for user input
        if ~ ishghandle(dialogH); return; end 
        uiwait(dialogH); 
        if ~ ishghandle(dialogH); return; end

        % act on user input
        if returnPressed % keyboard
            done=1;
        elseif escapePressed % keyboard
            done=1;
        elseif otherKeyPressed % keyboard
            ;
        else
            h=get(dialogH,'CurrentObject');

            if find(h==editH) % one of the edit boxes
                s=get(editH(1),'string');
                try
                   th=annotation('textbox', [.5 .04 .1 .1], 'string', 'evaluating...', 'fontsize', 14, 'edgecolor', [1 1 1]);
                   drawnow
                   if ~isempty(p) && isfield(p, 'matrix_modify'); 
                       p.matrix_modify=[p.matrix_modify, s];
                   else
                       p.matrix_modify=s;
                   end
                   eval(s);
                   set(th,'string','');
                   clf
                catch
                   beep;
                   warning(['incorrect indexing string: >',s,'<']);
                end
            elseif h==buttonReturnH 
                done=1;
            elseif h==buttonAssignH 
                assignin('base','p',p);
                assignin('base','data',data);
                done=1;
            else
                disp(num2str(h));
                error('unexpected handle')
            end
        end
        %clf
    end
    
    set(gcf,'userdata',p);
    if ~isempty(editH); delete(editH); end
    if ~isempty(buttonReturnH); delete(buttonReturnH); end
    if ishghandle(dialogH); delete(dialogH); end
    %set(buttonH, 'string', 'Recreate','Callback',@doButtonPress2);

    
    % return data - or not
    if nargout==0; clear data p; end
    
    if VERBOSE; disp('returning from nt_dataview...'); end
    %return;
 


%%%%%%%%%%%%%%%%%%%%%%%%   STRUCT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif isa(data, 'struct')
    % struct in workspace
    if VERBOSE; disp('struct'); end
    if isempty(struct); error('''struct'' is empty!'); end
    field_names=fieldnames(data);
    if isempty(field_names); error('structure is empty!'); end
    field_sizes=zeros(numel(field_names),1);
    for k=1:numel(field_names);
        field=getfield(data,field_names{k});
        field_sizes(k)=round(prod(size(field))/1024);
    end
    clear field;
    a=repmat(' (', numel(field_names),1);
    b=cellstr(num2str(field_sizes, '%9d'));
    b=char(b);
    c=[repmat(' Mb)', numel(field_names),1)];
    i=listdlg('liststring',cellstr([char(field_names),a,b,c]),...
        'name', 'Select field in struct:', ...
        'listsize', [600 300], ...
        'OKstring','Select',...
        'PromptString',p.data_name);
    
    % call this function on the selected field
    if i
        data=getfield(data,field_names{i});
        [p,data]=nt_dataview(data,p); 
        if nargout==0; data=[]; end
        return
    end
    

    
elseif isa(data, 'char') && ...
        ( exist(data,'file')==2  ||  ...
        (numel(data)>3 && all(data(end-2:end)=='.ds'))) % treat as file

%%%%%%%%%%%%%%%%%%%%%%%%   FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    fname=data;
    if VERBOSE; disp('file'); end
    
    if 0
    % intercept mat files
    if numel(fname)>4 & fname(end-3:end)=='.mat'
        if VERBOSE; disp('mat file'); end
        S=whos('-file',fname);
        var_names=char(S.name);
        var_sizes=round([S.bytes]/1024)';
        a=repmat(' (', size(var_names,1),1);
        b=cellstr(num2str(var_sizes, '%9d'));
        b=char(b);
        c=[repmat(' Mb)', size(var_names,1),1)];
        i=listdlg('liststring',cellstr([var_names,a,b,c]),...
            'name', 'Select variable in file:', ...
            'listsize', [600 300], ...
            'OKstring','Select',...
            'PromptString',p.data_name);
        if i
            X=load(fname,var_names(i,:));
            [p,data]=nt_dataview(X,p);
        end
        if nargout==0; data=[]; end
        return
    end
    end
        
    % hand over to data file reader
    [p,data]=nt_read_data(fname);
    [p,data]=nt_dataview(data,p);
    return
    
%%%%%%%%%%%%%%%%%%%%%%%%   DIRECTORY  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif isa(data, 'char') && exist(data, 'dir')==7
    % directory
    if VERBOSE; disp('directory'); end
    
    d=dir(data);
    if numel(d)==0
        error(['directory >',data,'< not found']);
    end
    fnames=char(d.name);
    idx=find(fnames(:,1)~='.');  % remove '.' and '..' and invisible files
    d=d(idx);
    fnames=fnames(idx,:);
    
    % separate directories and files
    didx=find([d.isdir]);
    fidx=find(~[d.isdir]);
    fnames=fnames([didx, fidx],:);
    
   % count files within the directories
    nfiles=zeros(numel(didx),1);
    for k=1:numel(didx)
        dd=dir([data,'/',d(didx(k)).name]);
        fns=char(dd.name);
        idx=find(fns(:,1)~='.');  % remove '.' and '..' and invisible files
        nfiles(k)=numel(idx);
    end
    
    % size of the files
    mbytes=[d(fidx).bytes]'/1024;
   
    % string arrays to put in dialog list
    a=repmat(' (', numel(d),1);
    if numel(didx)>0
        b=cellstr(num2str(nfiles, '%9d'));
    else
        b=[]; % stupid matlab!
    end
    if numel(fidx)>0
        b=[b;cellstr(num2str(mbytes,'%0.1f'))];
    end
    b=char(b);
    c=[repmat(' files)', numel(didx),1); repmat(' Mb)   ', numel(fidx),1)];
     
    % which directory or file is user interested in?
    
    i=listdlg('liststring',cellstr([fnames,a,b,c]),...
        'name', 'Select file:', ...
        'listsize', [300 300], ...
        'OKstring','Select',...
        'PromptString',p.data_name);
    
    % call this function on that file
    if i
        data=strcat(data,'/',fnames(i,:));
        [p,data]=nt_dataview(data,p);   
    end
    return
else
    disp([p.data_name,' not found']); 
    if nargout==0; data=[]; end
    return
end


%h=data;
if nargout==0; 
    disp('hereiam');
    clear data p 
end




%%%%%%%%%%%%%%%%%%%%  LOCAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function doFigureKeyPress(obj, evd)  %#ok
        switch(evd.Key)
            case {'return','space'}
                returnPressed = true;
            case 'escape'
                escapePressed=true;
            otherwise
                otherKeyPressed=true;
        end
        uiresume(gcbf)
    end

    function doDelete(varargin)
        delete(dialogH);
    end

    function editCallback(obj,evd);
        %editString = get(obj,'String');
        uiresume(gcbf);
    end

    function doButtonPress(obj,evd);
        whichButton=obj;
        uiresume(gcbf);
    end
    function doButtonPress2(obj,evd);
        whichButton=obj;
        p=get(gcf,'userdata');
        evalin('caller',p.recreate);
    end



end % this file's main function
