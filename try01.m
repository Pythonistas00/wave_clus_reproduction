function varargout = try01(varargin)
% TRY01 MATLAB code for try01.fig
%      TRY01, by itself, creates a new TRY01 or raises the existing
%      singleton*.
%
%      H = TRY01 returns the handle to a new TRY01 or the handle to
%      the existing singleton*.
%
%      TRY01('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRY01.M with the given input arguments.
%
%      TRY01('Property','Value',...) creates a new TRY01 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before try01_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to try01_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help try01

% Last Modified by GUIDE v2.5 24-Mar-2022 18:00:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @try01_OpeningFcn, ...
                   'gui_OutputFcn',  @try01_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before try01 is made visible.
function try01_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to try01 (see VARARGIN)

% Choose default command line output for try01
handles.output = hObject;

set(handles.plot_average_button,'value',1);
set(handles.plot_all_button,'value',0);

% Update handles structure
guidata(hObject, handles);

if nargin>3 && ischar(varargin{1})
      load_data_button_Callback('w_arg',varargin{1},handles)
end
% UIWAIT makes try01 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = try01_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

clus_colors = [[0.0 0.0 1.0];[1.0 0.0 0.0];[0.0 0.5 0.0];[0.620690 0.0 0.0];
    [0.413793 0.0 0.758621];[0.965517 0.517241 0.034483];[0.448276 0.379310 0.241379];
    [1.0 0.103448 0.724138];[0.545 0.545 0.545];[0.586207 0.827586 0.310345];
    [0.965517 0.620690 0.862069];[0.620690 0.758621 1.]];


set(0,'DefaultAxesColorOrder',clus_colors)


% --- Executes on button press in load_data_button.
function load_data_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_data_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% -------------------------select file----------------------------
if ischar(eventdata)
    if isempty(fileparts(eventdata))
        filename = eventdata;
        pathname = [pwd filesep];
    else
        [pathname,filename,ext]=fileparts(eventdata);
        pathname = [pathname filesep];
        filename = [filename ext];
    end
    hObject = handles.load_data_button;
else
    [filename, pathname] = uigetfile('*.*','Select file'); % Use only the supported extensions can bring case-sensitive related problems.
end

% if any file was selected, cancel the loading
if ~ischar(pathname)
    return
end

set(handles.file_name,'string',['Loading: ' pathname filename]); drawnow

cla(handles.cont_data);
clear functions     % reset functions, force to reload set_parameters next
handles.par = set_parameters();

cd(pathname);
handles.par.filename = filename;

% --------------------------spike detection-----------------------------
data_handler = readInData(handles.par);
handles.par = data_handler.par;

handles.par.fname_in = 'tmp_data_wc';                     % temporary filename used as input for SPC
handles.par.fname = ['data_' data_handler.nick_name];
handles.par.nick_name = data_handler.nick_name;
handles.par.fnamesave = handles.par.fname;                % filename if "save clusters" button is pressed
handles.par.fnamespc = 'data_wc';
handles.par.spikes_file = data_handler.spikes_file;  
handles.par = data_handler.update_par(handles.par);
check_WC_params(handles.par)

if data_handler.with_results    % data have _times files
    [clu, tree, spikes, index, inspk, ipermut, classes, forced, temp] = data_handler.load_results();
    rejected = data_handler.load_rejected();
    handles.setclus = 1;
    if isempty(ipermut)
        handles.par.permut = 'n';
    end

else

    if data_handler.with_spikes     % data have some time of _spikes files
        [spikes, index] = data_handler.load_spikes();
        if ~data_handler.with_wc_spikes
            [spikes] = spike_alignment(spikes,handles.par);
        end
    else
        set(handles.file_name,'string','Detecting spikes ...'); drawnow
        index = [];
        spikes = [];
        for n = 1:data_handler.max_segments
            x = data_handler.get_segment();
            [new_spikes, temp_aux_th, new_index]  = amp_detect(x, handles.par);
            index = [index data_handler.index2ts(new_index)];    % new_index to ms
            spikes = [spikes; new_spikes];
        end
        handles.par.detection_date =  datestr(now);
    end

    if size(spikes,1) < 15
        	ME = MException('MyComponent:notEnoughSpikes', 'Less than 15 spikes detected');
            throw(ME)
    end

    set(handles.file_name,'string','Calculating spike features ...'); drawnow
    [inspk] = wave_features(spikes,handles.par);              % Extract spike features.
    handles.par.inputs = size(inspk,2);                       % number of inputs to the clustering

    if handles.par.permut == 'y'
        if handles.par.match == 'y'
            naux = min(handles.par.max_spk,size(inspk,1));
            ipermut = randperm(length(inspk));
            ipermut(naux+1:end) = [];
        else
            ipermut = randperm(length(inspk));
        end
        inspk_aux = inspk(ipermut,:);
    else
        if handles.par.match == 'y'
            naux = min(handles.par.max_spk,size(inspk,1));
            inspk_aux = inspk(1:naux,:);
        else
            inspk_aux = inspk;
        end
    end

    % Interaction with SPC
    set(handles.file_name,'string','Running SPC ...'); drawnow
    fname_in = handles.par.fname_in;
    save(fname_in,'inspk_aux','-ascii');               % Input file for SPC

    [clu,tree] = run_cluster(handles.par);
    forced = false(size(spikes,1) ,1);
    rejected = false(1, size(spikes,1));
    handles.setclus = 2;    % uses min cluster size but doesn't reset force
end

handles.par.file_name_to_show = [pathname filename];

% draw raw data
if (data_handler.with_raw || data_handler.with_psegment) && handles.par.cont_segment       % raw exists
    [xd_sub, sr_sub] = data_handler.get_signal_sample();
    Plot_continuous_data(xd_sub, sr_sub, handles); drawnow
    clear xd_sub
end
% above run well 

% -----------------------------cluster-------------------------------
% Fixing lost elements of clu . Skiped elements will be  class -1 because in
% all the uses of clu are like: clu(temp,3:end)+1
if handles.par.permut == 'y' && ~isempty(clu)
    if isempty(ipermut)          % load from old result without ipermut or par, but par.permut=='y'
       naux =  size(clu,2)-2;
	   ipermut = 1:naux;
    end
    clu_aux = zeros(size(clu,1),2 + size(spikes,1)) -1; % when update classes from clu, not selected elements go to cluster 0
    clu_aux(:,ipermut+2) = clu(:,(1:length(ipermut))+2);
    clu_aux(:,1:2) = clu(:,1:2);
    clu = clu_aux;
    clear clu_aux
elseif ~isempty(clu)
    naux = size(clu,2)-2;
    clu_aux = zeros(size(clu,1),2 + size(spikes,1)) -1; % when update classes from clu, not selected elements go to cluster 0
    clu_aux(:,(1:naux)+2) = clu(:,(1:naux)+2);
    clu_aux(:,1:2) = clu(:,1:2);
    clu = clu_aux;
    clear clu_aux
end
USER_DATA = get(handles.projections,'userdata');
USER_DATA{1} = handles.par;
USER_DATA{2} = spikes;
USER_DATA{3} = index;
USER_DATA{4} = clu;
USER_DATA{5} = tree;
USER_DATA{7} = inspk;
if exist('ipermut','var') && ~isempty(ipermut)
    USER_DATA{12} = ipermut;
end
USER_DATA{13} = forced;
USER_DATA{14} = forced;
USER_DATA{15} = rejected;  % the clusters numbers are sorted
USER_DATA{16} = rejected;  % the clusters numbers are sorted

% set(handles.min_clus_edit,'string',num2str(handles.par.min_clus));
% setappdata(handles.temperature_plot,'auto_sort_info',[]);

if  data_handler.with_gui_status
    [saved_gui_status, current_temp,auto_sort_info] = data_handler.get_gui_status();
%     setappdata(handles.temperature_plot,'auto_sort_info',auto_sort_info);
    clustering_results = zeros(length(classes),4);
    clustering_results(:,1) = repmat(current_temp,length(classes),1);
    for i=1:max(classes)
      clustering_results(classes==i,3)  = temp(i);
    end

    clustering_results(:,2) = classes'; % GUI classes
    clustering_results(:,4) = saved_gui_status;
    handles.undo = 1;

elseif data_handler.with_results
    current_temp = 1;

    clustering_results(:,1) = repmat(current_temp,length(classes),1); % GUI temperatures
    clustering_results(:,2) = classes';   % GUI classes
    if exist('temp','var')
        for i=1:max(classes)
            clustering_results(classes==i,3)  = temp(i);
        end
    else
        for i=1:max(classes)
            clustering_results(classes==i,3)  = current_temp;
        end
    end
    clustering_results(:,4) = classes';   % original classes
    handles.undo = 1;

else
    % Selects temperature.
    [clust_num,temp,auto_sort] = find_temp(tree, clu, handles.par);
    current_temp = max(temp);
    classes = zeros(1,size(clu,2)-2);
    for c =1: length(clust_num)
        aux = clu(temp(c),3:end) +1 == clust_num(c);
        classes(aux) = c;
    end

    if handles.par.permut == 'n'
        classes = [classes zeros(1,max(size(spikes,1)-size(clu,2)-2,0))];
    end
    Temp = [];
    % Classes should be consecutive numbers
    classes_names = nonzeros(sort(unique(classes)));
    for i= 1:length(classes_names)
       c = classes_names(i);
       if c~= i
           classes(classes == c) = i;
       end
       Temp(i) = temp(i);
    end
%     setappdata(handles.temperature_plot,'auto_sort_info',auto_sort);
    % definition of clustering_results
    clustering_results = [];
    clustering_results(:,1) = repmat(current_temp,length(classes),1); % GUI temperatures
    clustering_results(:,2) = classes';    % GUI classes

    for i=1:max(classes)
      clustering_results(classes==i,3) = temp(i);
      clustering_results(classes==i,4) = clust_num(i); % original classes
    end

    handles.undo = 0;
end

% if data_handler.with_results && ~ data_handler.with_spc
%     temp = -1;                               % This will work as a flag for not drawing the temperature diagram
%     cla(handles.temperature_plot,'reset');
% end

clustering_results(:,5) = repmat(handles.par.min_clus,length(classes),1);   % minimum number of clusters
USER_DATA{6} = classes(:)';
if exist('current_temp','var')
    USER_DATA{8} = current_temp;
else
    USER_DATA{8} = temp(1);
end
USER_DATA{10} = clustering_results;
USER_DATA{11} = clustering_results;
handles.force = 0;
handles.merge = 0;

handles.minclus = handles.par.min_clus;

set(handles.projections,'userdata',USER_DATA);

clear clustering_results classes rejected spikes
% mark clusters when new data is loaded
guidata(hObject, handles);    % this is need for plot the isi histograms

plot_spikes(handles);   % This function edits userdata
USER_DATA = get(handles.projections,'userdata');
set(handles.projections,'userdata',USER_DATA);

% if isfield(handles,'force_unforce_button') && (nnz(forced)>0)
% 	set(handles.force_unforce_button,'Value',1)
%     set(handles.force_unforce_button,'String','FORCED')
%     %set(handles.change_temperature_button,'enable','off');
% elseif isfield(handles,'force_unforce_button')
%     set(handles.force_unforce_button,'Value',0)
%     set(handles.force_unforce_button,'String','Force')
% end
% if isfield(handles,'edit_max_force_dist')
%     set(handles.edit_max_force_dist,'string',num2str(handles.par.template_sdnum));
% end

set(handles.file_name,'string',handles.par.file_name_to_show);

% Update handles structure
guidata(hObject, handles);


function file_name_Callback(hObject, eventdata, handles)
% hObject    handle to file_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of file_name as text
%        str2double(get(hObject,'String')) returns contents of file_name as a double


% --- Executes during object creation, after setting all properties.
function file_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_result_button.
function save_result_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_result_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
USER_DATA = get(handles.projections,'userdata');
classes = USER_DATA{6};
clustering_results = USER_DATA{10};
rejected = USER_DATA{15};

USER_DATA{16} =  rejected;     %update bk of rejected spikes
USER_DATA{15} = false(size(rejected));
set(handles.projections,'userdata',USER_DATA);
handles.setclus = 1;
handles.force = 0;
handles.merge = 0;

handles.undo = 0;
plot_spikes(handles)

USER_DATA = get(handles.projections,'userdata');
spikes = USER_DATA{2};
used_par = USER_DATA{1};
index = USER_DATA{3};
classes = USER_DATA{6};
gui_classes_data = USER_DATA{10};

% Classes should be consecutive numbers
classes_names = nonzeros(sort(unique(classes)));
for i= 1:length(classes_names)
   c = classes_names(i);
   if c~= i
       classes(classes == c) = i;
   end
end

% Saves clusters
cluster_class = zeros(size(spikes,1),2);
cluster_class(:,1) = classes(:);
cluster_class(:,2) = USER_DATA{3}';

outfile=['times_' used_par.nick_name];

par = struct;
par = update_parameters(par,used_par,'relevant');
par = update_parameters(par,used_par,'batch_plot');
par.sorting_date = datestr(now);

gui_status = struct();
gui_status.current_temp =  gui_classes_data(1,1);
gui_status.original_classes = gui_classes_data(1:end,4);

Temp = zeros(length(classes_names),1);
for i = 1:length(classes_names)
    Temp(i) = gui_classes_data(find(classes==i,1,'first'),3);
end
forced = USER_DATA{13};

var_list = 'cluster_class'',''par'',''gui_status'', ''forced'', ''Temp''';

if ~isempty(USER_DATA{7})
    inspk = USER_DATA{7};
    var_list = strcat(var_list , ' ,''inspk''');
end

if ~isempty(USER_DATA{12})
    ipermut = USER_DATA{12};
    var_list = strcat(var_list , ' ,''ipermut''');
end

% if developer_mode
% 	var_list = strcat(var_list , ' ,''rejected''');
% end

if isempty(handles.par.spikes_file)
    var_list = strcat(var_list , ' ,''spikes''');
else
    spikes_file = handles.par.spikes_file;
    var_list = strcat(var_list , ' ,''spikes_file''');
end

ver = '';
currentver = version;
if currentver(1) >= 7
    ver = ',''-v6''';
end
try
	eval(['save( ''' outfile ''',''' var_list '' ver ');']);
catch
	eval(['save( ''' outfile ''',''' var_list ',''-v7.3'');']);
end
if exist([handles.par.fnamespc '.dg_01.lab'],'file')
    movefile([handles.par.fnamespc '.dg_01.lab'], [handles.par.fnamesave '.dg_01.lab'], 'f');
    movefile([handles.par.fnamespc '.dg_01'], [handles.par.fnamesave '.dg_01'], 'f');
end

%Save figures
fig_names = {'figure','aux','aux1','aux2','aux3','aux4','aux5'};
file_names = {'','a','b','c','d','e','f'};

h_figs = get(0,'children');
for i=1:length(fig_names)
	h_fig =  findobj(h_figs,'tag',['wave_clus_' fig_names{i}]);
    new_file_name = ['fig2print_' outfile(7:end) file_names{i} '.png'];
	if ~isempty(h_fig)
        figure(h_fig); set(gcf, 'PaperUnits', 'inches', 'PaperType', 'A4', 'PaperPositionMode', 'auto','PaperOrientation','portrait');
        print(h_fig,'-dpng',new_file_name,'-r300');
    else
        if exist(new_file_name, 'file')==2
            delete(new_file_name);
        end
	end
end

set(hObject,'value',0);

% --- Executes on button press in plot_average_button.
function plot_average_button_Callback(hObject, eventdata, handles)
% hObject    handle to plot_average_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(gcbo,'value',1);
set(handles.plot_all_button,'value',0);
USER_DATA = get(handles.projections,'userdata');
cluster_results = USER_DATA{10};
handles.setclus = 1;
handles.force = 0;
handles.merge = 0;

handles.undo = 0;
handles.minclus = cluster_results(1,5);
plot_spikes(handles);
% Hint: get(hObject,'Value') returns toggle state of plot_average_button


% --- Executes on button press in Plot_all_projections_button.
function Plot_all_projections_button_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_all_projections_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
USER_DATA = get(handles.projections,'userdata');
par = USER_DATA{1};
if par.channels > 1
    Plot_amplitudes(handles)
else
    Plot_all_features(handles)
end


% --- Executes on button press in plot_all_button.
function plot_all_button_Callback(hObject, eventdata, handles)
% hObject    handle to plot_all_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'value',1);
set(handles.plot_average_button,'value',0);
USER_DATA = get(handles.projections,'userdata');
cluster_results = USER_DATA{10};
handles.setclus = 1;
handles.force = 0;
handles.merge = 0;
handles.undo = 0;
handles.minclus = cluster_results(1,5);
plot_spikes(handles);
% Hint: get(hObject,'Value') returns toggle state of plot_all_button


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
