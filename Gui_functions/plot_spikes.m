function plot_spikes(handles)
set(handles.file_name,'string','Plotting...'); 
drawnow;
if exist('groot','builtin')
    if isprop(handles.projections,'GraphicsSmoothing')
        set(handles.projections,'GraphicsSmoothing','off');
    end
    try
        set(groot,'defaultfiguregraphicssmoothing','off');
        set(groot,'DefaultAxesFontSize',8)
    end
end

USER_DATA = get(handles.projections,'userdata');
par = USER_DATA{1};
spikes = USER_DATA{2};
classes = USER_DATA{6};
classes = classes(:)';
inspk = USER_DATA{7};
temp = USER_DATA{8};
ls = size(spikes,2);
minclus = handles.minclus;
clustering_results = USER_DATA{10};

% Extract spike features if needed
if get(handles.plot_average_button,'value') == 0
    if isempty(inspk) || (length(inspk)~=size(spikes,1))
        [inspk] = wave_features(spikes,handles);
        USER_DATA{7} = inspk;
    end
end

% Classes should be consecutive numbers
classes_names = sort(unique(classes));
classes_names = classes_names(classes_names>0);

% updates 'clustering_results_bk'
USER_DATA{11} = clustering_results; 

for i = 1:length(classes_names)
   c = classes_names(i);
   if c~= i
       classes(classes == c) = i;
   end
end

% Defines nclusters
cluster_sizes = zeros(1,33);

for i=1:size(cluster_sizes,2)
    cluster_sizes(i) = nnz(classes==i);
end

sizemin_clus = minclus;
clusn = find(cluster_sizes >= sizemin_clus);
nclusters = length(clusn);

% Defines classes
non_clustered = ones(1,size(spikes,1));
nclusters = 0;
for i = clusn
    class_temp = find(classes == i);
    if (length(class_temp) >= sizemin_clus)
        nclusters = nclusters+1;
        eval(['class' num2str(nclusters) '= class_temp;'])
        non_clustered(class_temp) = 0;
    end
end
rejected = USER_DATA{15};
class0 = find(non_clustered & ~rejected);
clear non_clustered rejected

% Redefines classes
classes = zeros(size(spikes,1),1);
for i = 0:nclusters
    if ~ (isempty(class0) && i==0)
        eval(['classes(class' num2str(i) ') = ' num2str(i) ';']);
    end
end

% Saves new classes
USER_DATA{6} = classes;

clustering_results(:,1) = temp;    % GUI temperature
clustering_results(:,5) = minclus; % GUI minimum cluster

% update new classes
clustering_results(:,2) = classes;

% Updates clustering_results in USER_DATA
USER_DATA{10} = clustering_results; 

for i=20:52
    USER_DATA{i} = [];
end

for i=4:33
    eval(['par.fix' num2str(i) '=0;']);
end

set(handles.projections,'userdata',USER_DATA)

cla(handles.projections)
hold(handles.projections,'on')

% Plot clusters
ylimit = [];
colors = [[0.0 0.0 1.0];[1.0 0.0 0.0];[0.0 0.5 0.0];[0.620690 0.0 0.0];[0.413793 0.0 0.758621];[0.965517 0.517241 0.034483];
    [0.448276 0.379310 0.241379];[1.0 0.103448 0.724138];[0.545 0.545 0.545];[0.586207 0.827586 0.310345];
    [0.965517 0.620690 0.862069];[0.620690 0.758621 1.]]; 
maxc = size(colors,1);
forced = USER_DATA{13};
figs_num = 6;
opened_figs = cell(1,figs_num);
spikes_num = [];
for i = 0:nclusters
    tmpy = [];     % as a flag to don't make the same vector twice
    if ~ (isempty(class0) && i==0)
        % PLOTS SPIKES OR PROJECTIONS
        class_i = eval(['class' num2str(i)]);
        sup_spikes = length(class_i);
        max_spikes = min(sup_spikes, par.max_spikes_plot);
        permut = randperm(sup_spikes);
        permut = permut(1:max_spikes);
        xlim(handles.projections,'manual');
        if get(handles.plot_all_button,'value') ==1
            % optimizing for speed:
            tmpy=spikes(class_i(permut),:);
            tmpn=size(tmpy,1);
            tmpx=repmat([1:ls NaN]',1,tmpn);
            tmpx=reshape(tmpx,numel(tmpx),1);
            tmpy=[tmpy'; repmat(NaN,1,tmpn)];
            tmpy=reshape(tmpy,numel(tmpy),1);
            line(tmpx,tmpy,'color',colors(mod(i-1,maxc)+1,:)*(i~=0),'Parent',handles.projections,'Visible','off');
			xlim(handles.projections, [1 ls])
        elseif get(handles.plot_average_button,'value') ==1
            av  = mean(spikes(class_i,:));
            plot(handles.projections,1:ls,av,'color',colors(mod(i-1,maxc)+1,:)*(i~=0),'linewidth',2);
            xlim(handles.projections,[1 ls])
%             legend
%             num_cluster = [[];{strcat('cluster ',num2str(i))}];
%             legend(num_cluster(1:i))
        else
            plot(handles.projections,inspk(class_i,1),inspk(class_i,2),'.','Color',colors(mod(i-1,maxc)+1,:)*(i~=0),'markersize',.5);
            axis(handles.projections,'auto');
        end
        
        eval(['aux=num2str(length(class' num2str(i) '));']);
        spikes_num = [spikes_num;{strcat('Cluster  ',num2str(i), ':  #  ', aux)}];
%         set(handles.listbox1,'string',['Cluster ' num2str(i) ':  # ' aux ' (' num2str(nnz(clustering_results(:,2)==i & ~forced(:))) ')']);
        set(handles.listbox1,'string',spikes_num)
        drawnow
    end
end

set(handles.file_name,'string', par.file_name_to_show);

set(allchild(handles.projections),'Visible','on')

% Resize axis
if ~isempty(ylimit)
    ymin = min(ylimit(:,1));
    ymax = max(ylimit(:,2));
    ylim(handles.spikes0,[ymin ymax]);
    if get(handles.plot_average_button,'value') ==1
        ylim(handles.projections,[ymin ymax]);
    end
    linkaxes(ax_v,'xy');     % drawnow inside
    ylim(ax_v(1),[ymin ymax]);
else
    drawnow
end

% for i =1:figs_num
%     if ~isempty(opened_figs{i})  
%         set(opened_figs{i},'units','normalized','outerposition',[0 0 1 1])
%         set(opened_figs{i},'Visible', 'on'); 
%     end
% end

if exist('groot','builtin')
    set(groot,'defaultfiguregraphicssmoothing','remove')
    set(groot,'DefaultAxesFontSize','remove')
end