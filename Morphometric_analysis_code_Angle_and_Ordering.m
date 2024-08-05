xlinker_id = 0;

for xlink = [0 0.25 0.5 0.75 1 1.25 1.5 1.75 2]
        xlinker_id = xlinker_id + 1;
        name = strcat('case_rest_xperiodic_100_11_',num2str(xlink),'_');
        list = get_directory_names(pwd,name);

        n_polymers = 100;
        n_links = 10;

        xdim = 40;
        ydim = 10;

        tot_links = n_polymers * n_links;

        angles_in_network = cell(length(list),1);
        mean_std_angles_in_network = cell(length(list),1);
        order_parameter_in_network = cell(length(list),1);

    for file = 1:length(list)

                txt = strcat(pwd,'/',list{file});
                fn = getfn(txt, 'links.txt');
                if(length(fn)>1)
                    s = dir(fn{2});
                    if(s.bytes==0)
                        continue;
                    else
                        links_res = readmatrix(fn{2});
                    end
                else
                    links_res = readmatrix(fn{1});
                end
                links_res = links_res(:,1:end);  
                links = links_res(all(~isnan(links_res),2),:);
                
                angles = zeros(tot_links,1);
                order_param_per_link = zeros(tot_links,1);

            for s = length(links)/tot_links
                min1 = (s * tot_links) - tot_links + 1;
                max1 = s * tot_links;
                count = 1;

                for i = min1:max1
                    a = links(i,3);
                    b = links(i,4);
                    theta = atand(b/a);
                    angles(count,1) = theta;
                    count = count+1;
                end
            end

            mean_angle = mean(angles); % returns mean angle of filaments in simulated network
            angle_variation = std(angles); % returns angular variation of filaments in simulated network
            
            for s = length(links)/tot_links
                min1 = (s * tot_links) - tot_links + 1;
                max1 = s * tot_links;
                count = 1;

                for i = min1:max1
                    a = links(i,3);
                    b = links(i,4);
                    theta = atand(b/a);
                    order_param_per_link(count,1) = (2*cosd(theta-mean_angle)*cosd(theta-mean_angle)) - 1;
                    count = count+1;
                end
            end
            
            order_param = mean(order_param_per_link); % returns the order parameter of simulated network

            angles_in_network{file,xlinker_id} = angles;
            mean_std_angles_in_network{file,xlinker_id} = [mean_angle,angle_variation];
            order_parameter_in_network{file,xlinker_id} = order_param;
    end 
end

parallelness = get_parallelness_value(angles_in_network); %parallelness calculation

function filenames = getfn(mydir, pattern)
%GETFN Get filenames in directory and subdirectories.
%
%   FILENAMES = GETFN(MYDIR, PATTERN)
%
% Example: Get all files that end with 'txt' in the current directory and
%          all subdirectories 
%
%    fn = getfn(pwd, 'txt$')
%
%   Thorsten.Hansen@psychol.uni-giessen.de  2016-07-06
if nargin == 0
  mydir = pwd;
end
% computes common variable FILENAMES: get all files in MYDIR and
% recursively traverses subdirectories to get all files in these
% subdirectories: 
getfnrec(mydir) 
% if PATTERN is given, select only those files that match the PATTERN:                 
if nargin > 1 
  idx = ~cellfun(@isempty, regexp(filenames, pattern));
  filenames = filenames(idx);
end
    function getfnrec(mydir)
    % nested function, works on common variable FILENAMES
    % recursively traverses subdirectories and returns filenames
    % with path relative to the top level directory
      d = dir(mydir);
      filenames = {d(~[d.isdir]).name};
      filenames = strcat(mydir, filesep, filenames); 
      dirnames = {d([d.isdir]).name};
      dirnames = setdiff(dirnames, {'.', '..'});  
      for i = 1:numel(dirnames)
        fulldirname = [mydir filesep dirnames{i}];
        filenames = [filenames, getfn(fulldirname)];
      end  
    end % nested function
end

function [ dirList ] = get_directory_names(dir_name,match)

    %get_directory_names; this function outputs a cell with directory names (as
    %strings), given a certain dir name (string)
    %from: http://stackoverflow.com/questions/8748976/list-the-subfolders-
    %in-a-folder-matlab-only-subfolders-not-files

    dd = dir(dir_name);
    isub = [dd(:).isdir]; %# returns logical vector
    dirList = {dd(isub).name}';
    dirList(ismember(dirList,{'.','..'})) = [];
    idx = find(contains(dirList,match));
    for i = 1:length(idx)
        dirList2{i} = dirList{idx(i)};
    end
    dirList = dirList2;
end

function para = get_parallelness_value(all_angles_in_network) %function to calculate the parallelness from the angles of filament links
parallelness = cell(size(all_angles_in_network,2),1);
for i = 1:size(all_angles_in_network,2)

parallel_per_traj = zeros(size(all_angles_in_network,1),1);

    for k = 1:size(all_angles_in_network,1)
        p = all_angles_in_network{k,i}+90;
        [bincounts, binedges] = histcounts(p,'BinWidth',0.5);
        n0 = 0; n45 = 0; n90 = 0; n135 = 0;

        for edge = 1:length(binedges)-1
            if(binedges(edge+1)<23 || (binedges(edge+1)>157 && binedges(edge+1)<=180))
                n0 = n0 + bincounts(edge);
            end
            if(binedges(edge+1)>=23 && binedges(edge+1)<=67.5)
                n45 = n45 + bincounts(edge);
            end
            if(binedges(edge+1)>=68 && binedges(edge+1)<=112.5)
                n90 = n90+ bincounts(edge);
            end
            if(binedges(edge+1)>=113 && binedges(edge+1)<=157)
                n135 = n135 + bincounts(edge);
            end
        end
            total = n0 + n45 + n90 + n135;
            a1 = abs(n0-n90);
            a2 = abs(n45-n135);
            para = (a1+a2)/total;
            parallel_per_traj(k,1) = para;
    end

    parallelness{i,1} = parallel_per_traj;

end
    para = parallelness;
end