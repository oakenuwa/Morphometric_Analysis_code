xlinker_id = 0;
all_median_per_site = zeros(50,1);
all_occup_per_site = zeros(50,1);

for xlink = [0 0.25 0.5 0.75 1 1.25 1.5 1.75 2]
        xlinker_id = xlinker_id + 1;

    for file = 1:50

        name = strcat(pwd,'/','case_rest_xperiodic_100_11_',num2str(xlink),'_40_10_',num2str(file),'.txt');
        
        n_polymers = 100;
        n_links = 10;
        npoints = 1000;

        nframes = 1001;
        xdim = 40;
        ydim = 10;

        tot_links = n_polymers * n_links;
        count = 0;

        unit_vector_x = zeros(tot_links,npoints+1);
        unit_vector_y = zeros(tot_links,npoints+1);


        links_res = readmatrix(name);
        links = links_res(all(~isnan(links_res),2),:);
               
        s = 1;
        
        min1 = (s * tot_links) - tot_links + 1;
        max1 = s * tot_links;

        for i = min1:max1
            a = links(i,3);
            b = links(i,4);
            c = links(i,1);
            d = links(i,2);
            e = links(i,5);
            dist = sqrt((a*a)+(b*b));
            fdist = dist/npoints;
            t = sqrt((fdist*fdist)/((a*a)+(b*b)));
            num = i-((s-1)*tot_links);
            for j = 1:npoints+1
                px = c+(t*(j-1)*a)+(xdim/2);
                py = d+(t*(j-1)*b)+(ydim/2);
                if(px < 0)
                    px = px + xdim;
                end
                if(py < 0)
                    py = py + ydim;
                end
                if(px > xdim)
                    px = px-xdim;
                end
                if(py > ydim)
                    py = py-ydim;
                end
                unit_vector_x(num,j) = px;
                unit_vector_y(num,j) = py;
            end
            unit_vector_x(num,j+1) = e;
            unit_vector_y(num,j+1) = e;
        end
        
        ntpoints = 1;
        xresolution = 0.02122; % x dimension of voxel size chosen to match pixel size of pseudo fluorescence images
        yresolution = 0.02118; % y dimension of voxel size chosen to match pixel size of pseudo fluorescence images

        lowDimX = (xresolution/2):xresolution:(xdim-xresolution/2);
        lowDimY = (yresolution/2):yresolution:(ydim-yresolution/2);
        kernelDensity = cell(length(lowDimY),length(lowDimX));

            for m = 1:size(unit_vector_x,1)
                fil_id = unit_vector_x(m,end);
                for n = 1:size(unit_vector_x,2)-1

                    a = unit_vector_x(m,n); b = unit_vector_y(m,n);

                    lox=ceil(a/xresolution);
                    loy=ceil(b/yresolution);

                    if(lox < 1)
                        lox = 1;
                    end
                    if(loy < 1)
                        loy = 1;
                    end
                    if(lox > length(lowDimX))
                        lox = length(lowDimX);
                    end
                    if(loy > length(lowDimY))
                        loy = length(lowDimY);
                    end

                      kernelDensity{loy,lox}(end+1) = fil_id;

                end
            end

            newKernelDensity = cellfun(@length,kernelDensity,'UniformOutput',false);
            bmap = cell2mat(newKernelDensity);
            bmap=bmap~=0;

            distmap = double(bwdist(bmap,'euclidean'));
            distmap(distmap == 0) = NaN;
            med_dist = median(distmap,'all','omitnan');

        all_median_per_site(file,xlinker_id) = med_dist*0.0212; %Distance for each for each simulation trajectory for a given crosslinker density
        all_occup_per_site(file,xlinker_id) = nnz(bmap)/(size(bmap,1)*size(bmap,2));%Occupancy for each simulation trajectory for a given crosslinker density
     end
end


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