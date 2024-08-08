xlinker_id = 0;
all_bunds_per_site = cell(50,1); 
LFB = zeros(50,1);
xdim = 40;
ydim = 10;
mask_file_name = '~/Morphometric_Actin_Masks/'; %change folder location to location of skeletonized image files

for xlink = [0.25 0.5 0.75 1 1.25 1.5 1.75 2] %crosslinker number density (per um^2)
        xlinker_id = xlinker_id + 1;

        %xl_mask_name = strcat(mask_file_name,num2str(xlink*xdim*ydim),'_XL_masks'); %Filaname for the skeletonized image
        %mask_list = getfn(xl_mask_name,'.tif');

    for file = 1:50

        name = strcat(pwd,'/','case_rest_xperiodic_100_11_',num2str(xlink),'_40_10_',num2str(file),'.txt');

        n_polymers = 100; %number of filaments
        n_links = 10; %number of links per filament
        npoints = 5; %number of points to be placed along a single filament link
        nframes = 1001; %number of simulation timesteps
        xdim = 40; %x-y dimensions of simulation domain
        ydim = 10;

        resolution = 0.025; %size of voxels used in local filament bundling calculation 
        tot_links = n_polymers * n_links;
        count = 0;

        unit_vector_x = zeros(tot_links,npoints+1);
        unit_vector_y = zeros(tot_links,npoints+1);

        bunds_per_site = []; %array to hold the number of filaments in each voxel

        links_res = readmatrix(name);
        links = links_res(all(~isnan(links_res),2),:); %matrix containing xy positions of both ends of filament links for all timeframes simulated 

        s = 1; %final timeframe of simulation

        min1 = (s * tot_links) - tot_links + 1; %returns the row positions of the filament links at the final simulation timeframe
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
            for j = 1:npoints+1 % This loop places points at uniform intervals along the filament links
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

        lowDimX = (resolution/2):resolution:(xdim-resolution/2);
        lowDimY = (resolution/2):resolution:(ydim-resolution/2);
        
        kernelDensity = cell(length(lowDimY),length(lowDimX)); %divide the system into voxels
       
        for m = 1:size(unit_vector_x,1)
            fil_id = unit_vector_x(m,end);
            for n = 1:size(unit_vector_x,2)-1
                a = unit_vector_x(m,n); b = unit_vector_y(m,n);
                lox=ceil(a/resolution);
                loy=ceil(b/resolution);
                
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
              
                kernelDensity{loy,lox}(end+1) = fil_id; %Finds filaments that are in each voxel.
          
            end
        end

        newKernelDensity = cellfun(@unique,kernelDensity,'UniformOutput',false); %returns new cell array with only unique filaments that are in each voxel 
        
        site_count = 10; %resolution*site_count*2 gives the side length of the square voxel. For this work, side length = 0.5 um.
        
        [actin_mask_row,actin_mask_col] = find(~cellfun(@isempty,newKernelDensity)); %return all voxels in skeletonized image that contain filaments
        
         %actin_mask = imread(mask_list{file});
         %[actin_mask_row,actin_mask_col] = find(actin_mask>0); %return all voxels in skeletonized image that contain filaments

           for non_zero_pixels = 1:length(actin_mask_col)

                ny = actin_mask_row(non_zero_pixels); %using voxel location from skeletonized image as array location of voxels in raw data
                nx = actin_mask_col(non_zero_pixels);

                lo_nx=nx - site_count;
                lo_ny=ny - site_count;
                hi_nx=nx + site_count;
                hi_ny=ny + site_count;

                if(lo_nx < 1)
                    lo_nx = 1;
                end
                if(lo_ny < 1)
                    lo_ny = 1;
                end
                if(hi_nx > length(lowDimX))
                    hi_nx = length(lowDimX);
                end
                if(hi_ny > length(lowDimY))
                    hi_ny = length(lowDimY);
                end

                bundle = []; %array to hold the unique filaments in the 0.5 um x 0.5 um region

                 for x = lo_nx:hi_nx
                    for y = lo_ny:hi_ny
                        fils_at_site = newKernelDensity{y,x};
                        for mem = 1:length(fils_at_site)
                            if(~ismember(fils_at_site(mem),bundle))
                                bundle(end+1) = fils_at_site(mem); %unique filaments within close proximity of filament-occupied voxel.
                            end
                        end
                    end
                 end

                 bundle = unique(bundle);

                 if(~isempty(bundle))&&(length(bundle)>3)
                    bunds_per_site(end+1,1) = length(bundle); %number of unique filaments in occupied 0.5 um x 0.5 um region
                 end
           end
         all_bunds_per_site{file,xlinker_id} = bunds_per_site; %return the LFB for each simulation trajectory for a given crosslinker density
         LFB(file,xlinker_id) = mean(bunds_per_site);
    end
   
 end


%%% Accessory functions to read data from all files for a given crosslinker density 
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