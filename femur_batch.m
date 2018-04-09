function [] = femur_batch(input_list)
% function [] = femur_batch(input_list)
%
% Batch processing of femur cross-section images. Input parameters are read
% from a csv file, with comma-separated parameters for each image in one
% row. The expected format is one header line, followed by parameters:
%
%   file, scale, threshold, radius, bandwidth
%   
% Where:
%   file = image file name
%   scale = image to world coordinate scaling, [microns/pixel]
%   threshold = minimum threshold for image data, [intensity]
%   radius  = morphological filter structuring element radius, [microns]
%   bandwidth = kernel smoother bandwidth, [microns]
% %

% constant parameters
num_pts = 1000; 
sigma_threshold = 0.9; 

% sanity check
validateattributes(input_list, {'char'}, {'vector'});
assert(exist(input_list, 'file') == 2);

% loop over lines in input parameter list
fp = fopen(input_list, 'r');
fgetl(fp); % skip header
while 1
   
    %% read input data
    
    % read line, exit if EOF
    line = fgetl(fp);
    if line == -1; break; end
    
    % extract parameters for this image
    line_parts = strsplit(line, ',');
    image_file = line_parts{1};
    microns_per_pixel = str2double(line_parts{2});
    bone_threshold = str2double(line_parts{3});
    disk_radius = str2double(line_parts{4});
    bandwidth = str2double(line_parts{5});
    
    [~, base] = fileparts(image_file);
    mkdir(base);
    
    %% prepare image: read, convert to grayscale, and threshold
    
    image = imread(image_file);
    if ndims(image) == 3
        image = rgb2gray(image);
    end    
    image(image<bone_threshold) = 0;
    image = double(image);
    
    %% get local coordinates in microns 
    
    [dist, xp, yp] = ...
        femur_local_coord(image, microns_per_pixel, disk_radius, sigma_threshold, 1);
        
    saveas(gcf, fullfile(base, [base, '_coord.fig']));
    close(gcf);
    
    %% compute bone/total area using kernel estimation
    
    %...for whole cross section
    mask = ones(size(image));
    [dd_all, bata_all, bone_adf_all, total_adf_all] = ...
        femur_kernel_bata(image, mask, microns_per_pixel, dist, ...
            bone_threshold, num_pts, bandwidth, 1); %#ok       
    saveas(gcf, fullfile(base, [base, '_bata_all.fig']));
    close(gcf);
    
    %...for each half
    mask = yp>=0;
    [dd_n, bata_n, bone_adf_n, total_adf_n] = ...
        femur_kernel_bata(image, mask, microns_per_pixel, dist, ...
            bone_threshold, num_pts, bandwidth, 1); %#ok
    saveas(gcf, fullfile(base, [base, '_bata_n.fig']));
    close(gcf);
        
    mask = yp<=0;
    [dd_s, bata_s, bone_adf_s, total_adf_s] = ...
        femur_kernel_bata(image, mask, microns_per_pixel, dist, ...
            bone_threshold, num_pts, bandwidth, 1); %#ok
    saveas(gcf, fullfile(base, [base, '_bata_s.fig']));
    close(gcf);
    
    mask = xp>=0;
    [dd_e, bata_e, bone_adf_e, total_adf_e] = ...
        femur_kernel_bata(image, mask, microns_per_pixel, dist, ...
            bone_threshold, num_pts, bandwidth, 1); %#ok
    saveas(gcf, fullfile(base, [base, '_bata_e.fig']));
    close(gcf);
    
    mask = xp<=0;
    [dd_w, bata_w, bone_adf_w, total_adf_w] = ...
        femur_kernel_bata(image, mask, microns_per_pixel, dist, ...
            bone_threshold, num_pts, bandwidth, 1); %#ok
    saveas(gcf, fullfile(base, [base, '_bata_w.fig']));
    close(gcf);
                
    %...for each quadrant
    mask = xp>=0 & yp>=0;
    [dd_ne, bata_ne, bone_adf_ne, total_adf_ne] = ...
        femur_kernel_bata(image, mask, microns_per_pixel, dist, ...
            bone_threshold, num_pts, bandwidth, 1); %#ok
    saveas(gcf, fullfile(base, [base, '_bata_ne.fig']));
    close(gcf);
    
    mask = xp>=0 & yp<=0;
    [dd_se, bata_se, bone_adf_se, total_adf_se] = ...
        femur_kernel_bata(image, mask, microns_per_pixel, dist, ...
            bone_threshold, num_pts, bandwidth, 1); %#ok
    saveas(gcf, fullfile(base, [base, '_bata_se.fig']));
    close(gcf);
    
    mask = xp<=0 & yp>=0;
    [dd_nw, bata_nw, bone_adf_nw, total_adf_nw] = ...
        femur_kernel_bata(image, mask, microns_per_pixel, dist, ...
            bone_threshold, num_pts, bandwidth, 1); %#ok
    saveas(gcf, fullfile(base, [base, '_bata_nw.fig']));
    close(gcf);
    
    mask = xp<=0 & yp<=0;
    [dd_sw, bata_sw, bone_adf_sw, total_adf_sw] = ...
        femur_kernel_bata(image, mask, microns_per_pixel, dist, ...
            bone_threshold, num_pts, bandwidth, 1); %#ok       
    saveas(gcf, fullfile(base, [base, '_bata_sw.fig']));
    close(gcf);
    
    %% save results
    
    save(fullfile(base, [base, '.mat']));

end
fclose(fp);

