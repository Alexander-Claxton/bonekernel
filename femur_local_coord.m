function [dist, x_principal, y_principal] = femur_local_coord(image, microns_per_pixel, disk_radius, sigma_threshold, show)
% function [dist, x_principal, y_principal] = femur_local_coord(image, disk_radius, sigma_threshold, show)
%
% Return image data for femur CT-scan cross sections in local coordinate
% systems, including:
%   - distance from outer boundary (1D)
%   - principal coordinates (major and minor axes of the bone footprint)
%
% All distance units are in pixels (image coordinates)
% 
% About segmentation:
%   DETAILS NEEDED 
%
% About distance coordinate:
%   DETAILS NEEDED 
%
% About principal coordinates: 
%   DETAILS NEEDED 
%
% Arguments:
% 
%   image = 2D matrix, grayscale CT-scan cross section image, cropped to 
%       remove any annotation.
%
%   microns_per_pixel = Scalar, unit conversion factor
%
%   disk_radius = (optional) Scalar, integer. Radius of structuring element
%       used for image segmentation by morphological closing, default = 20.
%
%   sigma_threshold = (optional) Scalar. Parameter used to exlude rotation
%       for bone images that are nearly circular. In these cases, the
%       initial orientation is almost certainly better than fitting the
%       noise in the bone footprint. The maximum ratio of minor principle
%       axis standard deviation to major principle axis standard deviation.
%       Must be <= 1 and >= 0. Default is 0.9.
%
%   show = (Optional) Scalar, binary, flag to generate plots (1) or not (0).
%       Default is 1.
% 
%   dist = 2D matrix, size(img), distance to the bone boundary for all
%       pixels within the bone footprint, NaN for all pixels outside the
%       bone footprint.
%
%   x_principal, y_principal = 2D matrix, size(image), location of each
%       image pixel in principal coordinate system.
% %

%% init inputs 

% set defaults, if needed
narginchk(1,5);
if nargin<3; disk_radius = 20;      end
if nargin<4; sigma_threshold = 0.9; end
if nargin<5; show = 1;              end

% sanity check
validateattributes(image, {'numeric'}, {'2d', 'real'});
validateattributes(microns_per_pixel, {'numeric'}, {'scalar'});
validateattributes(sigma_threshold, {'numeric'}, {'scalar', '<=', 1, '>=', 0});
validateattributes(disk_radius, {'numeric'}, {'scalar', 'integer'});

%% main

% get femur roi

%...convert disk to pixels
disk_radius = ceil(disk_radius/microns_per_pixel);

%...pad to avoid imclose edge effects
image = padarray(image, disk_radius*[1,1], 0, 'both'); 

%...morphological filters
disk = strel('disk', disk_radius);
roi = imclose(image>0, disk); % close (dilate, then erode)
roi = imfill(roi, 'holes');

%...remove pad
image = image(disk_radius+1:end-disk_radius, disk_radius+1:end-disk_radius); 
roi = roi(disk_radius+1:end-disk_radius, disk_radius+1:end-disk_radius); 

%...keep only 1 largest object in ROI (the bone, one hopes)
[label, numobj] = bwlabel(roi);
numelem = zeros(numobj,1);

for ii = 1:numobj
    numelem(ii) = sum(label(:) == ii);
end

[~, biggest] = max(numelem);
roi = label==biggest;

% compute distace from the bone boundary
bnd = bwperim(roi, 8);
dist = bwdist(bnd); % kept for plotting
dist(~roi) = NaN;

% get principal coordinates

% ...bone centroid in image coordinates (x = col, y = row, reversed)
[x_img, y_img] = meshgrid(1:size(image,2), size(image,1):-1:1);
x_bone_img = x_img(roi);
y_bone_img = y_img(roi);
centroid_x_bone_img = mean(x_bone_img);
centroid_y_bone_img = mean(y_bone_img);

%...principal axes 
x_img_centered = x_img-centroid_x_bone_img;
y_img_centered = y_img-centroid_y_bone_img;
[eig_vec, eig_val] = eig(cov([x_img_centered(roi), y_img_centered(roi)]));

sigma = sqrt(sort([eig_val(1,1), eig_val(2,2)]));
if sigma(1)/sigma(2) > sigma_threshold
    %...don't rotate if std() of principal axes is nearly equal (i.e. bone is circular)
    x_principal3 = x_img_centered;
    y_principal3 = y_img_centered;

    theta1 = -45;
    rot1 = [cosd(theta1), -sind(theta1); sind(theta1), cosd(theta1)];
    tmp = [x_img_centered(:), y_img_centered(:)]*rot1;
    x_principal1 = reshape(tmp(:,1), size(image));
    y_principal1 = reshape(tmp(:,2), size(image));
    
else
    %...select minimum rotation from image coordinates to principal coordinates,
    eig_vec_y_all = [eig_vec(2,:), -eig_vec(2,:)]; % get y component of all 4 principal basis vectors
    theta = acosd(max(eig_vec_y_all));
    theta1 = theta-45;
    %...generate principal coordinate matrices
    rot = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
    %rot1 = [cosd(theta1), -sind(theta1); sind(theta1), cosd(theta1)];
    
    tmp = [x_img_centered(:), y_img_centered(:)]*rot;
    x_principal3 = reshape(tmp(:,1), size(image));
    y_principal3 = reshape(tmp(:,2), size(image));
    
    rot1 = [cosd(theta1), -sind(theta1); sind(theta1), cosd(theta1)];
    tmp1 = [x_img_centered(:), y_img_centered(:)]*rot1;
    x_principal1 = reshape(tmp1(:,1), size(image));
    y_principal1 = reshape(tmp1(:,2), size(image));
end

% convert units from pixels to microns
dist = dist*microns_per_pixel;
x_principal = x_principal1*microns_per_pixel;
y_principal = y_principal1*microns_per_pixel;

%% show

if show
    
    set(0,'DefaultSurfaceEdgeColor','none');
    
    figure
    set(gcf, 'Name', 'femur_local_coord', ...
        'Units', 'Normalized');
    pos = get(gcf, 'Position');
    pos(1) = 0;
    pos(3) = 1;
    set(gcf, 'Position', pos);

    subplot(1,5,1)    
    pcolor(x_img, y_img, image);
    xlabel('x [pixel]');
    ylabel('y [pixel]');
    colormap(gca, gray);
    colorbar
    axis equal
    grid on
    title('original image');    
    
    subplot(1,5,2)  
    pcolor(x_img, y_img, double(roi));
    xlabel('x [pixel]');
    ylabel('y [pixel]');    
    colormap(gca, gray);
    colorbar
    axis equal
    grid on
    title('bone footprint');  
    
    subplot(1,5,3)    
    pcolor(x_img, y_img, dist); 
    xlabel('x [pixel]');
    ylabel('y [pixel]');    
    colormap(gca, 'default')
    colorbar
    axis equal
    grid on
    title('distance to bone boundary')

    subplot(1,5,4)    
    bone_image = double(image); 
    bone_image(~roi) = NaN;
    pcolor(x_principal3, y_principal3, bone_image); 
    xlabel('x [\mu{}m]');
    ylabel('y [\mu{}m]');    
    colorbar
    axis equal
    grid on
    title('principal coordinates');
    
    subplot(1,5,5)    
    bone_image = double(image); 
    bone_image(~roi) = NaN;
    pcolor(x_principal1, y_principal1, bone_image); 
    xlabel('x [\mu{}m]');
    ylabel('y [\mu{}m]');    
    colorbar
    axis equal
    grid on
    title('principal coordinates rotated');
    
end