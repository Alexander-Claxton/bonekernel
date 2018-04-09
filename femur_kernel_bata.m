function [xx, bata, bone_adf, total_adf] = ...
    femur_kernel_bata(image, mask, microns_per_pixel, dist, bone_thresh, ...
        num_pts, bandwidth, show)
%
% Compute bone area, total area, and their ratio as a function of distance
% from the boundary from a cross section of a CT-scanned femur bone.
%
% Input Arguments:
% 
% image = 2D matrix, grayscale CT-scan cross section image, cropped to 
%   remove any annotation.
%
% mask = 2D matrix, logical mask indicating if each pixel should be included
%   in the analysis (1) or not (0)
%
% microns_per_pixel = Unit conversion factor for pixels in image
%
% dist = 2D matrix, distance to the bone boundary for all pixels within the
%   bone footprint, NaN for all pixels outside the bone footprint, [micron]
%
% bone_thresh = Scalar, threshold intensity value above which pixels are
%   classified as "bone", below which pixels are classified as "not bone".
%
% num_pts = (optional) Scalar, integer. Number of points at which to
%   estimate the outputs, default = 100.
%
% bandwidth = (optional) Scalar, >0. Bandwidth for kernel density PDF
%   estimation, in [microns], default = 50.
%
% show = (optional) Scalar, true/false, flag turning the plotting routine
%   on (true) or off (false), default = false.
%
% Output Arguments:
%
% xx = Vector, double, length = num_pts. Distance to the outer boundary
%   of the bone, serves as the independant variable (e.g. x-axis) for bata,
%   bone_adf, and total_adf.
%
% bata = Vector, double, length = num_pts. Bone area / total area,
%   dimensionless, the primary output of this function.
%
% bone_adf = Vector, double, length = num_pts. The "area density function"
%   for bone pixels, computed by estimating the probability density
%   function of distance from the boundary for all bone pixels, and
%   multiplying by the total area of bone. Units are pixel^2/pixel = pixel.
%
% total_adf = Vector, double, length = num_pts. Same as bone_adf except the
%   calculation includes all points within the bone (pore space as well as
%   bone)
%
% %

warning('off', 'images:initSize:adjustingMag');

%% parse inputs

% set defaults
narginchk(4, 8);
if nargin < 5 || isempty(bone_thresh); bone_thresh = 0; end
if nargin < 6 || isempty(num_pts);     num_pts = 100;   end
if nargin < 7 || isempty(bandwidth);   bandwidth = 50; end
if nargin < 8 || isempty(show);        show = 1;        end

validateattributes( image,             {'numeric'}, {'2d', 'real'}             );
validateattributes( mask,              {'numeric', 'logical'}, {'2d', 'binary'});
validateattributes( microns_per_pixel, {'numeric'}, {'scalar', 'positive'}     );
validateattributes( dist,              {'numeric'}, {'2d', 'real'}             );
validateattributes( bone_thresh,       {'numeric'}, {'scalar'}                 );
validateattributes( num_pts,           {'numeric'}, {'scalar', 'integer'}      );
validateattributes( bandwidth,         {'numeric'}, {'scalar', '>', 0}         );
validateattributes( show,              {'numeric', 'logical'}, {'binary'}      );

assert(all(size(image) == size(mask)));
assert(all(size(image) == size(dist)));

%% main 

roi = mask & ~isnan(dist);
bone = roi & image>bone_thresh;

% estimate probability density functions for bony pixels and all pixels
xx_min = min(dist(roi));
xx_max = max(dist(roi));
xx = linspace(xx_min, xx_max, num_pts);

get_pdf = @(d) ksdensity(d, xx, ...
                         'kernel', 'epanechnikov', ...
                         'function', 'pdf', ...
                         'bandwidth', bandwidth);
                     
total_pdf = get_pdf(dist(roi)); % [1/micron]
bone_pdf = get_pdf(dist(bone)); % [1/micron]

% convert probability density to "area density" in [micron^2/micron]
%... if area = A(x), then area density = (d/dx) A(x)
total_adf = total_pdf*sum(roi(:))*microns_per_pixel^2;
bone_adf = bone_pdf*sum(bone(:))*microns_per_pixel^2;

% estimate (bone area)/(total area), dimensionless
bata = bone_adf./total_adf;

%% plot

if show
   
    % plot bone area / total area results
    figure
    set(gcf, 'Name', 'Bone Area Results', ...
        'Units', 'Normalized');
    pos = get(gcf, 'Position');
    pos(2) = 0;
    pos(4) = 1;
    set(gcf, 'Position', pos);

    subplot(2,1,1)
    plot(xx, total_adf, xx, bone_adf)
    xlabel('distance from boundary [\mu{}m]');
    ylabel('area density [\mu{}m^{2}/\mu{}m]');
    legend('total', 'bone')
    
    subplot(2,1,2)
    plot(xx, bata)
    xlabel('distance from boundary [\mu{}m]');
    ylabel('bone area / total area [1]');

end
