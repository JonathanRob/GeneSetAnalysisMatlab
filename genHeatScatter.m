function [] = genHeatScatter(sizedata,colordata,colnames,rownames,cluster_by,cluster_dist,col_map,sizebounds,colorbounds)
%genHeatScatter  Generate a heat-scatter plot with varing cell size/color.
%
% Usage:
%
%   genHeatScatter(sizedata, colordata, colnames, rownames, cluster_by, ...
%                  cluster_dist, col_map, sizebounds, colorbounds);
%
% Inputs:
%
%   sizedata     Matrix of data corresponding to the size of the markers.
%
%   colordata    Matrix of data corresponding to the colors of the markers.
%
%   colnames     A list of names of data columns. If not provided, column
%                index numbers will be used.
% 
%   rownames     A list of names of data rows. If not provided, row index
%                numbers will be used.
%
%   cluster_by  'none'  - the data will be plotted as provided (DEFAULT)
%               'color' - cluster/rearrange rows and columns based on their
%                         distances in COLORDATA
%               'size'  - cluster/rearrange rows and columns based on their
%                         distances in SIZEDATA
%               'both'  - cluster/rearrange rows and columns based on their
%                         distances in a matrix obtained by normalizing
%                         COLORDATA and SIZEDATA from 0 to 1 and averaging
%
%   cluster_dist   Distance metric to be used for clustering, ignored if
%                  CLUST_DIM is 'none'. Options are the same as those for
%                  distance in, e.g., PDIST ('euclidean', 'hamming', etc.).
%                  (DEFAULT = 'euclidean')
%
%   col_map      Colormap to which COLORDATA will be mapped.
%               (DEFAULT = custom 'redblue' colormap)
%
%   sizebounds  A 2-element vector indicating the min and max values in
%               SIZEDATA to which the point sizes will be scaled. Values at
%               or below the min value in SIZEBOUNDS will be set to the
%               minimum point size, whereas values at or above the max 
%               value in SIZEBOUNDS will be set to the maximum point size.
%               (Default = min and max of SIZEDATA)
%
%              *** NOTE: To adjust the min and max point sizes, change
%                        the SIZE_LIMS in the "additional settings"
%                        hard-coded below.
%
%   colorbounds A 2-element vector indicating the min and max values in
%               COLORDATA to which the point colors will be scaled. Values 
%               at or below the min value in COLORBOUNDS will be set to the
%               first color in the colormap, whereas values at or above the 
%               max value in COLORBOUNDS will be set to the last color.
%               (Default = min and max of COLORDATA)
%
%
% The following additional parameters can be modified within the function
% file itself:
%
%       size_lims           range of marker sizes
%       marker              shape of markers (square, circle, star, etc.)
%       marker_edge_color   color of marker edges (none, black, red, etc.)
%       grid_color          color of grid (none, black, gray, etc.)
%       linkage_method      algorithm used for hierarchical clustering
%
%
% Jonathan Robinson, 2019-03-07
%


%************************** ADDITIONAL SETTINGS ***************************
% Some additional minor settings to avoid too many function inputs.

size_lims = [20 200];  % adjusts the range of point sizes
marker = 'o';  % shape of marker to plot ('o' circle, 's' square, etc.)
marker_edge_color = 'k';  % e.g., 'k' for black, 'none' for no edges
grid_color = 'none';  % e.g., 'k' for black, 'none' for no grid
linkage_method = 'average';  % linkage algorithm for hierarchical clustering (see linkage function for more options)

%**************************************************************************


% handle input arguments and assign default values
if nargin < 3 || isempty(colnames)
    colnames = (1:size(colordata,2))';
end
if nargin < 4 || isempty(rownames)
    rownames = (1:size(colordata,1))';
end
if nargin < 5 || isempty(cluster_by)
    cluster_by = 'none';
end
if nargin < 6 || isempty(cluster_dist)
    cluster_dist = 'euclidean';
end
if nargin < 7 || isempty(col_map)
    col_map = 'redblue';
end
if nargin < 8
    sizebounds = [];
end
if nargin < 9
    colorbounds = [];
end

% verify input dimensions
if ~isequal(size(sizedata),size(colordata))
    error('"sizedata" and "colordata" must have the same dimensions.');
end
if length(colnames) ~= size(colordata,2)
    error('"colnames" must have the same numer of entries as columns in "sizedata" and "colordata".');
end
if length(rownames) ~= size(colordata,1)
    error('"rownames" must have the same numer of entries as rows in "sizedata" and "colordata".');
end

% determine data upon which to perform clustering 
switch cluster_by
    case 'size'
        clust_data = sizedata;
    case 'color'
        clust_data = colordata;
    case 'both'
        % normalize and average size and color data
        sizedata_norm = (sizedata - min(sizedata(:)))./var(sizedata(:),'omitnan');
        sizedata_norm(isnan(sizedata_norm)) = 0;  % remove NaNs
        colordata_norm = (colordata - min(colordata(:)))./var(colordata(:),'omitnan');
        colordata_norm(isnan(colordata_norm)) = 0;  % remove NaNs
        clust_data = (sizedata_norm + colordata_norm)./2;
    otherwise
        clust_data = [];
end
clust_data(isnan(clust_data)) = 0;  % not ideal, but one way to deal with NaNs

% perform hierarchical clustering to sort rows and columns (if specified)
if ~isempty(clust_data)
    L = linkage(clust_data,linkage_method,cluster_dist);
    row_ind = optimalleaforder(L,pdist(clust_data,cluster_dist));
    L = linkage(clust_data',linkage_method,cluster_dist);
    col_ind = optimalleaforder(L,pdist(clust_data',cluster_dist));
else
    row_ind = 1:size(sizedata,1);
    col_ind = 1:size(sizedata,2);
end

% reorder data matrix according to clustering results
sizedata_sort = sizedata(row_ind,col_ind);
colordata_sort = colordata(row_ind,col_ind);
sortrows = rownames(row_ind);
sortcols = colnames(col_ind);

% check if data is square matrix with identical row and column names
if (length(colnames) == length(rownames)) && all(strcmp(colnames,rownames))
    % flip data so the diagonal is from upper left to lower right
    sizedata_sort = fliplr(sizedata_sort);
    colordata_sort = fliplr(colordata_sort);
    sortcols = flipud(sortcols);
end

% scale data for plotting
if ~isempty(sizebounds)
    sizedata_sort(sizedata_sort < sizebounds(1)) = sizebounds(1);
    sizedata_sort(sizedata_sort > sizebounds(2)) = sizebounds(2);
else
    sizebounds = [min(sizedata_sort(:)),max(sizedata_sort(:))];
end
sizedata_sort = range(size_lims)/range(sizebounds)*(sizedata_sort - sizebounds(1))+size_lims(1);

% generate point locations
[ny,nx] = size(sizedata_sort);
x = repmat(1:nx,ny,1); x = x(:);
y = repmat(1:ny,nx,1)'; y = y(:);

% generate grid data
xgrid = repmat(2:nx,2,1)' - 0.5;
ygrid = repmat(2:ny,2,1)' - 0.5;


% generate plot
a = axes;
set(a,'YAxisLocation','Right');
set(a,'TickLength',[0 0],'XLim',[0.5 nx+0.5],'YLim',[0.5 ny+0.5]);
hold on

% draw grid lines
if ~isempty(grid_color) && ~strcmpi(grid_color,'none')
    line(xgrid,[0.5,ny+0.5],'color',grid_color);
    line([0.5,nx+0.5],ygrid,'color',grid_color);
end

% draw points
if strcmpi(marker_edge_color,'none')
    scatter(x,y,sizedata_sort(:),colordata_sort(:),marker,'filled');
else
    scatter(x,y,sizedata_sort(:),colordata_sort(:),marker,'filled','MarkerEdgeColor',marker_edge_color,'LineWidth',0.5);
end
set(gca,'XTick', 1:nx);
set(gca,'YTick', 1:ny);
set(gca,'XTickLabels',sortcols,'YTickLabels',sortrows);
set(gca,'XTickLabelRotation',90);

% bad form to have try/catch, but it works for now
try
    colormap(col_map);  % check if Matlab colormap exists
catch
    colormap(custom_cmap(col_map));  % check if custom colormap exists
end

if ~isempty(colorbounds)
    caxis(colorbounds);
end

% draw box around plot area (probably not necessary)
% xl = get(gca,'XLim');
% yl = get(gca,'YLim');
% plot(xl([1,1,2,2,1]),yl([1,2,2,1,1]),'w','LineWidth',1);

% Matlab will sometimes break up figures with polygons, which can add
% strange cuts through the figure. This command will prevent that.
set(gcf,'Renderer','painters');


