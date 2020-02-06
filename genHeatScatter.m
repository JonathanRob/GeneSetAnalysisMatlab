function [] = genHeatScatter(sizedata,colordata,varargin)
%genHeatScatter  Generate a heatscatter plot with varing marker size/color
%
% Usage:
%
%   genHeatScatter(sizedata, colordata, ...);
%
% Input:
%
%   sizedata     Matrix of data corresponding to the size of the markers.
%
%   colordata    Matrix of data corresponding to the colors of the markers.
%
% Optional Settings:
%
%   'colNames'     A list of names of data columns.
%                  (DEFAULT = column index numbers)
% 
%   'rowNames'     A list of names of data rows.
%                  (DEFAULT = row index numbers)
%
%   'clusterDim'  'none'  - the data will be plotted as provided
%                 'color' - cluster/rearrange rows and columns based on
%                           their distances in COLORDATA
%                 'size'  - cluster/rearrange rows and columns based on
%                           their distances in SIZEDATA
%                 'both'  - (DEFAULT) cluster/rearrange rows and columns
%                           based on their distances in a matrix obtained
%                           by mean-centering and scaling COLORDATA and
%                           SIZEDATA to unit variance
%
%   'clusterDist'   Distance metric to be used for clustering, ignored if
%                   'clusterDim' is 'none'. Options are the same as those
%                   for distance in PDIST ('euclidean', 'hamming', etc.).
%                   (DEFAULT = 'euclidean')
%
%   'linkage'       String specifying the linkage method to be used for
%                   hierarchical clustering (see built-in LINKAGE function
%                   function for options).
%                   (DEFAULT = 'average')
%
%   'sizeBounds'    A 2-element vector specifying the min and max values in
%                   SIZEDATA to which the point sizes will be scaled.
%                   Values at or below the min value in SIZEBOUNDS will be
%                   set to the minimum point size, whereas values at or
%                   above the max value in SIZEBOUNDS will be set to the
%                   maximum point size. 
%                   (DEFAULT = min and max of SIZEDATA)
%
%   'colorBounds'   A 2-element vector specifying the min and max values in
%                   COLORDATA to which the point colors will be scaled.
%                   Values at or below the min value in COLORBOUNDS will be
%                   set to the first color in the colormap, whereas values
%                   at or above the max value in COLORBOUNDS will be set to
%                   the last color. 
%                   (DEFAULT = min and max of COLORDATA)
%
%   'colorMap'      Colormap to which COLORDATA will be mapped.
%                   (DEFAULT = custom 'whitemagma' colormap)
%
%   'gridColor'     Color of the grid.
%                   (DEFAULT = 'none')
%
%   'marker'        Marker shape.
%                   (DEFAULT = 'o')
%
%   'markerEdgeColor'   Color of marker edges (can also be 'none')
%                       (DEFAULT = 'k')
%
%   'markerSizeRange'   A 2-element vector specifying the minimum and
%                       maximum sizes between which the markers will vary.
%                       (DEFAULT = [20 200])
%
%
% Jonathan Robinson, 2020-02-06


% set defaults
opt.colnames = (1:size(colordata,2))';
opt.rownames = (1:size(colordata,1))';
opt.clusterdim  = 'both';
opt.clusterdist = 'euclidean';
opt.linkage     = 'average';
opt.sizebounds  = [min(sizedata(:)), max(sizedata(:))];
opt.colorbounds = [min(colordata(:)), max(colordata(:))];
opt.colormap    = 'whitemagma';
opt.gridcolor   = 'none';
opt.marker      = 'o';
opt.markeredgecolor = 'k';
opt.markersizerange = [20 200];

% overwrite defaults with input settings (if provided)
opt = modifyOptSettings(opt,varargin);

% verify input dimensions
if ~isequal(size(sizedata),size(colordata))
    error('"sizedata" and "colordata" must have the same dimensions.');
end
if length(opt.colnames) ~= size(colordata,2)
    error('"colnames" must have the same numer of entries as columns in "sizedata" and "colordata".');
end
if length(opt.rownames) ~= size(colordata,1)
    error('"rownames" must have the same numer of entries as rows in "sizedata" and "colordata".');
end

% determine data upon which to perform clustering 
switch opt.clusterdim
    case 'size'
        clustdata = sizedata;
    case 'color'
        clustdata = colordata;
    case 'both'
        % normalize and average size and color data
        sizedata_norm = (sizedata - min(sizedata(:)))./var(sizedata(:),'omitnan');
        sizedata_norm(isnan(sizedata_norm)) = 0;  % remove NaNs
        colordata_norm = (colordata - min(colordata(:)))./var(colordata(:),'omitnan');
        colordata_norm(isnan(colordata_norm)) = 0;  % remove NaNs
        clustdata = (sizedata_norm + colordata_norm)./2;
    otherwise
        clustdata = [];
end
clustdata(isnan(clustdata)) = 0;  % not ideal, but one way to deal with NaNs

% perform hierarchical clustering to sort rows and columns (if specified)
if ~isempty(clustdata)
    L = linkage(clustdata, opt.linkage, opt.clusterdist);
    row_ind = optimalleaforder(L, pdist(clustdata, opt.clusterdist));
    L = linkage(clustdata', opt.linkage, opt.clusterdist);
    col_ind = optimalleaforder(L, pdist(clustdata', opt.clusterdist));
else
    row_ind = 1:size(sizedata,1);
    col_ind = 1:size(sizedata,2);
end

% reorder data matrix according to clustering results
sizedata_sort = sizedata(row_ind,col_ind);
colordata_sort = colordata(row_ind,col_ind);
sortrows = opt.rownames(row_ind);
sortcols = opt.colnames(col_ind);

% check if data is square matrix with identical row and column names
if (length(opt.colnames) == length(opt.rownames)) && all(strcmp(opt.colnames,opt.rownames))
    % flip data so the diagonal is from upper left to lower right
    sizedata_sort = fliplr(sizedata_sort);
    colordata_sort = fliplr(colordata_sort);
    sortcols = flipud(sortcols);
end

% scale data for plotting
sizedata_sort(sizedata_sort < opt.sizebounds(1)) = opt.sizebounds(1);
sizedata_sort(sizedata_sort > opt.sizebounds(2)) = opt.sizebounds(2);
sizedata_sort = range(opt.markersizerange) / range(opt.sizebounds) * (sizedata_sort - opt.sizebounds(1)) + opt.markersizerange(1);

% generate point locations
[ny,nx] = size(sizedata_sort);
x = repmat(1:nx,ny,1); x = x(:);
y = repmat(1:ny,nx,1)'; y = y(:);

% generate grid data
xgrid = repmat(2:nx,2,1)' - 0.5;
ygrid = repmat(2:ny,2,1)' - 0.5;


% generate plot
figure;
a = axes;
set(a,'YAxisLocation','Right');
set(a,'TickLength',[0 0],'XLim',[0.5 nx+0.5],'YLim',[0.5 ny+0.5]);
hold on

% draw grid lines
if ~isempty(opt.gridcolor) && ~strcmpi(opt.gridcolor,'none')
    line(xgrid,[0.5,ny+0.5],'color',opt.gridcolor);
    line([0.5,nx+0.5],ygrid,'color',opt.gridcolor);
end

% draw points
if strcmpi(opt.markeredgecolor,'none')
    scatter(x,y,sizedata_sort(:),colordata_sort(:),opt.marker,'filled');
else
    scatter(x,y,sizedata_sort(:),colordata_sort(:),opt.marker,'filled','MarkerEdgeColor',opt.markeredgecolor,'LineWidth',0.5);
end
set(gca,'XTick', 1:nx);
set(gca,'YTick', 1:ny);
set(gca,'XTickLabels',sortcols,'YTickLabels',sortrows);
set(gca,'XTickLabelRotation',90);

% bad form to have try/catch, but it works for now
try
    colormap(opt.colormap);  % check if Matlab colormap exists
catch
    colormap(custom_cmap(opt.colormap));  % check if custom colormap exists
end

if ~isempty(opt.colorbounds)
    caxis(opt.colorbounds);
end

% Matlab will sometimes break up figures with polygons, which can add
% strange cuts through the figure. This command will prevent that.
set(gcf,'Renderer','painters');


