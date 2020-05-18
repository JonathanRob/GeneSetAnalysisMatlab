function h = genHeatMap(data,varargin)
%genHeatMap  Generate a heatmap for a given matrix of data.
%
% Usage:
%
%   genHeatMap(data, 'OPTION1', 'VALUE1', 'OPTION2', 'VALUE2', ...);
%
%   For example, to use the Hamming distance as the cluster distance metric
%   and a black grid color:
%
%   genHeatMap(data, 'clusterDist', 'hamming', 'gridColor', 'k');
%
%
% Input:
%
%   data        Numerical matrix.
%
%
% Additional Settings:
%
%   'colNames'     A list of names of data columns.
%                  (DEFAULT = column index numbers)
% 
%   'rowNames'     A list of names of data rows.
%                  (DEFAULT = row index numbers)
%
%   'clusterDim'   'none'  the data will be plotted as provided
%                  'rows'  cluster/rearrange the rows based on distance
%                  'cols'  cluster/rearrange the columns based on distance
%                  'both'  (DEFAULT) cluster/rearrange rows and columns
%                          based on distance
%
%   'clusterDist'  Distance metric to be used for clustering, ignored if
%                  'clusterDim' is 'none'. Options are the same as those
%                  in, e.g., PDIST ('euclidean', 'hamming', etc.).
%                  (DEFAULT = 'euclidean')
%
%   'linkage'      String specifying the linkage method to be used for
%                  hierarchical clustering (see built-in LINKAGE function
%                  for options).
%                  (DEFAULT = 'average')
%
%   'colorBounds'  A 2-element vector with min and max values to manually
%                  set the bounds of the colormap.
%                  (DEFAULT = min/max of data)
%
%   'colorMap'     Colormap provided as string (e.g. 'parula', 'hot', etc.)
%                  or an Nx3 RGB matrix of N colors.
%                  (DEFAULT = 'whitemagma')
%
%   'gridColor'    String or 1x3 RGB vector specifying the color of grid
%                  lines.
%                  (DEFAULT = 'none')
%
% Output:
%
%   h           (Optional) Handle to a pcolor plotting object. The heatmap
%               will be plotted automatically.
%
%
% Jonathan Robinson, 2020-05-18


% set defaults
opt.colnames = (1:size(data,2))';
opt.rownames = (1:size(data,1))';
opt.clusterdim  = 'both';
opt.clusterdist = 'euclidean';
opt.colormap    = 'whitemagma';
opt.colorbounds = [min(data(:)),max(data(:))];
opt.gridcolor   = 'none';
opt.linkage     = 'average';

% overwrite defaults with input settings (if provided)
opt = modifyOptSettings(opt,varargin);

% further verify some inputs
if ~ismember(opt.clusterdim,{'none','rows','cols','both'})
    error('"%s" is not a valid "clusterDim" option. Choose "none", "rows", "cols", or "both".',opt.clusterdim);
end

% perform hierarchical clustering to sort rows (if specified)
if ismember(opt.clusterdim,{'rows','both'})
    L = linkage(data,opt.linkage,opt.clusterdist);
    row_ind = optimalleaforder(L,pdist(data,opt.clusterdist));
else
    row_ind = 1:size(data,1);
end
% perform hierarchical clustering to sort columns (if specified)
if ismember(opt.clusterdim,{'cols','both'})
    L = linkage(data',opt.linkage,opt.clusterdist);
    col_ind = optimalleaforder(L,pdist(data',opt.clusterdist));
else
    col_ind = 1:size(data,2);
end

% reorder data matrix according to clustering results
sortdata = data(row_ind,col_ind);
sortrows = opt.rownames(row_ind);
sortcols = opt.colnames(col_ind);

% check if data is square matrix with identical row and column names
if (length(opt.colnames) == length(opt.rownames)) && all(strcmp(opt.colnames,opt.rownames))
    % flip data so the diagonal is from upper left to lower right
    sortdata = fliplr(sortdata);
    sortcols = flipud(sortcols);
end

% pad data matrix with zeros (pcolor cuts off last row and column)
sortdata(end+1,end+1) = 0;

% generate pcolor plot
figure;
a = axes;
set(a,'YAxisLocation','Right','XTick',[],'YTick', (1:size(sortdata,1))+0.5,'YTickLabels',sortrows);
set(a,'TickLength',[0 0],'XLim',[1 size(sortdata,2)],'YLim',[1 size(sortdata,1)]);
hold on

h = pcolor(sortdata);
set(h,'EdgeColor',opt.gridcolor);
set(gca,'XTick', (1:size(sortdata,2))+0.5);
set(gca,'YTick', (1:size(sortdata,1))+0.5);
set(gca,'XTickLabels',sortcols,'YTickLabels',sortrows);
set(gca,'XTickLabelRotation',90);

% bad form to have try/catch, but it works for now
try
    colormap(opt.colormap);
catch
    colormap(custom_cmap(opt.colormap));
end

if ~isempty(opt.colorbounds)
    caxis(opt.colorbounds);
end


