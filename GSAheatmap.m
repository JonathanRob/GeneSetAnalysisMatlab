function [] = GSAheatmap(GSAres,adjusted,filterMethod,filterThresh,colorMax)
% GSAheatmap  Generate a heatmap to visualize GSA results.
%
%
% Usage:
%
%   GSA_heatmap(GSAres,adjusted,filterMethod,filterThresh,colorMax);
%
%
% Inputs:
%
%   GSAres          GSA results structure obtained from the geneSetAnalysis
%                   function.
%
%   adjusted        If TRUE, use the adjusted p-values from the GSA.
%                   If FALSE, use the non-adjusted p-values.
%                   (Default = TRUE)
%
%   filterMethod    Method for filtering out which gene sets to show in the
%                   heatmap.
%
%                   'none'      no filtering - all gene sets will be shown
%
%                   'pval'      keep gene sets that contain at least one
%                               directionality p-value below filterThresh
%
%                   'top each'  keep the X lowest-pvalue gene sets for each
%                               directionality type, where X = filterThresh
%                               (Default)
%
%                   'top all'   keep the X lowest-pvalue gene sets overall
%                               among all directionality types, where 
%                               X = filterThresh
%
%   filterThresh    Filter threshold to determine which gene sets to show
%                   in the heatmap.
%
%              Defaults:  filterMethod   filterThresh
%                         'none'         (not used)
%                         'pval'         0.05
%                         'top each'     10
%                         'top all'      30
%                   
%   colorMax        The -log10(p value) corresponding to the darkest color
%                   on the heatmap colorbar.
%                   (Default = maximum -log10 pvalue)
%
%
% Jonathan Robinson, 2020-01-19


%% Initialize and organize data

% assign default inputs
if nargin < 2 || isempty(adjusted)
    adjusted = true;
end
if nargin < 3 || isempty(filterMethod)
    filterMethod = 'top each';
end
if nargin < 4 || isempty(filterThresh)
    switch lower(filterMethod)
        case {'pval','p-val','pvalue','p-value'}
            filterThresh = 0.05;
        case 'top each'
            filterThresh = 10;
        case 'top all'
            filterThresh = 30;
        otherwise
            filterThresh = [];
    end
end
if nargin < 5
    colorMax = [];
end

% extract p-value data from GSAres
pValCols = {'p_distdn';'p_mixdn';'p_nondir';'p_mixup';'p_distup'};
pValColNames = {'dist-down';'mix-down';'non-dir';'mix-up';'dist-up'};
if ( adjusted )
    pValCols = regexprep(pValCols,'p_','padj_');
end
[keep,pValColInds] = ismember(pValCols,GSAres(1,:));
pValColInds = pValColInds(keep);  % in case some p-value types are missing
pValColNames = pValColNames(keep);
pData = cell2mat(GSAres(2:end,pValColInds));

% set max color value if not specified
if isempty(colorMax)
    colorMax = max(-log10(pData(:)));
end

% get rownames (gene-set names)
pValRowNames = GSAres(2:end,1);


% remove rows from data according to filtering method
switch filterMethod
    case 'pval'
        % keep rows that contain at least one entry below threshold p-val
        keep_ind = find(any(pData < filterThresh,2));
    case 'top each'
        % keep the X lowest-pvalue-rows of each column, where X = filterThresh
        [~,sort_ind] = sort(pData);
        keep_ind = unique(sort_ind(1:filterThresh,:));
    case 'top all'
        % keep the X lowest-pvalue-rows overall, where X = filterThresh
        [~,sort_ind] = sort(min(pData,[],2));
        keep_ind = sort_ind(1:filterThresh);
    case 'none'
        % keep all rows
        keep_ind = 1:size(pData,1);
    otherwise
        error('Invalid filtering method.');
end
if isempty(keep_ind)
    error('No gene sets passed the filter - try relaxing the constraints.');
end
pData = pData(keep_ind,:);
pValRowNames = pValRowNames(keep_ind);


% log-transform p-values
pData = -log10(pData);

% calculate "directionality score" of each row
% dir.score = (p.mix.up*(p.non.dir + p.dist.up) - p.mix.dn*(p.non.dir + p.dist.dn)) / (2*max(p)^2)
if size(pData,2) == 5
    dir_score = (pData(:,4).*(pData(:,3) + pData(:,5)) - pData(:,2).*(pData(:,3) + pData(:,1)))./(2*max(pData(:)).^2);
elseif all(strcmp(pValColNames,{'dist-down';'non-dir';'dist-up'})) || ...
        all(strcmp(pValColNames,{'mix-down';'non-dir';'mix-up'}))
    % if mix-directional or distinct-directional p-values are missing
    dir_score = (pData(:,3) - pData(:,1)).*pData(:,2);
end
[~,sort_ind] = sort(dir_score);

% sort rows of pData by directionality score
pData = pData(sort_ind,:);
pValRowNames = pValRowNames(sort_ind);


%% Generate heatmap

% define custom colormap
hotcmap = flipud(hot(120));
cmap_down = custom_cmap([0.1 0.1 0.8;0 0.7 0.9;1 1 1],100,[0.25;0.4;0.35]);
cmap_nondir = custom_cmap([0 0 0;1 1 1]);
cmap_up = hotcmap(1:100,:);

% Columns will be colored according to their directionality.
% Columns 1 and 2 (p dist.dir.dn and mix.dir.dn) are colored according
% to the first colormap in CMAP, column 3 (p non.dir) is colored
% according to the second colormap in CMAP, and columns 4 and 5
% (p mix.dir.up and dist.dir.up) are colored according to the third
% color in CMAP.
cmap = [cmap_down; cmap_nondir; cmap_up];
if length(pValColNames) == 5
    pData(:,1:2) = min(pData(:,1:2),0.999999*colorMax);
    pData(:,3)   = min(pData(:,3) + 1.000001*colorMax,1.999999*colorMax);
    pData(:,4:5) = min(pData(:,4:5) + 2.000001*colorMax,2.999999*colorMax);
elseif length(pValColNames) == 3
    pData(:,1) = min(pData(:,1),0.999999*colorMax);
    pData(:,2) = min(pData(:,2) + 1.000001*colorMax,1.999999*colorMax);
    pData(:,3) = min(pData(:,3) + 2.000001*colorMax,2.999999*colorMax);
else
    error('Incorrect number of columns.');
end

pData(end,end) = 3*colorMax;  % necessary for proper color mapping
colorBounds = [0,3*colorMax];


%% Generate heatmap

% use custom function to generate heatmap
genHeatMap(pData, pValColNames, pValRowNames, 'none', [], cmap, colorBounds, 'k');

% scale plot horizontally to deal with too long or short row names
maxLength = max(cellfun(@length, pValRowNames));
scaleX = min(max(0.5, 1-maxLength/100), 0.8);  % scale between 0.5x - 0.8x
pos = get(gca,'Position');
set(gca,'Position',pos .* [1 1 scaleX 1]);


%% add custom colorbars

% reduce figure height to make room for the colorbars
scaleY = min(0.85, 0.75 + size(pData,1)*0.002);
plotPos = get(gca,'Position') .* [1 1 1 scaleY];
set(gca,'Position', plotPos);

% determine the position of the colorbars based on the plot position
cbar_width = max(0.02, 0.04 - size(pData,1)*0.001);  % scale width based on number of rows
c1Pos = [plotPos(1), plotPos(2)+plotPos(4)+0.02, plotPos(3), cbar_width];
c2Pos = [plotPos(1), plotPos(2)+plotPos(4)+0.02+1.5*cbar_width, plotPos(3), cbar_width];
c3Pos = [plotPos(1), plotPos(2)+plotPos(4)+0.02+3.0*cbar_width, plotPos(3), cbar_width];

% add colorbars on top of plot
c1 = colorbar('Location','NorthOutside','Position',c1Pos,'TickLabels',[],'Limits',[0 colorMax]);
c2 = colorbar('Location','NorthOutside','Position',c2Pos,'TickLabels',[],'Limits',[0 colorMax]);
c3 = colorbar('Location','NorthOutside','Position',c3Pos,'FontSize',10,'Limits',[0 colorMax]);

% label colorbar axis
xlabel(c3,'-log_1_0(p-value)','FontSize',10);

% map each colorbar to a different p-value direction type
colormap(c1,repmat(cmap_down,3,1));
colormap(c2,repmat(cmap_nondir,3,1));
colormap(c3,repmat(cmap_up,3,1));



