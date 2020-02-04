function [] = GSAheatmap(GSAres,adjusted,filterMethod,filterThresh,colorMax,dirType)
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
%   dirType         String specifying the type of gene set p-values to use
%                   in the heatmap.
%                   
%                   'nondir'  use non-directional p-values
%
%                   'dir'     use distinct-directional p-values
%
%                   'both'    (Default) use both non- and distinct-
%                             directional p-values (creates 2 figures)
%
%                   NOTE: This input is only used when GSAres is a cell
%                   array of multiple GSA result structures.
%
%
% Jonathan Robinson, 2020-01-22


%% Handle inputs

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
if nargin < 6
    dirType = 'both';
elseif ismember(dirType,{'dir','distdir','dist-dir','dist.dir','directional'})
    dirType = 'dir';
elseif ismember(dirType,{'nondir','non-dir','non.dir','nondirectional','non-directional'})
    dirType = 'nondir';
else
    error('dirType not recognized. Must be "nondir", "dir", or "both".');
end


%% Extract and prepare p-value data

if istable(GSAres)  % single GSA result table
    
    % extract p-value data from GSAres
    table_cols = {'p_distdn';'p_mixdn';'p_nondir';'p_mixup';'p_distup'};
    colnames = {'dist-down';'mix-down';'non-dir';'mix-up';'dist-up'};
    if ( adjusted )
        table_cols = regexprep(table_cols,'p_','padj_');
    end
    [keep,col_inds] = ismember(table_cols,GSAres.Properties.VariableNames);
    col_inds = col_inds(keep);  % in case some p-value types are missing
    colnames = colnames(keep);
    pData = table2array(GSAres(:,col_inds));
    
    % set max color value if not specified
    if isempty(colorMax)
        colorMax = max(-log10(pData(:)));
    end
    
    % assign rownames (gene-set names)
    rownames = GSAres.GS_name;
    
    % filter rows from p-value matrix
    keep_rows = filter_pData(pData,filterMethod,filterThresh);
    if isempty(keep_rows)
        error('No gene sets passed the filter - try relaxing the constraints.');
    end
    pData = pData(keep_rows,:);
    rownames = rownames(keep_rows);
    
    % log-transform p-values
    log_pData = -log10(pData);
    
    % calculate "directionality score" of each row
    % dir.score = (p.mix.up*(p.non.dir + p.dist.up) - p.mix.dn*(p.non.dir + p.dist.dn)) / (2*max(p)^2)
    if size(log_pData,2) == 5
        dir_score = (log_pData(:,4).*(log_pData(:,3) + log_pData(:,5)) - log_pData(:,2).*(log_pData(:,3) + log_pData(:,1)))./(2*max(log_pData(:)).^2);
    elseif all(strcmp(colnames,{'dist-down';'non-dir';'dist-up'})) || ...
            all(strcmp(colnames,{'mix-down';'non-dir';'mix-up'}))
        % if mix-directional or distinct-directional p-values are missing
        dir_score = (log_pData(:,3) - log_pData(:,1)).*log_pData(:,2);
    end
    [~,sort_ind] = sort(dir_score);
    
    % sort rows of pData by directionality score
    log_pData = log_pData(sort_ind,:);
    rownames = rownames(sort_ind);
    
elseif iscell(GSAres)  % array of multiple GSA result tables
    
    % find gene sets that are not present in all GSAres tables
    GSnames = getFieldContentsFromArray(GSAres,'GS_name');
    GSnames_intersect = multiIntersect(GSnames{:});
    GSnames_union = unique(vertcat(GSnames{:}));
    GSremoved = setdiff(GSnames_union, GSnames_intersect);
    if ~isempty(GSremoved)
        fprintf('\nWarning! The following gene sets are not present in all GSAres tables and will therefore be ignored:\n');
        fprintf('\t%s\n',GSremoved{1:min(5,numel(GSremoved))});  % only print first 5
        if numel(GSremoved) > 5
            fprintf('\t...and %u more.\n\n',numel(GSremoved)-5);
        end
    end
    
    [pData_nondir,pData_distdn,pData_distup] = deal(ones(numel(GSnames_intersect),numel(GSAres)));
    for i = 1:numel(GSAres)
        [~,ind] = ismember(GSnames_intersect, GSAres{i}.GS_name);
        if ~(adjusted)
            pData_nondir(:,i) = GSAres{i}.p_nondir(ind);
            pData_distdn(:,i) = GSAres{i}.p_distdn(ind);
            pData_distup(:,i) = GSAres{i}.p_distup(ind);
        else
            pData_nondir(:,i) = GSAres{i}.padj_nondir(ind);
            pData_distdn(:,i) = GSAres{i}.padj_distdn(ind);
            pData_distup(:,i) = GSAres{i}.padj_distup(ind);
        end
    end
    
    % log-transform p-values, and merge directional p-values
    log_pData_nondir = -log10(pData_nondir);
    pData_distdir = min(pData_distdn, pData_distup);
    log_pData_distdir = -log10(pData_distdir);
    
    % if p_distdn < p_distup, then switch the sign (make negative)
    log_pData_distdir(pData_distdn < pData_distup) = -log_pData_distdir(pData_distdn < pData_distup);
    
    % set max color value if not specified
    if isempty(colorMax)
        colorMax = max([log_pData_nondir(:); abs(log_pData_distdir(:))]);
    end
    
    % assign rownames (gene-set names)
    rownames = GSnames_intersect;
    
    % filter rows from non-directional p-value matrix
    if ~strcmp(dirType,'dir')
        keep_rows = filter_pData(pData_nondir,filterMethod,filterThresh);
        if isempty(keep_rows)
            if strcmp(dirType,'both')
                warning('No gene sets passed the filter on non-directional p-values - plot will not be shown!');
            elseif strcmp(dirType,'nondir')
                error('No gene sets passed the filter on non-directional p-values - try relaxing the constraints.');
            end
        end
    else
        keep_rows = [];
    end
    log_pData_nondir = log_pData_nondir(keep_rows,:);
    rownames_nondir = rownames(keep_rows);
    
    % filter rows from distinct-directional p-value matrix
    if ~strcmp(dirType,'nondir')
        keep_rows = filter_pData(pData_distdir,filterMethod,filterThresh);
        if isempty(keep_rows)
            if strcmp(dirType,'both')
                warning('No gene sets passed the filter on distinct-directional p-values - plot will not be shown!');
            elseif strcmp(dirType,'nondir')
                error('No gene sets passed the filter on distinct-directional p-values - try relaxing the constraints.');
            end
        end
    else
        keep_rows = [];
    end
    log_pData_distdir = log_pData_distdir(keep_rows,:);
    rownames_distdir = rownames(keep_rows);
    
else
    error('GSAres is in an unrecognized format. Should be a table, or a cell array of tables.');
end


%% Generate heatmap

if istable(GSAres)
    
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
    if length(colnames) == 5
        log_pData(:,1:2) = min(log_pData(:,1:2),0.999999*colorMax);
        log_pData(:,3)   = min(log_pData(:,3) + 1.000001*colorMax,1.999999*colorMax);
        log_pData(:,4:5) = min(log_pData(:,4:5) + 2.000001*colorMax,2.999999*colorMax);
    elseif length(colnames) == 3
        log_pData(:,1) = min(log_pData(:,1),0.999999*colorMax);
        log_pData(:,2) = min(log_pData(:,2) + 1.000001*colorMax,1.999999*colorMax);
        log_pData(:,3) = min(log_pData(:,3) + 2.000001*colorMax,2.999999*colorMax);
    else
        error('Incorrect number of columns.');
    end
    colorBounds = [0,3*colorMax];
    
    % trim very long gene set names
    maxChar = 75;
    longNames = cellfun(@(s) length(s) > maxChar, rownames);
    rownames(longNames) = cellfun(@(s) [s(1:maxChar) '...'], rownames(longNames), 'UniformOutput', false);
    
    % generate heatmap
    genHeatMap(log_pData, colnames, rownames, 'none', [], cmap, colorBounds, 'k');
    
    % scale plot horizontally to deal with too long or short gene set names
    maxLength = max(cellfun(@length, rownames));
    scaleX = min(max(0.5, 1-maxLength/50), 0.8);  % scale between 0.5x - 0.8x
    pos = get(gca,'Position');
    set(gca,'Position',pos .* [1 1 scaleX 1]);
    
    % add custom colorbars
    % reduce figure height to make room for the colorbars
    scaleY = min(0.85, 0.75 + size(log_pData,1)*0.002);
    plotPos = get(gca,'Position') .* [1 1 1 scaleY];
    set(gca,'Position', plotPos);
    
    % determine the position of the colorbars based on the plot position
    cbar_width = max(0.02, 0.04 - size(log_pData,1)*0.001);  % scale width based on number of rows
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
    
elseif iscell(GSAres)
    
    % retrieve or create column names
    colnames = getFieldContentsFromArray(GSAres,'Properties.Description');
    if all(cellfun(@isempty, colnames))
        % if GSAres tables are not labeled, just use numbers (1, 2, 3, ...)
        colnames = arrayfun(@num2str, (1:numel(GSAres))');
    end

    % plot non-directional p-value heatmap
    if ismember(dirType,{'nondir','both'}) && ~isempty(log_pData_nondir)
        % trim very long gene set names
        maxChar = 75;
        longNames = cellfun(@(s) length(s) > maxChar, rownames_nondir);
        rownames_nondir(longNames) = cellfun(@(s) [s(1:maxChar) '...'], rownames_nondir(longNames), 'UniformOutput', false);
        
        % generate heatmap
        genHeatMap(log_pData_nondir, colnames, rownames_nondir, 'both', ...
            'euclid', custom_cmap('magma'), [0,colorMax], 'k');
        
        % scale plot horizontally to deal with too long or short gene set names
        maxLength = max(cellfun(@length, rownames_nondir));
        scaleX = min(max(0.5, 1-maxLength/50), 0.8);  % scale between 0.5x - 0.8x
        pos = get(gca,'Position');
        set(gca,'Position',pos .* [1 1 scaleX 1]);
        
        % add colorbar
        c = colorbar('Location','NorthOutside');
        xlabel(c,'-log_1_0(p.non.dir)','FontSize',10);
    end
    
     % plot distinct-directional p-value heatmap
    if ismember(dirType,{'dir','both'}) && ~isempty(log_pData_distdir)
        % trim very long gene set names
        maxChar = 75;
        longNames = cellfun(@(s) length(s) > maxChar, rownames_distdir);
        rownames_distdir(longNames) = cellfun(@(s) [s(1:maxChar) '...'], rownames_distdir(longNames), 'UniformOutput', false);
        
        % generate heatmap
        genHeatMap(log_pData_distdir, colnames, rownames_distdir, 'both', ...
            'euclid', custom_cmap('redblue'), [-colorMax,colorMax], 'k');
        
        % scale plot horizontally to deal with too long or short gene set names
        maxLength = max(cellfun(@length, rownames_distdir));
        scaleX = min(max(0.5, 1-maxLength/50), 0.8);  % scale between 0.5x - 0.8x
        pos = get(gca,'Position');
        set(gca,'Position',pos .* [1 1 scaleX 1]);
        
        % add colorbar
        c = colorbar('Location','NorthOutside');
        xlabel(c,'Signed log_1_0(p.dist.dir)','FontSize',10);
    end
end


end


%% Additional functions

function contents = getFieldContentsFromArray(array,field)
% Retreive the contents of a given field from each structure in an array.
% This function is necessary because dot indexing is not supported for cell
% arrays of structures.
if contains(field,'.')
    % only works for one level of subfield (e.g., "field.subfield")
    field = strsplit(field,'.');
    contents = cellfun(@(x) {x.(field{1}).(field{2})}, array);
else
    contents = cellfun(@(x) {x.(field)}, array);
end
end


function intersection = multiIntersect(varargin)
% Determine intersection of two or more input vectors, arrays, etc.
intersection = varargin{1};
for i = 2:numel(varargin)
    intersection = intersect(intersection, varargin{i});
end
end


function keep_rows = filter_pData(pData,filterMethod,filterThresh)
% Determine which rows of a p-value matrix to keep based on filter criteria

switch filterMethod
    case 'pval'
        % keep rows that contain at least one entry below threshold p-val
        keep_rows = find(any(pData < filterThresh,2));
    case 'top each'
        % keep the X lowest-pvalue-rows of each column, where X = filterThresh
        [~,sort_ind] = sort(pData);
        keep_rows = unique(sort_ind(1:filterThresh,:));
    case 'top all'
        % keep the X lowest-pvalue-rows overall, where X = filterThresh
        [~,sort_ind] = sort(min(pData,[],2));
        keep_rows = sort_ind(1:filterThresh);
    case 'none'
        % keep all rows
        keep_rows = 1:size(pData,1);
    otherwise
        error('Invalid filtering method.');
end

end

