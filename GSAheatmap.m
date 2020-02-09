function [] = GSAheatmap(GSAres,varargin)
% GSAheatmap  Generate a heatmap to visualize GSA results.
%
%
% Usage:
%
%   GSA_heatmap(GSAres, ...);
%
%
% Inputs:
%
%   GSAres          GSA results structure obtained from the geneSetAnalysis
%                   function.
%
%
% Additional Settings:
%
%   'adjusted'        If TRUE, use the adjusted p-values from the GSA.
%                     If FALSE, use the non-adjusted p-values.
%                     (DEFAULT = TRUE)
%
%   'filterMethod'    Method for filtering out which gene sets to show in
%                     the heatmap.
%
%                     'none'      no filtering - all gene sets are shown
%
%                     'pval'      keep gene sets that contain at least one
%                                 p-value below filterThreshold
%
%                     'top each'  (DEFAULT) keep the X lowest-pvalue gene
%                                 sets for each directionality type, where
%                                 X = filterThreshold
%
%                     'top all'   keep the X lowest-pvalue gene sets
%                                 overall among all directionality types,
%                                 where X = filterThreshold
%
%   'filterThreshold' Filter threshold used to determine which gene sets to
%                     show in the heatmap.
%
%                     DEFAULTS:  filterMethod   filterThreshold
%                                'none'         (not used)
%                                'pval'         0.05
%                                'top each'     10
%                                'top all'      30
%
%   'dirType'         String specifying the directional class of gene set
%                     p-values to show in the heatmap.
%
%                     NOTE: If "GSAres" is a single GSA result table, the
%                           non-directional p-values will ALWAYS be shown.
%                           If "GSAres" is an array of multiple GSA result
%                           tables, the "mixed-directional" p-values will
%                           NOT be shown.
%                           
%                     'all'       (DEFAULT) show all directional classes of
%                                 p-values (non-, mixed-, and distinct-
%                                 directional).
%
%                     'nondir'    show only non-directional p-values.
%
%                     'mixed'     show mixed-directional p-values.
%
%                     'distinct'  show distinct-directional p-values.
%
%   'colorMax'        The -log10(p value) corresponding to the darkest
%                     color on the heatmap colorbar.
%                     (DEFAULT = maximum -log10 pvalue)
%
%
% Jonathan Robinson, 2020-02-08


%% Handle inputs

% set defaults
opt.adjusted = true;
opt.dirtype = 'all';
opt.colormax = [];  % will be set based on the data later on

% since the filterThreshold default depends on the filterMethod, we need to
% already check this in the inputs
ind = find(strcmpi('filtermethod',varargin));
if isempty(ind)
    opt.filtermethod = 'top each';  % default
else
    opt.filtermethod = varargin{ind+1};
end

% allow abbreviation of "filterThresh"
varargin(strcmpi('filterthresh',varargin)) = {'filterthreshold'};

% assign filterThreshold default based on filterMethod
switch lower(opt.filtermethod)
    case {'pval','p-val','pvalue','p-value'}
        opt.filterthreshold = 0.05;
    case {'top each','topeach'}
        opt.filterthreshold = 10;
    case {'top all','topall'}
        opt.filterthreshold = 30;
    case 'none'
        opt.filterthreshold = [];
    otherwise
        error('"%s" is not a recognized filterMethod. Valid options are "none", "pval", "top each", or "top all".', opt.filtermethod);
end

% overwrite defaults with input settings (if provided)
opt = modifyOptSettings(opt,varargin);

% verify dirType input
if contains(opt.dirtype,'mix')
    if iscell(GSAres)
        error('The mixed-directional dirType option is not available for multi-GSA results.');
    end
    opt.dirtype = 'mix';
elseif contains(opt.dirtype,'dist')
    opt.dirtype = 'dist';
elseif contains(opt.dirtype,'non')
    opt.dirtype = 'non';
elseif ~strcmp(opt.dirtype,'all')
    error('"%s" is not a recognized dirType. Valid options are "all", "nondir", "mixed", or "distinct".', opt.dirtype);
end


%% Extract and prepare p-value data

if istable(GSAres)  % single GSA result table
    
    % determine which p-value types to plot
    table_cols = {'p_distdn';'p_mixdn';'p_nondir';'p_mixup';'p_distup'};
    colnames = {'dist-down';'mix-down';'non-dir';'mix-up';'dist-up'};
    if ~strcmp(opt.dirtype, 'all')
        keep_cols = contains(colnames, [{'non'}, opt.dirtype]);  % always show non-directional
        table_cols = table_cols(keep_cols);
        colnames = colnames(keep_cols);
    end
    if ( opt.adjusted )
        table_cols = regexprep(table_cols,'p_','padj_');
    end
    [keep,col_inds] = ismember(table_cols,GSAres.Properties.VariableNames);
    col_inds = col_inds(keep);  % in case some p-value classes are missing
    colnames = colnames(keep);
    pData = table2array(GSAres(:,col_inds));
    
    % set max color value if not specified
    if isempty(opt.colormax)
        opt.colormax = max(-log10(pData(:)));
    end
    
    % assign rownames (gene-set names)
    rownames = GSAres.GS_name;
    
    % filter rows from p-value matrix
    keep_rows = filter_pData(pData, opt.filtermethod, opt.filterthreshold);
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
    elseif size(log_pData,2) == 3
        % if mix-directional or distinct-directional p-values are missing
        dir_score = (log_pData(:,3) - log_pData(:,1)).*max(log_pData(:,2), 0.1);
    else
        dir_score = log_pData;  % non-directional only
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
        if ~(opt.adjusted)
            pData_nondir(:,i) = GSAres{i}.p_nondir(ind);
            if ismember('p_distdn', GSAres{i}.Properties.VariableNames)
                pData_distdn(:,i) = GSAres{i}.p_distdn(ind);
                pData_distup(:,i) = GSAres{i}.p_distup(ind);
            end
        else
            pData_nondir(:,i) = GSAres{i}.padj_nondir(ind);
            if ismember('padj_distdn', GSAres{i}.Properties.VariableNames)
                pData_distdn(:,i) = GSAres{i}.padj_distdn(ind);
                pData_distup(:,i) = GSAres{i}.padj_distup(ind);
            end
        end
    end
    
    % log-transform p-values, and merge directional p-values
    log_pData_nondir = -log10(pData_nondir);
    pData_distdir = min(pData_distdn, pData_distup);
    log_pData_distdir = -log10(pData_distdir);
    
    % if p_distdn < p_distup, then switch the sign (make negative)
    log_pData_distdir(pData_distdn < pData_distup) = -log_pData_distdir(pData_distdn < pData_distup);
    
    % set max color value if not specified
    if isempty(opt.colormax)
        opt.colormax = max([log_pData_nondir(:); abs(log_pData_distdir(:))]);
    end
    
    % assign rownames (gene-set names)
    rownames = GSnames_intersect;
    
    % filter rows from non-directional p-value matrix
    if ~strcmp(opt.dirtype,'dist')
        keep_rows = filter_pData(pData_nondir, opt.filtermethod, opt.filterthreshold);
        if isempty(keep_rows)
            if strcmp(opt.dirtype,'all')
                warning('No gene sets passed the filter on non-directional p-values - plot will not be shown!');
            elseif strcmp(opt.dirtype,'non')
                error('No gene sets passed the filter on non-directional p-values - try relaxing the constraints.');
            end
        end
    else
        keep_rows = [];
    end
    log_pData_nondir = log_pData_nondir(keep_rows,:);
    rownames_nondir = rownames(keep_rows);
    
    % filter rows from distinct-directional p-value matrix
    if ~strcmp(opt.dirtype,'non')
        keep_rows = filter_pData(pData_distdir, opt.filtermethod, opt.filterthreshold);
        if isempty(keep_rows)
            if strcmp(opt.dirtype,'all')
                warning('No gene sets passed the filter on distinct-directional p-values - plot will not be shown!');
            elseif strcmp(opt.dirtype,'dist')
                error('No gene sets passed the filter on distinct-directional p-values - try relaxing the constraints.');
            end
        end
    else
        keep_rows = [];
    end
    log_pData_distdir = log_pData_distdir(keep_rows,:);
    rownames_distdir = rownames(keep_rows);
    
else
    error('GSAres is in an unrecognized format. It should be a table, or a cell array of tables.');
end


%% Generate heatmap

if istable(GSAres)
    
    % define custom colormap
    hotcmap = flipud(hot(120));
    cmap_down = custom_cmap([0.1 0.1 0.8;0 0.7 0.9;1 1 1],100,[0.25;0.4;0.35]);
    cmap_nondir = custom_cmap([0 0 0;1 1 1]);
    cmap_up = hotcmap(1:100,:);
    
    % Columns will be colored according to their directionality.
    % Distinct- and mixed-directional-down are colored according to the
    % first colormap in CMAP, non-directional is colored according to the
    % second colormap in CMAP, and distinct- and mixed-directional-up are
    % are colored according to the third color in CMAP.
    cmap = [cmap_down; cmap_nondir; cmap_up];
    if length(colnames) == 5
        log_pData(:,1:2) = min(log_pData(:,1:2), 0.999999*opt.colormax);
        log_pData(:,3)   = min(log_pData(:,3)   + 1.000001*opt.colormax, 1.999999*opt.colormax);
        log_pData(:,4:5) = min(log_pData(:,4:5) + 2.000001*opt.colormax, 2.999999*opt.colormax);
        colorBounds = [0, 3 * opt.colormax];
    elseif length(colnames) == 3
        log_pData(:,1) = min(log_pData(:,1), 0.999999*opt.colormax);
        log_pData(:,2) = min(log_pData(:,2) + 1.000001*opt.colormax, 1.999999*opt.colormax);
        log_pData(:,3) = min(log_pData(:,3) + 2.000001*opt.colormax, 2.999999*opt.colormax);
        colorBounds = [0, 3 * opt.colormax];
    elseif length(colnames) == 1
        % use a different heatmap if only non-directional will be shown,
        % otherwise it is pretty boring
        cmap = custom_cmap('whitemagma',100);
        log_pData = min(log_pData, opt.colormax);
        colorBounds = [0, opt.colormax];
    else
        error('Incorrect number of columns.');  % should never see this error
    end
    
    
    % trim very long gene set names
    maxChar = 75;
    longNames = cellfun(@(s) length(s) > maxChar, rownames);
    rownames(longNames) = cellfun(@(s) [s(1:maxChar) '...'], rownames(longNames), 'UniformOutput', false);
    
    % generate heatmap
    genHeatMap(log_pData, 'colNames', colnames, 'rowNames', rownames, ...
        'clusterDim', 'none', 'colorMap', cmap, 'colorBounds', colorBounds, ...
        'gridColor', 'k');
    
    % scale plot horizontally to deal with too long or short gene set names
    maxLength = max(cellfun(@length, rownames));
    scaleX = min(max(0.5, 1-maxLength/50), 0.8);  % scale between 0.5x - 0.8x
    pos = get(gca,'Position');
    set(gca,'Position',pos .* [1 1 scaleX 1]);
    
    % add custom colorbars
    if length(colnames) == 1  % non-directional only

        % add colorbar
        c = colorbar('Location','NorthOutside');
        xlabel(c,'log_1_0(p.non.dir)','FontSize',10);
                
    else
        
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
        c1 = colorbar('Location','NorthOutside','Position',c1Pos,'TickLabels',[],'Limits',[0 opt.colormax]);
        c2 = colorbar('Location','NorthOutside','Position',c2Pos,'TickLabels',[],'Limits',[0 opt.colormax]);
        c3 = colorbar('Location','NorthOutside','Position',c3Pos,'FontSize',10,'Limits',[0 opt.colormax]);
        
        % label colorbar axis
        xlabel(c3,'-log_1_0(p-value)','FontSize',10);
        
        % map each colorbar to a different p-value direction type
        colormap(c1,repmat(cmap_down,3,1));
        colormap(c2,repmat(cmap_nondir,3,1));
        colormap(c3,repmat(cmap_up,3,1));
        
    end
    
elseif iscell(GSAres)
    
    % retrieve or create column names
    colnames = getFieldContentsFromArray(GSAres,'Properties.Description');
    if all(cellfun(@isempty, colnames))
        % if GSAres tables are not labeled, just use numbers (1, 2, 3, ...)
        colnames = arrayfun(@num2str, (1:numel(GSAres))');
    end

    % plot non-directional p-value heatmap
    if ismember(opt.dirtype,{'non','all'}) && ~isempty(log_pData_nondir)
        % trim very long gene set names
        maxChar = 75;
        longNames = cellfun(@(s) length(s) > maxChar, rownames_nondir);
        rownames_nondir(longNames) = cellfun(@(s) [s(1:maxChar) '...'], rownames_nondir(longNames), 'UniformOutput', false);
        
        % generate heatmap
        genHeatMap(log_pData_nondir, 'colNames', colnames, 'rowNames', ...
            rownames_nondir, 'clusterDim', 'both', 'clusterDist', 'euclid', ...
            'colorMap', custom_cmap('whitemagma'), 'colorBounds', [0,opt.colormax], ...
            'gridColor','k');
        
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
    if ismember(opt.dirtype,{'dist','all'}) && ~isempty(log_pData_distdir)
        % trim very long gene set names
        maxChar = 75;
        longNames = cellfun(@(s) length(s) > maxChar, rownames_distdir);
        rownames_distdir(longNames) = cellfun(@(s) [s(1:maxChar) '...'], rownames_distdir(longNames), 'UniformOutput', false);
        
        % generate heatmap
        genHeatMap(log_pData_distdir, 'colNames', colnames, 'rowNames', ...
            rownames_distdir, 'clusterDim', 'both', 'clusterDist', 'euclid', ...
            'colorMap', custom_cmap('redblue'), 'colorBounds', [-opt.colormax,opt.colormax], ...
            'gridColor','k');
        
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


function keep_rows = filter_pData(pData,filterMethod,filterThreshold)
% Determine which rows of a p-value matrix to keep based on filter criteria

switch filterMethod
    case 'pval'
        % keep rows that contain at least one entry below threshold p-val
        keep_rows = find(any(pData < filterThreshold,2));
    case 'top each'
        % keep the X lowest-pvalue-rows of each column, where X = filterThreshold
        [~,sort_ind] = sort(pData);
        keep_rows = unique(sort_ind(1:filterThreshold,:));
    case 'top all'
        % keep the X lowest-pvalue-rows overall, where X = filterThreshold
        [~,sort_ind] = sort(min(pData,[],2));
        keep_rows = sort_ind(1:filterThreshold);
    case 'none'
        % keep all rows
        keep_rows = 1:size(pData,1);
    otherwise
        error('Invalid filtering method.');
end

end

