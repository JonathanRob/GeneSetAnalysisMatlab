function [GSAres,gene_sets_proc] = geneSetAnalysis( ...
    genes,pvals,dirs,gene_sets,method,nperms,GS_size_bounds,stat_type)
%geneSetAnalysis  Perform a gene set analysis (GSA).
%
% Performs a GSA given gene-level statistics (pvals) and directionality
% (dirs) for each gene, using gene sets defined by gene_sets and the
% specified method of calculating the test-statistic. GeneSetAnalysis
% returns the calculated gene set sizes and p-values in GSAres, as well as
% the processed list of gene sets (gene_sets_proc) if requested.
%
%
% Usage:
%
%   [GSAres,gene_sets_proc] = geneSetAnalysis( ...
%       genes,pvals,dirs,gene_sets,method,nperms,GS_size_bounds,stat_type);
%
%
% Input:
%
%   genes       Column of gene names or identifiers corresponding to data
%               in PVALS and DIRS.
%
%   pvals       Column vector of gene-level statistics (e.g., p-values)
%               corresponding to GENES.
%
%   dirs        Column vector of directions associated with each gene.
%               (up = 1, down = -1, no change = 0)
%               Continuous data (such as fold-change values) will be
%               converted to ternary (+1,-1,0) values.
%               
%               *** NOTE: leave empty [] if PVALS are signed. ***
%
%   gene_sets   Nx2 cell array of gene set information, where the first
%               column contains the gene set names, and the second column
%               contains the gene names associated with each gene set.
%
%   method      Satistical method used to combine gene-level statistics:
%                   'fisher'    Fisher's method.
%                 'reporter'    Reporter method.
%                 'stouffer'    Stouffer's method.
%                 'wilcoxon'    Wilcoxon rank-sum test.
%                     'mean'    Mean (arithmetic) value.
%                  'geomean'    Mean (geometric) value.
%                   'median'    Median value.
%                     'GSEA'    Gene Set Enrichment Analysis (not complete)
%
%   nperms      Number of permutations to perform when calculating the 
%               significance of each gene set. (DEFAULT = 10,000)
% 
%   GS_size_bounds  The min and max size allowed for a gene set to be
%                   included in the analysis. (DEFAULT = [5,Inf]) 
%
%   stat_type   The type of gene-level statistic provided in PVALS:
%
%                   'p'     (DEFAULT) p-values ranging between 0 and 1,
%                           where a lower value indicates a greater score
%                           or "significance".
%
%               'other'     Some other non-p-like gene-level statistic.
%                           Numbers with a greater value will be considered
%                           to indicate a greater score or "significance".
%
%
% Output:
%
%   GSAres      A cell array containing the GSA results, including gene set
%               names and sizes, and their associated p-values (raw and
%               adjusted) for each of the relevant directionality classes.
%
%
% Jonathan Robinson, 2019-04-14


%% Handle input arguments

if strcmpi(method,'GSEA')
    % still need to code the significance calculation portion of GSEA
    error('The GSEA method is not yet available.');
end

if nargin < 8
    stat_type = 'p';
end
if nargin < 7 || isempty(GS_size_bounds)
    GS_size_bounds = [5,Inf];
end
if nargin < 6 || isempty(nperms)
    nperms = 10000;
end

if ~ismember(stat_type,{'p','other'})
    error('Invalid STAT_TYPE. Valid options are "p" or "other".');
end

if length(GS_size_bounds) == 1
    % if only one value is provided for GS_size_bounds, assume it is the minimum set size allowed
    GS_size_bounds(2) = Inf;
end
[minGSsize,maxGSsize] = deal(GS_size_bounds(1),GS_size_bounds(2));

if ismember(lower(method),{'fisher','reporter','stouffer'}) && ~strcmpi(stat_type,'p')
    error('Fisher, Reporter, and Stouffer methods are only valid for p-like statistics (STAT_TYPE = "P").');
end


% convert inputs to columns in case they are supplied as row vectors
genes = genes(:);
pvals = pvals(:);
dirs = dirs(:);


if strcmpi(stat_type,'p')
    % handle p-values equal to 0 or 1
    ind = pvals == 0;
    if ~isempty(ind)
        fprintf('*** WARNING: p-values equal to ZERO will be set to 0.5x the lowest nonzero p-value ***\n');
        pvals(ind) = 0.5*min(pvals(~ind));
    end
    ind = pvals == 1;
    if ~isempty(ind)
        fprintf('*** WARNING: p-values equal to ONE will be set to 1 - eps (2^-52) ***\n');
        pvals(ind) = 1-eps;
    end
end


if ~isempty(dirs)
    % convert DIRS to +1, -1, 0
    dirs = sign(dirs);
    fprintf('\nNOTE: genes with directional value of zero WILL be included in the analysis!\n');
    
    % alternative option: remove entries with DIRS equal to zero
    % ind = dirs == 0;
    % genes(ind) = []; pvals(ind) = []; dirs(ind) = [];
    % fprintf('\nRemoved %u genes with directional value of zero.\n',sum(ind));
end


%% Process gene sets

% remove rows of GENE_SETS that do not contain genes in provided list
fprintf('Checking for empty gene sets... ');
a = length(unique(gene_sets(:,1)));  % check number of sets before removal
ind = ~ismember(gene_sets(:,2),genes);
gene_sets(ind,:) = [];
b = length(unique(gene_sets(:,1)));  % check number of sets after removal
fprintf('Removed %u empty sets.\n',a-b);

% remove repeated rows in GENE_SETS
fprintf('Checking for duplicated rows in GENE_SETS... ');
a = size(gene_sets,1);  % check number of rows before removal
[~,gs_ind] = ismember(gene_sets,gene_sets);
[~,uniq_ind] = unique(gs_ind,'rows');
gene_sets = gene_sets(uniq_ind,:);
b = size(gene_sets,1);  % check number of rows after removal
fprintf('Removed %u duplicated rows.\n',a-b);

% determine gene set sizes and remove those not satisfying constraints
fprintf('Checking gene set sizes... ');
[GSnames,GSsizes] = cellfreq(gene_sets(:,1));
ind = (GSsizes < minGSsize) | (GSsizes > maxGSsize);
gene_sets(ismember(gene_sets(:,1),GSnames(ind)),:) = [];
% GSnames(ind) = [];
% GSsizes(ind) = [];
fprintf('Removed %u gene sets not satisfying size limits.\n',sum(ind));
fprintf('Final number of gene sets remaining: %u\n',length(unique(gene_sets(:,1))));

% handle duplicated gene names
[uGenes,freq_uGenes] = cellfreq(genes); 
dup_genes = uGenes(freq_uGenes > 1);  % find duplicated gene names
if ~isempty(dup_genes)
    fprintf('Handling duplicated gene names... ');
    for i = 1:length(dup_genes)
        
        % append repeated gene names with '_1', '_2', etc.
        gene_ind = ismember(genes,dup_genes(i));
        genes(gene_ind) = strcat(repmat(dup_genes(i),sum(gene_ind),1),'_',arrayfun(@num2str,[1:sum(gene_ind)]','UniformOutput',false));
        
        % add new GS-gene associations for each of the newly labeled genes
        GS_ind = ismember(gene_sets(:,2),dup_genes(i));
        GSrep = repmat(gene_sets(GS_ind,1),1,sum(gene_ind))';
        gene_sets = [gene_sets;[GSrep(:),repmat(genes(gene_ind),sum(GS_ind),1)]];
        
        % remove the old GS-gene associations
        gene_sets(GS_ind,:) = [];
        
    end
    fprintf('Done.\n');
end

% convert GENES and GENE_SETS to numeric arrays to speed up calculations
gene_sets_proc = gene_sets;  % save processed GENE_SETS
GSnames = unique(gene_sets(:,1),'stable');  % save list of gene set names

[~,gsnums] = ismember(gene_sets(:,1),gene_sets(:,1));  % index gene sets
[~,gnums] = ismember(gene_sets(:,2),genes);  % index genes
gene_sets = [gsnums,gnums];  % numeric gene_set

[~,genes] = ismember(genes,genes);  % numeric genes

[GSnums,GSsizes] = cellfreq(gene_sets(:,1));  % recalc gene set sizes


%% Calculate test-statistics

% generate subsets of data for mixed-directional calculations
if ~isempty(dirs)
    fprintf('Subsetting data for mixed-directional calculations... ');
    ind_up = dirs > 0;
    genes_mixup = genes(ind_up);
    pvals_mixup = pvals(ind_up);
    
    ind_dn = dirs < 0;
    genes_mixdn = genes(ind_dn);
    pvals_mixdn = pvals(ind_dn);
    
    ind_up = ismember(gene_sets(:,2),genes_mixup);
    gene_sets_mixup = gene_sets(ind_up,:);
    [GSnums_mixup,GSsizes_mixup] = cellfreq(gene_sets_mixup(:,1));
    
    ind_dn = ismember(gene_sets(:,2),genes_mixdn);
    gene_sets_mixdn = gene_sets(ind_dn,:);
    [GSnums_mixdn,GSsizes_mixdn] = cellfreq(gene_sets_mixdn(:,1));
    fprintf('Done.\n');
end

% initialize directional p-values
[pvals_distup,pvals_distdn] = deal(pvals);

fprintf('Calculating test statistic... ');
switch lower(method)
    
    case 'fisher'
        % Fisher method test statistic: -2*sum(log(p))
        statvals_nondir = arrayfun(@(nums) fisher(pvals(ismember(genes,gene_sets(ismember(gene_sets(:,1),nums),2)))),GSnums);
        
        if ~isempty(dirs)
            statvals_mixup = arrayfun(@(nums) fisher(pvals_mixup(ismember(genes_mixup,gene_sets_mixup(ismember(gene_sets_mixup(:,1),nums),2)))),GSnums_mixup);
            statvals_mixdn = arrayfun(@(nums) fisher(pvals_mixdn(ismember(genes_mixdn,gene_sets_mixdn(ismember(gene_sets_mixdn(:,1),nums),2)))),GSnums_mixdn);
        end
            
    case {'stouffer','reporter'}
        % Stouffer method test statistic: sum(norminv(1-p))/sqrt(length(p))
        % Reporter method uses the same test statistic, but significance is
        % calculated differently later on.
        statvals_nondir = arrayfun(@(nums) stouffer(pvals(ismember(genes,gene_sets(ismember(gene_sets(:,1),nums),2)))),GSnums);
        
        if ~isempty(dirs)
            statvals_mixup = arrayfun(@(nums) stouffer(pvals_mixup(ismember(genes_mixup,gene_sets_mixup(ismember(gene_sets_mixup(:,1),nums),2)))),GSnums_mixup);
            statvals_mixdn = arrayfun(@(nums) stouffer(pvals_mixdn(ismember(genes_mixdn,gene_sets_mixdn(ismember(gene_sets_mixdn(:,1),nums),2)))),GSnums_mixdn);
            
            statvals_distup = arrayfun(@(nums) stouffer_dir(pvals(ismember(genes,gene_sets(ismember(gene_sets(:,1),nums),2))),dirs(ismember(genes,gene_sets(ismember(gene_sets(:,1),nums),2)))),GSnums);
            statvals_distdn = -statvals_distup;
        end
    
    case 'wilcoxon'
        % Wilcoxon rank-sum test statistic: sum(tiedrank(p))
        
        if ~isempty(dirs)
            
            if strcmpi(stat_type,'p')
                % convert p-values into rank numbers
                
                %  *Distinct directional ranks are generated in the complicated
                %   manner below to avoid rounding errors caused by subtracting
                %   very small p-values (i.e., less than 2^-52) from one.
                pvals_distup(dirs >= 0) = tiedrank(-pvals(dirs >= 0)) + sum(dirs < 0);
                pvals_distup(dirs < 0) = tiedrank(pvals(dirs < 0));
                pvals_distdn(dirs < 0) = tiedrank(-pvals(dirs < 0)) + sum(dirs >= 0);
                pvals_distdn(dirs >= 0) = tiedrank(pvals(dirs >= 0));
                
                % rank mixed-directional p-values separately
                pvals_mixup = tiedrank(-pvals_mixup);
                pvals_mixdn = tiedrank(-pvals_mixdn);
                
            elseif strcmpi(stat_type,'other')
                
                % convert gene-level statistics to ranks
                pvals_distup = tiedrank(pvals.*dirs);
                pvals_distdn = tiedrank(-pvals.*dirs);
            
                pvals_mixup = tiedrank(pvals_mixup);
                pvals_mixdn = tiedrank(pvals_mixdn);
                
            end
        
            % calculate directional test statistics
            statvals_distup = arrayfun(@(nums) sum(pvals_distup(ismember(genes,gene_sets(ismember(gene_sets(:,1),nums),2))),1),GSnums);
            statvals_distdn = arrayfun(@(nums) sum(pvals_distdn(ismember(genes,gene_sets(ismember(gene_sets(:,1),nums),2))),1),GSnums);
            
            statvals_mixup = arrayfun(@(nums) sum(pvals_mixup(ismember(genes_mixup,gene_sets_mixup(ismember(gene_sets_mixup(:,1),nums),2))),1),GSnums_mixup);
            statvals_mixdn = arrayfun(@(nums) sum(pvals_mixdn(ismember(genes_mixdn,gene_sets_mixdn(ismember(gene_sets_mixdn(:,1),nums),2))),1),GSnums_mixdn);
            
        end
        
        % rank non-directional test statistics
        if strcmpi(stat_type,'p')
            pvals = tiedrank(-pvals);
        elseif strcmpi(stat_type,'other')
            pvals = tiedrank(pvals);
        end
        
        % calculate non-directional test statistics
        statvals_nondir = arrayfun(@(nums) sum(pvals(ismember(genes,gene_sets(ismember(gene_sets(:,1),nums),2))),1),GSnums);
        
        
%     case 'gsea'
%         % Gene Set Enrichment Analysis (GSEA)
%         
%         if ~(GSEA_override)
%             % convert PVALS to t-like statistic
%             pvals_distup = -log(pvals).*dirs;
%         end
%         [pvals_distup,sort_ind] = sort(-pvals_distup);  % sort by ascending t-stat value
%         pvals_distup = abs(pvals_distup);  % use absolute value of t-stat
%         
%         % calculate test statistics
%         statvals_distup = arrayfun(@(nums) gsea(pvals_distup,ismember(genes(sort_ind),gene_sets(ismember(gene_sets(:,1),nums),2))),GSnums);
%         statvals_distdn = -statvals_distup;
%         
%         % NOTE: non-directional and mixed-directional classes are not available for GSEA
%         statvals_nondir = zeros(size(GSnums));  % zeros as placeholder
%         statvals_mixup = zeros(size(GSnums_mixup));  % zeros as placeholder
%         statvals_mixdn = zeros(size(GSnums_mixdn));  % zeros as placeholder
        

    otherwise
        error('Invalid METHOD requested.');
end
fprintf('Done.\n');


%% Calculate significance


if isnumeric(nperms)
    
    % randomly shuffle p-values and assemble into matrix
    fprintf('Calculating significance via gene shuffling... ');
    rng('shuffle');  % re-seed random number generator
    
    maxVal = max(GSsizes);
    indmat = arrayfun(@(x,mv) randperm(x,mv),ones(1,nperms)*length(genes),repmat(maxVal,1,nperms),'UniformOutput',false);
    indmat = reshape([indmat{:}],maxVal,nperms);
    
    if strcmpi(method,'gsea')
        p_permute = indmat;  % only need indices for GSEA approach
    else
        p_permute = pvals(indmat);
    end
    
    % get unique list of different sizes of gene sets
    uniq_sizes = unique(GSsizes);
    
    % shuffle directional p-values
    if ~isempty(dirs)
        p_permute_distup = pvals_distup(indmat);
        p_permute_distdn = pvals_distdn(indmat);
        dir_permute = dirs(indmat);
        
        maxVal_mixup = max(GSsizes_mixup);
        indmat_mixup = arrayfun(@(x,mv) randperm(x,mv),ones(1,nperms)*length(genes_mixup),repmat(maxVal_mixup,1,nperms),'UniformOutput',false);
        indmat_mixup = reshape([indmat_mixup{:}],maxVal_mixup,nperms);
        p_permute_mixup = pvals_mixup(indmat_mixup);
        
        maxVal_mixdn = max(GSsizes_mixdn);
        indmat_mixdn = arrayfun(@(x,mv) randperm(x,mv),ones(1,nperms)*length(genes_mixdn),repmat(maxVal_mixdn,1,nperms),'UniformOutput',false);
        indmat_mixdn = reshape([indmat_mixdn{:}],maxVal_mixdn,nperms);
        p_permute_mixdn = pvals_mixdn(indmat_mixdn);
        
        uniq_sizes_mixup = unique(GSsizes_mixup);
        uniq_sizes_mixdn = unique(GSsizes_mixdn);
    end
    
    
    % determine background distribution (test statistics for shuffled PVALS)
    switch lower(method)
        case 'fisher'
            
            bg_stat_nondir = arrayfun(@(s) fisher(p_permute(1:s,:)),uniq_sizes,'UniformOutput',false);
            
            if ~isempty(dirs)
                bg_stat_mixup = arrayfun(@(s) fisher(p_permute_mixup(1:s,:)),uniq_sizes_mixup,'UniformOutput',false);
                bg_stat_mixdn = arrayfun(@(s) fisher(p_permute_mixdn(1:s,:)),uniq_sizes_mixdn,'UniformOutput',false);
            end
            
            
        case {'stouffer','reporter'}
            
            % Reporter background test statistic is calculated in the same way
            % as Stouffer
            bg_stat_nondir = arrayfun(@(s) stouffer(p_permute(1:s,:)),uniq_sizes,'UniformOutput',false);
            
            if ~isempty(dirs)
                bg_stat_mixup = arrayfun(@(s) stouffer(p_permute_mixup(1:s,:)),uniq_sizes_mixup,'UniformOutput',false);
                bg_stat_mixdn = arrayfun(@(s) stouffer(p_permute_mixdn(1:s,:)),uniq_sizes_mixdn,'UniformOutput',false);
                
                bg_stat_distup = arrayfun(@(s) stouffer_dir(p_permute(1:s,:),dir_permute(1:s,:)),uniq_sizes,'UniformOutput',false);
                bg_stat_distup = reshape([bg_stat_distup{:}],nperms,length(uniq_sizes))';
                
                bg_stat_distdn = -bg_stat_distup;
            end
            
            
        case 'wilcoxon'
            
            bg_stat_nondir = arrayfun(@(s) sum(p_permute(1:s,:),1),uniq_sizes,'UniformOutput',false);
            
            if ~isempty(dirs)
                bg_stat_mixup = arrayfun(@(s) sum(p_permute_mixup(1:s,:),1),uniq_sizes_mixup,'UniformOutput',false);
                bg_stat_mixdn = arrayfun(@(s) sum(p_permute_mixdn(1:s,:),1),uniq_sizes_mixdn,'UniformOutput',false);
                
                bg_stat_distup = arrayfun(@(s) sum(p_permute_distup(1:s,:),1),uniq_sizes,'UniformOutput',false);
                bg_stat_distup = reshape([bg_stat_distup{:}],nperms,length(uniq_sizes))';
                
                bg_stat_distdn = arrayfun(@(s) sum(p_permute_distdn(1:s,:),1),uniq_sizes,'UniformOutput',false);
                bg_stat_distdn = reshape([bg_stat_distdn{:}],nperms,length(uniq_sizes))';
            end
            
            
        case 'gsea'
            
            pvals_rep = repmat(pvals_distup,1,nperms);
            bg_stat_distup = arrayfun(@(s) gsea_multi(pvals_rep,ismember(p_permute,1:s)),uniq_sizes,'UniformOutput',false);
            bg_stat_distup = reshape([bg_stat_distup{:}],nperms,length(uniq_sizes))';
            
            bg_stat_distdn = -bg_stat_distup;
            
            % NOTE: non-directional and mixed-directional classes are not available for GSEA
            bg_stat_nondir = {zeros(nperms,length(uniq_sizes))};  % zeros as placeholder
            bg_stat_mixup = {zeros(nperms,length(uniq_sizes_mixup))};  % zeros as placeholder
            bg_stat_mixdn = {zeros(nperms,length(uniq_sizes_mixdn))};  % zeros as placeholder
            
    end
    
    % reshape background distribution test statistic matrices
    bg_stat_nondir = reshape([bg_stat_nondir{:}],nperms,length(uniq_sizes))';
    if ~isempty(dirs)
        bg_stat_mixup = reshape([bg_stat_mixup{:}],nperms,length(uniq_sizes_mixup))';
        bg_stat_mixdn = reshape([bg_stat_mixdn{:}],nperms,length(uniq_sizes_mixdn))';
    end
    
    
    % calculate significance of test statistics and adjust p-values for FDR
    % Note: this needs to be performed differently for GSEA (INCOMPLETE)
    GS_pval_nondir = arrayfun(@(set_stat,set_size) (1 + sum(bg_stat_nondir(ismember(uniq_sizes,set_size),:) >= set_stat))/(nperms + 1),statvals_nondir,GSsizes);
    GS_padj_nondir = adjust_pvalues(GS_pval_nondir,'benjamini',1);
    if ~isempty(dirs)
        GS_pval_mixup = arrayfun(@(set_stat,set_size) (1 + sum(bg_stat_mixup(ismember(uniq_sizes_mixup,set_size),:) >= set_stat))/(nperms + 1),statvals_mixup,GSsizes_mixup);
        GS_pval_mixdn = arrayfun(@(set_stat,set_size) (1 + sum(bg_stat_mixdn(ismember(uniq_sizes_mixdn,set_size),:) >= set_stat))/(nperms + 1),statvals_mixdn,GSsizes_mixdn);
        
        GS_padj_mixup = adjust_pvalues(GS_pval_mixup,'benjamini',1);
        GS_padj_mixdn = adjust_pvalues(GS_pval_mixdn,'benjamini',1);
    end
    
    if ismember(lower(method),{'stouffer','wilcoxon','mean','median','gsea'}) && ~isempty(dirs)
        GS_pval_distup = arrayfun(@(set_stat,set_size) (1 + sum(bg_stat_distup(ismember(uniq_sizes,set_size),:) >= set_stat))/(nperms + 1),statvals_distup,GSsizes);
        GS_pval_distdn = arrayfun(@(set_stat,set_size) (1 + sum(bg_stat_distdn(ismember(uniq_sizes,set_size),:) >= set_stat))/(nperms + 1),statvals_distdn,GSsizes);
        
        GS_padj_distup = adjust_pvalues(GS_pval_distup,'benjamini',1);
        GS_padj_distdn = adjust_pvalues(GS_pval_distdn,'benjamini',1);
    end
    
    if strcmpi(method,'reporter')
        % Normalize gene-set test statistics by subtracting the mean and
        % dividing by the standard deviation of the test statistics of
        % randomized gene sets of the same size.
        statvals_nondir = arrayfun(@(set_stat,set_size) (set_stat - mean(bg_stat_nondir(ismember(uniq_sizes,set_size),:)))/std(bg_stat_nondir(ismember(uniq_sizes,set_size),:)),statvals_nondir,GSsizes);
        
        if ~isempty(dirs)
            statvals_mixup = arrayfun(@(set_stat,set_size) (set_stat - mean(bg_stat_mixup(ismember(uniq_sizes_mixup,set_size),:)))/std(bg_stat_mixup(ismember(uniq_sizes_mixup,set_size),:)),statvals_mixup,GSsizes_mixup);
            statvals_mixdn = arrayfun(@(set_stat,set_size) (set_stat - mean(bg_stat_mixdn(ismember(uniq_sizes_mixdn,set_size),:)))/std(bg_stat_mixdn(ismember(uniq_sizes_mixdn,set_size),:)),statvals_mixdn,GSsizes_mixdn);
            
            statvals_distup = arrayfun(@(set_stat,set_size) (set_stat - mean(bg_stat_distup(ismember(uniq_sizes,set_size),:)))/std(bg_stat_distup(ismember(uniq_sizes,set_size),:)),statvals_distup,GSsizes);
            statvals_distdn = arrayfun(@(set_stat,set_size) (set_stat - mean(bg_stat_distdn(ismember(uniq_sizes,set_size),:)))/std(bg_stat_distdn(ismember(uniq_sizes,set_size),:)),statvals_distdn,GSsizes);
        end
        
        % get p-values of the normalized test statistics from the normal distribution
        GS_pval_nondir = normcdf(-statvals_nondir);
        
        if ~isempty(dirs)
            GS_pval_mixup = normcdf(-statvals_mixup);
            GS_pval_mixdn = normcdf(-statvals_mixdn);
            
            GS_pval_distup = normcdf(-statvals_distup);
            GS_pval_distdn = normcdf(-statvals_distdn);
        end
        
        % calculate FDR-adjusted p-values
        GS_padj_nondir = adjust_pvalues(GS_pval_nondir,'benjamini',1);
        
        if ~isempty(dirs)
            GS_padj_mixup = adjust_pvalues(GS_pval_mixup,'benjamini',1);
            GS_padj_mixdn = adjust_pvalues(GS_pval_mixdn,'benjamini',1);
            
            GS_padj_distup = adjust_pvalues(GS_pval_distup,'benjamini',1);
            GS_padj_distdn = adjust_pvalues(GS_pval_distdn,'benjamini',1);
        end
    end
    
elseif strcmpi(nperms,'nulldist')
   
    % calculate significance based on the null distribution
    fprintf('Calculating significance based on the null distribution... ');
    
    switch lower(method)
        
        case 'fisher'
            
            % calculate p-values based on the chi^2 distribution
            GS_pval_nondir = chi2cdf(statvals_nondir,GSsizes*2,'upper');
            GS_padj_nondir = adjust_pvalues(GS_pval_nondir,'benjamini',1);
            
            if ~isempty(dirs)
                GS_pval_mixup = chi2cdf(statvals_mixup,GSsizes_mixup*2,'upper');
                GS_pval_mixdn = chi2cdf(statvals_mixdn,GSsizes_mixdn*2,'upper');
                
                GS_padj_mixup = adjust_pvalues(GS_pval_mixup,'benjamini',1);
                GS_padj_mixdn = adjust_pvalues(GS_pval_mixdn,'benjamini',1);
            end
            
        case 'stouffer'
            
            % calculate p-values based on the normal distribution
            GS_pval_nondir = normcdf(statvals_nondir,'upper');
            GS_padj_nondir = adjust_pvalues(GS_pval_nondir,'benjamini',1);
            
            if ~isempty(dirs)
                GS_pval_distup = normcdf(statvals_distup,'upper');
                GS_pval_distdn = normcdf(statvals_distdn,'upper');
                
                GS_padj_distup = adjust_pvalues(GS_pval_distup,'benjamini',1);
                GS_padj_distdn = adjust_pvalues(GS_pval_distdn,'benjamini',1);
                
                GS_pval_mixup = normcdf(statvals_mixup,'upper');
                GS_pval_mixdn = normcdf(statvals_mixdn,'upper');
                
                GS_padj_mixup = adjust_pvalues(GS_pval_mixup,'benjamini',1);
                GS_padj_mixdn = adjust_pvalues(GS_pval_mixdn,'benjamini',1);
            end
            
   
        otherwise
            error('The null distribution significance estimation is not available for the selected stat method.');
    end
    
else
    error('Invalid NPERMS entry. Must be an integer, or "nulldist".');
end

fprintf('Done.\n');


%% Organize results

% map mixed-directional results to full list of gene sets
if ~isempty(dirs)
    sv_up = zeros(size(statvals_nondir));
    sv_dn = zeros(size(statvals_nondir));
    p_up = ones(size(GS_pval_nondir)); pa_up = p_up;
    p_dn = ones(size(GS_pval_nondir)); pa_dn = p_dn;
    
    [~,ind] = ismember(GSnums_mixup,GSnums);
    sv_up(ind) = statvals_mixup;
    statvals_mixup = sv_up;
    p_up(ind) = GS_pval_mixup;
    GS_pval_mixup = p_up;
    pa_up(ind) = GS_padj_mixup;
    GS_padj_mixup = pa_up;
    
    [~,ind] = ismember(GSnums_mixdn,GSnums);
    sv_dn(ind) = statvals_mixdn;
    statvals_mixdn = sv_dn;
    p_dn(ind) = GS_pval_mixdn;
    GS_pval_mixdn = p_dn;
    pa_dn(ind) = GS_padj_mixdn;
    GS_padj_mixdn = pa_dn;
end


% collect results
GSAres = [{'GS_name', 'GS_size', 'stat_nondir', 'p_nondir', 'padj_nondir'};
    GSnames,num2cell([GSsizes, statvals_nondir, GS_pval_nondir, GS_padj_nondir])];

if ~isempty(dirs) || strcmpi(method,'gsea')
    
    GSAres = [GSAres, [{'stat_mixup', 'p_mixup', 'padj_mixup', ...
        'stat_mixdn', 'p_mixdn', 'padj_mixdn'};
        num2cell([statvals_mixup, GS_pval_mixup, GS_padj_mixup, ...
        statvals_mixdn, GS_pval_mixdn, GS_padj_mixdn])]];
    
    if ismember(lower(method),{'reporter','stouffer','wilcoxon','mean','median','gsea'})
        GSAres = [GSAres, [{'stat_distup', 'p_distup', 'padj_distup', ...
            'stat_distdn', 'p_distdn', 'padj_distdn'};
            num2cell([statvals_distup, GS_pval_distup, GS_padj_distup, ...
            statvals_distdn, GS_pval_distdn, GS_padj_distdn])]];
    end

end


end



%% Additional functions

function s = fisher(p)
% Fisher's method
s = -2*sum(log(p),1);
end

function s = stouffer(p)
% Stouffer's method
s = sum(-norminv(p),1)./sqrt(size(p,1));
end

function s = stouffer_dir(p,d)
% Stouffer's method, with directionality
s = sum(-d.*norminv(0.5*p),1)./sqrt(size(p,1));
end

function s = gsea(p,inGS)
% GSEA method
% note: p = abs(pvals), inGS = logical indicating which indices are in GS
p(~inGS) = 1/(length(inGS) - sum(inGS));
dev = cumsum(p.*inGS)./sum(p(inGS)) - cumsum(p.*~inGS);
[~,max_ind] = max(abs(dev));
s = dev(max_ind);
end

function s = gsea_multi(p,inGS)
% GSEA method (for multiple inputs simultaneously)
% note: p and inGS are both matrices of the same dimension
p(~inGS) = 1/(size(inGS,1)-sum(inGS(:,1)));
dev = cumsum(p.*inGS)./sum(p.*inGS) - cumsum(p.*~inGS);
[~,max_ind] = max(abs(dev));
s = dev(max_ind + (0:size(inGS,1):(numel(inGS)-1)));
end

