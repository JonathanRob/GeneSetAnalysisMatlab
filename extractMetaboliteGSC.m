function gsc = extractMetaboliteGSC(model,includeComps,outfile)
%extractMetaboliteGSC  Extract metabolite gene sets from a GEM.
%
% Construct a get set collection (GSC) from a genome-scale metabolic model
% (GEM), where each gene set is a metabolite in the model.
%
%
% Usage:
%
%   gsc = extractMetaboliteGSC(model,includeComps,outfile);
%
%
% Input:
%
%   model         Model structure containing gene-reaction associations.
%
%   includeComps  Logical indicating whether metabolite names should
%                 include compartment.
%                 (opt, Default = FALSE)
%
%   outfile       File name to which the GSC will be written. See the
%                 "exportGSC" function for more detail. 
%                 (opt, Default = No file will be written)
%
%
% Output:
%
%   gsc          Gene set collection information as a 2-column cell array.
%                The first column contains the names of the gene sets
%                (metabolites), and the second column contains the names of
%                the genes associated with each gene set.
%
%
% Jonathan Robinson, 2020-01-13


if nargin < 2 || isempty(includeComps)
    includeComps = false;
end
if nargin < 3
    outfile = [];
end

% obtain met-gene association matrix
M2Gmat = logical(logical(model.S) * model.rxnGeneMat);

% obtain indices of each met-gene pair
[met_ind,gene_ind] = find(M2Gmat);

% add compartments to metabolite names if requested
if includeComps
    if all(isfield(model,{'comps','metComps'}))
        % compartment info is in "comps" and "metComps" fields
        metNames = strcat(model.metNames,'[',model.comps(model.metComps),']');
    elseif all(endsWith(model.mets,']'))
        % metIDs end in [compartment]; e.g. '[c]'
        metComps = regexp(model.mets,'\[\w+\]$','match');
        metNames = strcat(model.metNames, [metComps{:}]');
    elseif ~any(cellfun(@isempty, regexp(model.mets,'_\w{1,2}$')))
        % metIDs end in _compartment; e.g. '_c'
        metComps = regexp(model.mets,'_\w{1,2}$','match');
        metComps = regexprep(metComps,'_','');
        metNames = strcat(model.metNames,'[',metComps,']');
    else
        error('Could not find metabolite compartment information.');
    end
else
    metNames = model.metNames;
end

% generate metabolite-gene association array
gsc = [metNames(met_ind), model.genes(gene_ind)];

% remove duplicated rows
[~,G2Mnum] = ismember(gsc,gsc);
[~,keep] = unique(G2Mnum,'rows');
gsc = gsc(keep,:);

% print GSC stats
fprintf('Gene set collection contains %u gene sets and %u unique genes.\n', ...
    length(unique(gsc(:,1))),length(unique(gsc(:,2))));

% write to file if requested
if ~isempty(outfile)
    fprintf('Writing GSC to file... ');
    exportGSC(gsc,outfile);
    fprintf('Done.\n');
end


