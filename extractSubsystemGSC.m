function gsc = extractSubsystemGSC(model,exclude)
%extractSubsystemGSC  Extract subsystem gene sets from a GEM.
%
% Construct a get set collection from a genome-scale metabolic model (GEM),
% where each gene set is a subsystem in the model.
% 
% NOTE: This function currently only works with models that have a 1:1
%       relationship between reactions and subsystems.
%
%
% Usage:
%
%   gsc = extractSubsystemGSC(model,exclude);
%
%
% Input:
%
%   model     Model structure containing gene-reaction associations.
%
%   exclude   A list of subsystems to exclude (case insensitive).
%             For example "Artificial reactions" or "Pool reactions".
%             (Opt, default none).
%
%
% Output:
%
%   gsc      Gene set collection information as a 2-column cell array.
%            The first column contains the names of the gene sets
%            (subsystems), and the second column contains the names of
%            the genes associated with each gene set.
%
%
% Jonathan Robinson, 2020-01-13


if nargin < 2
    exclude = [];
end

% get subsystem names
if iscell(model.subSystems{1})
    % need to convert to array of strings if it's an array of cells
    subSystems = cellfun(@(s) s, model.subSystems);
else
    subSystems = model.subSystems;
end

% obtain reaction-gene association matrix from model
rxnGeneMat = model.rxnGeneMat;

% exclude subsystems, if requested
if ~isempty(exclude)
    if any(~ismember(lower(exclude),lower(subSystems)))
        error('One more more subsystems to be excluded are not present in the model.');
    end
    rem_ind = ismember(lower(subSystems),lower(exclude));
    subSystems(rem_ind) = [];
    rxnGeneMat(rem_ind,:) = [];
end

% map subsystems to genes
[rxn_ind,gene_ind] = find(rxnGeneMat);
[uniq_subsys,~,subsys_groups] = unique(subSystems);
subsys_ind = subsys_groups(rxn_ind);

% get unique associations
S2G_ind = unique([subsys_ind,gene_ind],'rows');

% convert back to names
gsc = [uniq_subsys(S2G_ind(:,1)),model.genes(S2G_ind(:,2))];

% print GSC stats
fprintf('Gene set collection contains %u gene sets and %u unique genes.\n', ...
    length(unique(gsc(:,1))),length(unique(gsc(:,2))));


