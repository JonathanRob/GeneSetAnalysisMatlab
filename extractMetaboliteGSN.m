function gsn = extractMetaboliteGSN(model,includeComps,outfile)
% extractMetaboliteGSN  Extract metabolite geneset-geneset interaction file
% from a GEM.
%
% Construct a getset-geneset interaction (GSN) file from a genome-scale
% metabolic model (GEM).
%
%
% Usage:
%
%   gsn = extractMetaboliteGSN(model,includeComps,outfile);
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
%   outfile       File name to which the GSN will be written. See the
%                 "exportGSC" function for more detail. 
%                 (opt, Default = No file will be written)
%
%
% Output:
%
%   gsn          geneset-geneset interactions as a 2-column cell array.
%                Both columns contain the names of the metabolite gene sets
%
%



if nargin < 2 || isempty(includeComps)
    includeComps = false;
end
if nargin < 3
    outfile = [];
end

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

% generate metabolite geneset-geneset interaction array
Sbin = (model.S ~= 0);
Smets = triu(Sbin*Sbin' ~= 0) - eye(length(metNames));
[indx1, indx2] = find(Smets);
gsn = metNames([indx1, indx2]);

% remove duplicate and self-interactions
uniqMetNames = unique(metNames);
[~,gsn_indx] = ismember(gsn, uniqMetNames);
gsn_indx = unique(sort(gsn_indx, 2), 'rows');
gsn_indx(gsn_indx(:,1) == gsn_indx(:,2), :) = [];
gsn = uniqMetNames(gsn_indx);

% write to file if requested
if ~isempty(outfile)
    fprintf('Writing GSN to file... ');
    exportGSC(gsn,outfile,[]);
    fprintf('Done.\n');
end

