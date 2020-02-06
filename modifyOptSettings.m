function newOpt = modifyOptSettings(opt,newSettings)
% parseOptSettings  Assign parameter settings to options structure.
%
% Usage:
%
%   newOpt = modifyOptSettings(opt,newSettings);
%
% Input:
%
%   opt           Options structure with default settings.
%
%   newSettings   Cell array of new settings provided as NAME, VALUE pairs.
%
% Output:
%
%   newOpt        Options structure updated with new settings.
%
% 
% Jonathan Robinson, 2020-02-06


% confirm that varargin has an even number of elements
if mod(numel(newSettings),2) ~= 0
    error('Options must be provided as "NAME", "VALUE" pairs.');
end

% extract settings and confirm that option names are valid
validNames = fieldnames(opt);
optNames = lower(newSettings(1:2:end));
optVals = newSettings(2:2:end);
invalidNames = setdiff(optNames,validNames);
if ~isempty(invalidNames)
    error('"%s" is not a valid option.\n',invalidNames{:});
end

% assign new settings
for i = 1:numel(optNames)
    opt.(optNames{i}) = optVals{i};
end
newOpt = opt;

