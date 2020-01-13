function gsc = importGSC(filename)
%importGSC  Import a gene set collection (GSC) from a .gmt file.
%
% Usage:
%
%   gsc = importGSC(filename);
%
%
% Input:
%
%   filename  A string specifying the gene set collection (.gmt) file name.
%             Also include the path to the file if it is not in the current
%             working directory.
%
% Output:
%
%   gsc       Gene set collection information as a 2-column cell array.
%             The first column contains the names of the gene sets, and the 
%             second column contains the names of the genes associated with
%             each gene set.
%
%
% Jonathan Robinson, 2020-01-13


% import data
fid = fopen(filename);
x = textscan(fid,'%s','Delimiter','\n');
fclose(fid);

% process data
z = cellfun(@make_set,x{1},'UniformOutput',false);
gsc = vertcat(z{:});

% replace all underscores "_" in set names with dashes "-"
% (this is to avoid formatting issues if the set names are plotted)
gsc(:,1) = strrep(gsc(:,1),'_','-');

fprintf('Gene set collection contains %u gene sets and %u unique genes.\n',length(unique(gsc(:,1))),length(unique(gsc(:,2))));

end


% make_set sub-function
function gene_set = make_set(set_string)
    s = strsplit(set_string,'\t');
    gene_set = [repmat(s(1),length(s)-2,1),s(3:end)'];
end

