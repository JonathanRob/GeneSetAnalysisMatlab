function [] = exportGSC(gsc,filename,header)
%exportGSC  Export a gene set collection (GSC) to a .txt or .gmt file
%
% Usage:
%
%   exportGSC(gsc,filename);
%
%
% Input:
%
%   gsc       Gene set collection information as a 2-column cell array.
%             The first column contains the names of the gene sets, and the 
%             second column contains the names of the genes associated with
%             each gene set.
%
%   filename  A string specifying the name of the file to which the GSC
%             will be written.
%
%             If the filename extension is not provided or is '.gmt', it
%             will be written in the Gene Matrix Transposed (GMT) format.
%
%             If the filename extension is '.txt', it will be written as a
%             two-column tab-separated .txt file where first column
%             contains the names of the gene sets, and the second column
%             contains the genes associated with each gene set. The file
%             will also contain a header line with column tiles "gene set"
%             and "gene".
%
%   header    An Nx2 cell array of strings to be used as column headers for
%             the output file if the filename extension is '.txt'. Input
%             empty brackets ([]) for no header.
%             (Opt, Default = {'gene set','gene'} if filename extension is
%             '.txt', otherwise Default is no header).
%
%
% Jonathan Robinson, 2020-02-04


if nargin < 3
    header = {'gene set','gene'};
end

% check filename extension
filetype = 'gmt';
if ~contains(filename,'.')
    filename = strcat(filename,'.gmt');
elseif endsWith(lower(filename),'.txt')
    filetype = 'txt';
elseif ~endsWith(lower(filename),'.gmt')
    nameparts = strsplit(filename,'.');
    error('Filename extension "%s" not recognized. Valid options are ".gmt" or ".txt".',...
        nameparts{end});
end

% write to file
if strcmp(filetype,'txt')
    writecell([header; gsc], filename, 'Delimiter', '\t');
else
    % compress associated genes into nested cells
    [uniq_gs,~,gs_num] = unique(gsc(:,1));
    gsc_compress = [uniq_gs, arrayfun(@(i) gsc(gs_num==i,2), unique(gs_num), 'UniformOutput', false)];
    
    % write file line by line
    fid = fopen(filename,'w');
    for i = 1:numel(uniq_gs)
        % "gene set name; gene set description ('na'); gene1; gene2; ..."
        gs = [gsc_compress(i,1); {'na'}; gsc_compress{i,2}];
        formatSpec = regexprep(repmat('%s\t',1,numel(gs)),'\\t$','\\n');
        fprintf(fid,formatSpec,gs{:});
    end
    fclose(fid);
end


