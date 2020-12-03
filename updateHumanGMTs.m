function [] = updateHumanGMTs(model)
%updateHumanGMTs  Create updated gene set collection files from Human-GEM.
%
% NOTE: Requires functions from the Human-GEM GitHub repository!
%
% This function is used to update the gene set collection (GSC) files
% (.gmt files) derived from Human-GEM stored in the gsc/ subdirectory.
% GSC files will be generated for metabolites (with and without
% compartments) and subsystems, and for Ensembl (ENSG) IDs and gene
% symbols. The version of the Human-GEM model used will be included in the
% output filenames.
%


%% Ensembl IDs
% Human-GEM gene IDs are by default Ensembl IDs (ENSG)

% metabolites without compartments
outfile = ['gsc/HumanGEM_v' model.version '_ensembl_metabolites.gmt'];
extractMetaboliteGSC(model, false, outfile);

% metabolites with compartments
outfile = ['gsc/HumanGEM_v' model.version '_ensembl_metabolites_comps.gmt'];
extractMetaboliteGSC(model, true, outfile);

% subsystems (excluding some uninformative subsystems)
outfile = ['gsc/HumanGEM_v' model.version '_ensembl_subsystems.gmt'];
extractSubsystemGSC(model, {'Artificial reactions','Pool reactions','Miscellaneous','Isolated'}, outfile);


%% Gene symbols

% Convert Human-GEM gene IDs to gene symbols
% NOTE: Requires that the Human-GEM GitHub repository is on the MATLAB Path
[model.grRules, model.genes, model.rxnGeneMat] = translateGrRules(model.grRules, 'Name', 'ENSG');

% metabolites without compartments
outfile = ['gsc/HumanGEM_v' model.version '_symbols_metabolites.gmt'];
extractMetaboliteGSC(model, false, outfile);

% metabolites with compartments
outfile = ['gsc/HumanGEM_v' model.version '_symbols_metabolites_comps.gmt'];
extractMetaboliteGSC(model, true, outfile);

% subsystems (excluding some uninformative subsystems)
outfile = ['gsc/HumanGEM_v' model.version '_symbols_subsystems.gmt'];
extractSubsystemGSC(model, {'Artificial reactions','Pool reactions','Miscellaneous','Isolated'}, outfile);





