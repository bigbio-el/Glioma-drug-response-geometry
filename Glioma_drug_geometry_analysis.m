function gbm_drug_geometry_analysis_reviewer_ready

clear; clc; close all;
rng(1);

% settings
fprintf('\nGLIOMA DRUG-STATE GEOMETRY ANALYSIS\n');
fprintf('================================\n\n');

cfg.minCellsPatientControl = 100;
cfg.minCellsPerState       = 25;
cfg.minStatesForMetric     = 3;
cfg.confidentMargin        = 0.25;
cfg.minPatientsReplication = 3;
cfg.nBootstrap             = 2000;
cfg.eps                    = 1e-8;
cfg.nNullDraws             = 1000;
cfg.nullSeed               = 20260914;

cfg.nIndependentHVG        = 500;
cfg.hvgMinDetection        = 0.02;
cfg.hvgCellsPerRun         = 75;
cfg.hvgMaxRunsPerPatient   = 2;
cfg.hvgCacheVersion        = 'HVG_VALIDATION_V2_TISSUEID_20260908';

cfg.neftelCrossfitSplits       = 10;
cfg.neftelCrossfitSeed         = 19062019;
cfg.neftelMinGenesPerHalf      = 12;
cfg.neftelMinValidFoldFraction = 0.75;
cfg.neftelCacheVersion         = 'NEFTEL_CROSSFIT_V2_TISSUEID_20260908';

marker.AC  = {'CST3','GFAP','S100B','HOPX','SLC1A3','MLC1','AQP4'};
marker.OPC = {'PLP1','ALCAM','OLIG1','OLIG2','OMG','PLLP','PDGFRA'};
marker.NPC = {'SOX4','SOX11','DCX','CD24','DLL3','STMN2','STMN4', 'RND3','DLX5','DLX6-AS1'};
marker.MES = {'VIM','ANXA1','ANXA2','CHI3L1','CD44','DDIT3', 'ADM','HILPDA','ENO2','LDHA'};
states = {'AC','OPC','NPC','MES'};
featureGenes = unique([marker.AC marker.OPC marker.NPC marker.MES], 'stable');

% input
[fileName,filePath] = uigetfile({'*.loom','Loom files (*.loom)'}, 'Select the integrated glioma loom file');
if isequal(fileName,0)
    error('No loom file was selected.');
end
loomFile = fullfile(filePath,fileName);
fprintf('Input: %s\n\n', loomFile);

meta = inspectAndDetectLoom(loomFile, featureGenes);
fprintf('Detected gene field:      %s\n', meta.genePath);
fprintf('Detected sample field:   %s\n', meta.patientPath);
fprintf('Detected treatment field: %s\n', meta.treatmentPath);
fprintf('Detected cell-type field: %s\n', meta.celltypePath);
if ~isempty(meta.librarySizePath)
    fprintf('Detected library-size field: %s\n', meta.librarySizePath);
else
    fprintf('Library-size field:       not found (total counts will be computed from full matrix)\n');
end
fprintf('\n');

% metadata
genes = readTextDataset(loomFile, meta.genePath, meta.nGenes);
patientRaw   = readTextDataset(loomFile, meta.patientPath, meta.nCells);
treatmentRaw = readTextDataset(loomFile, meta.treatmentPath, meta.nCells);
celltypeRaw  = readTextDataset(loomFile, meta.celltypePath, meta.nCells);

try
    tissueRaw = readTextDataset(loomFile,'/col_attrs/TissueID',meta.nCells);
catch
    error('The released TissueID field is required for biological-unit inference.');
end

pwLabel   = normalizePatient(patientRaw);
patient   = upper(strtrim(string(tissueRaw)));
treatment = normalizeTreatment(treatmentRaw);
celltype  = lower(strtrim(celltypeRaw));

malignant = contains(celltype,'transform') | contains(celltype,'malig') | contains(celltype,'neoplast') | contains(celltype,'tumor') | contains(celltype,'glioma');

if mean(malignant) > 0.95
    error(['The automatically selected cell-type column marks >95%% of cells as tumor-like. ' 'This is suspicious. Run inspect_glioma_loom.m and send its output before analysis.']);
end
if nnz(malignant) < 500
    error(['Fewer than 500 transformed/neoplastic cells were detected. ' 'Run inspect_glioma_loom.m and send its output before analysis.']);
end

validMeta = patient ~= "" & treatment ~= "";
malignant = malignant & validMeta;

discoveryIDs=unique(patient(malignant & pwLabel=="PW030"));
discoveryIDs=discoveryIDs(discoveryIDs~="");
if numel(discoveryIDs)~=1
    error('PW030 maps to %d TissueID values; expected exactly one.',numel(discoveryIDs));
end
discoveryPatient=discoveryIDs(1);

fprintf('Cells in loom:                    %d\n', meta.nCells);
fprintf('Transformed cells kept:           %d\n', nnz(malignant));
fprintf('Biological TissueID units:        %d\n', numel(unique(patient(malignant))));
fprintf('Sequencing-name PW prefixes:      %d\n', numel(unique(pwLabel(malignant))));
fprintf('PW030 biological TissueID:        %s\n', discoveryPatient);
fprintf('Conditions represented:           %d\n\n', numel(unique(treatment(malignant))));

[foundGenes, geneIndex, missingGenes] = mapGenes(geneListUpper(genes), featureGenes);
if ~isempty(missingGenes)
    fprintf('Missing marker genes (%d): %s\n', numel(missingGenes), strjoin(missingGenes, ', '));
end
if numel(foundGenes) < 20
    error('Too few GBM state marker genes were found (%d).', numel(foundGenes));
end
fprintf('Marker genes found:       %d/%d\n\n', numel(foundGenes), numel(featureGenes));

X = readSelectedGenesFromLoom(loomFile, meta, geneIndex);
X = double(X);

if ~isempty(meta.librarySizePath)
    lib = readNumericDataset(loomFile, meta.librarySizePath);
    lib = double(lib(:));
    if numel(lib) ~= meta.nCells || any(~isfinite(lib)) || median(lib(lib>0)) <= 0
        warning('Library-size field was unusable. Recomputing total counts from the full loom matrix.');
        lib = getOrComputeLibrarySize(loomFile, meta);
    end
else
    lib = getOrComputeLibrarySize(loomFile, meta);
end
lib = double(lib(:));
lib(lib<=0 | ~isfinite(lib)) = 1;
Xlog = log1p(1e4 .* X ./ lib);

Z = nan(size(Xlog));
uniquePatients = unique(patient(malignant));
keepPatient = false(size(patient));

for ip = 1:numel(uniquePatients)
    p = uniquePatients(ip);
    idxAll = malignant & patient==p;
    idxCtrl = idxAll & treatment=="vehicle";
    if nnz(idxCtrl) < cfg.minCellsPatientControl
        continue;
    end
    mu = mean(Xlog(idxCtrl,:),1,'omitnan');
    sd = std(Xlog(idxCtrl,:),0,1,'omitnan');
    sd(sd < 1e-6 | ~isfinite(sd)) = 1;
    Z(idxAll,:) = (Xlog(idxAll,:) - mu) ./ sd;
    keepPatient(idxAll) = true;
end

analysisMask = malignant & keepPatient;
fprintf('Transformed cells with usable paired control reference: %d\n', nnz(analysisMask));
fprintf('Usable TissueID units: %d\n\n', numel(unique(patient(analysisMask))));

Score = nan(meta.nCells,4);
for k = 1:4
    gset = marker.(states{k});
    cols = find(ismember(upper(string(foundGenes)), upper(string(gset))));
    if numel(cols) < 3
        error('State %s has only %d available markers.', states{k}, numel(cols));
    end
    Score(:,k) = mean(Z(:,cols),2,'omitnan');
end

[bestScore,bestState] = max(Score,[],2);
Sorted = sort(Score,2,'descend');
margin = Sorted(:,1)-Sorted(:,2);
bestState(~analysisMask) = 0;
bestScore(~analysisMask) = NaN;
margin(~analysisMask) = NaN;

R = computeGeometryTable(patient,treatment,analysisMask,bestState,margin,Z, states,cfg,false);
Rconf = computeGeometryTable(patient,treatment,analysisMask,bestState,margin,Z, states,cfg,true);

if isempty(R)
    error('No biological-unit-drug pairs passed the minimum cell/state thresholds.');
end

rngState=rng;
rng(cfg.nullSeed);
fprintf('\nCalibrating C, K and M against within-unit null models (%d draws each)\n',cfg.nNullDraws);
NullSplit=computeVehicleSplitNull(patient,treatment,analysisMask,bestState,Z,R,states,cfg);
NullPerm=computeLabelPermutationNull(patient,treatment,analysisMask,bestState,Z,R,states,cfg);
NullSummary=summarizeNullCalibration(NullSplit,NullPerm);
CrossDrug=computeCrossDrugCoherence(patient,treatment,analysisMask,bestState,Z,R,states,cfg);
CrossDrugSummary=summarizeCrossDrugCoherence(CrossDrug);
[CrossDrugUnit,CrossDrugGlobal]=summarizeCrossDrugByUnit(CrossDrug,cfg);
RotNull=computeRotationNull(patient,treatment,analysisMask,bestState,Z,R,states,cfg);
RotNullSummary=summarizeRotationNull(RotNull);
DrugStratifiedSign=computeDrugStratifiedSignTest(R);
rng(rngState);

[hvgGenes,hvgIndex,HVGSelection] = selectIndependentHVGs( loomFile,meta,genes,lib,patient,treatment,analysisMask, featureGenes,cfg);

Rind = computeIndependentHVGGeometry( loomFile,meta,hvgIndex,lib,patient,treatment,analysisMask, bestState,R,states,cfg);

IndSummary = summarizeAcrossPatients(Rind,cfg);
IndValidation = comparePrimaryIndependentGeometry(R,Rind);

NeftelModules = getNeftelModules2019();

[ZNeftel,NeftelGeneNames,NeftelModuleCols,NeftelCoverage] = loadNeftelModuleExpression(loomFile,meta,genes,lib,patient,treatment, analysisMask,NeftelModules,cfg);

analysisIndex = find(analysisMask);
patientA = patient(analysisIndex);
treatmentA = treatment(analysisIndex);
compactStateA = bestState(analysisIndex);

[neftelFullState,~] = assignNeftelFourStates( ZNeftel,NeftelModuleCols,true(size(ZNeftel,2),1));

RneftelFull = computeGeometryLocal( patientA,treatmentA,neftelFullState,ZNeftel,states,cfg);

NeftelFullValidation = compareGeometryByKey(R,RneftelFull);
[NeftelAssignmentSummary,NeftelAssignmentConfusion] = compareStateAssignments(compactStateA,neftelFullState,patientA,states);

[CrossfitAll,CrossfitSummary,CrossfitValidation,CrossfitAssignment] = runNeftelCrossfit(ZNeftel,NeftelModuleCols,patientA,treatmentA, compactStateA,R,states,cfg);

CrossfitEligible = CrossfitSummary(CrossfitSummary.eligible_validation,:);
if isempty(CrossfitEligible)
    CrossfitDrugSummary = table();
else
    CrossfitDrugSummary = summarizeAcrossPatients(CrossfitEligible,cfg);
end

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
outDir = fullfile(scriptDir,'results');
if ~exist(outDir,'dir'), mkdir(outDir); end

writetable(R, fullfile(outDir,'patient_drug_geometry.csv'));
writetable(Rconf, fullfile(outDir,'patient_drug_geometry_confident_cells.csv'));

writetable(NullSplit, fullfile(outDir,'null_vehicle_split.csv'));
writetable(NullPerm, fullfile(outDir,'null_label_permutation.csv'));
writetable(NullSummary, fullfile(outDir,'null_calibration_summary.csv'));
writetable(CrossDrug, fullfile(outDir,'cross_drug_coherence.csv'));
writetable(CrossDrugSummary, fullfile(outDir,'cross_drug_coherence_summary.csv'));
writetable(CrossDrugUnit, fullfile(outDir,'cross_drug_unit_summary.csv'));
writetable(CrossDrugGlobal, fullfile(outDir,'cross_drug_global_inference.csv'));
writetable(RotNull, fullfile(outDir,'null_rotation.csv'));
writetable(RotNullSummary, fullfile(outDir,'null_rotation_summary.csv'));
writetable(DrugStratifiedSign, fullfile(outDir,'drug_stratified_sign_test.csv'));

PatientAudit = buildPatientAudit(patient,patientRaw,tissueRaw,treatment, malignant,analysisMask,R);
writetable(PatientAudit, fullfile(outDir,'patient_cohort_audit.csv'));

writetable(HVGSelection, fullfile(outDir,'independent_HVG_genes.csv'));
writetable(Rind, fullfile(outDir,'independent_HVG_geometry.csv'));
writetable(IndSummary, fullfile(outDir,'independent_HVG_drug_summary.csv'));
writetable(IndValidation, fullfile(outDir,'independent_HVG_validation_summary.csv'));

writeNeftelModules(fullfile(outDir,'neftel_2019_meta_modules_used.csv'),NeftelModules);
writetable(NeftelCoverage, fullfile(outDir,'neftel_module_gene_coverage.csv'));
writetable(RneftelFull, fullfile(outDir,'neftel_full_module_geometry.csv'));
writetable(NeftelFullValidation, fullfile(outDir,'neftel_full_module_validation_summary.csv'));
writetable(NeftelAssignmentSummary, fullfile(outDir,'neftel_state_assignment_summary.csv'));
writetable(NeftelAssignmentConfusion, fullfile(outDir,'neftel_state_assignment_confusion.csv'));
writetable(CrossfitAll, fullfile(outDir,'neftel_crossfit_all_folds.csv'));
writetable(CrossfitSummary, fullfile(outDir,'neftel_crossfit_patient_drug_summary.csv'));
writetable(CrossfitValidation, fullfile(outDir,'neftel_crossfit_validation_summary.csv'));
writetable(CrossfitAssignment, fullfile(outDir,'neftel_crossfit_assignment_stability.csv'));
if ~isempty(CrossfitDrugSummary)
    writetable(CrossfitDrugSummary, fullfile(outDir,'neftel_crossfit_drug_summary.csv'));
end

Disc = R(R.patient==discoveryPatient,:);
writetable(Disc, fullfile(outDir,'PW030_discovery_geometry.csv'));

Summary = summarizeAcrossPatients(R,cfg);
writetable(Summary, fullfile(outDir,'replicated_drug_summary.csv'));

Robust = outerjoin(R,Rconf,'Keys',{'patient','drug'},'MergeKeys',true, 'LeftVariables',{'patient','drug','coherence_C','contraction_K','magnitude_M','n_states'}, 'RightVariables',{'coherence_C','contraction_K','magnitude_M','n_states'});
writetable(Robust, fullfile(outDir,'robustness_primary_vs_confident.csv'));

PairSummary = summarizeMatchedDrugPairs(R,cfg);
writetable(PairSummary, fullfile(outDir,'matched_drug_pair_summary.csv'));

RobustSummary = summarizeRobustness(Robust);
writetable(RobustSummary, fullfile(outDir,'robustness_summary.csv'));

RelationshipSummary = summarizeMetricRelationships(R);
writetable(RelationshipSummary, fullfile(outDir,'metric_relationship_summary.csv'));

MagnitudeAdjusted = summarizeMagnitudeAdjustedK(R);
writetable(MagnitudeAdjusted, fullfile(outDir,'magnitude_adjusted_K.csv'));

writeMarkerFile(fullfile(outDir,'marker_genes_used.txt'), marker, missingGenes);
writeDetectionFile(fullfile(outDir,'loom_fields_detected.txt'), meta);
writeKeyResults(fullfile(outDir,'key_results.txt'),R,Disc,Summary, PairSummary,RobustSummary,RelationshipSummary,MagnitudeAdjusted, Rind,IndSummary,IndValidation,PatientAudit, NeftelFullValidation,NeftelAssignmentSummary,CrossfitSummary, CrossfitValidation,CrossfitAssignment,NullSplit,NullPerm,NullSummary, CrossDrug,CrossDrugSummary,RotNull,RotNullSummary,DrugStratifiedSign,discoveryPatient);

fprintf('\nPRIMARY BIOLOGICAL-UNIT-DRUG METRICS\n');
fprintf('============================\n');
disp(R(:,{'patient','drug','n_states','coherence_C','contraction_K','magnitude_M','state_shift_TV'}));

fprintf('\nREPLICATED DRUG SUMMARY (>= %d TissueID units)\n',cfg.minPatientsReplication);
fprintf('=========================================\n');
disp(Summary);

fprintf('\nNULL CALIBRATION OF C, K AND M\n');
fprintf('==============================\n');
disp(NullSummary);

fprintf('\nWITHIN-UNIT CROSS-DRUG COHERENCE\n');
fprintf('================================\n');
disp(CrossDrugSummary);

fprintf('\nTISSUEID-LEVEL CROSS-DRUG SPECIFICITY\n');
fprintf('=====================================\n');
disp(CrossDrugUnit);
fprintf('\nTISSUEID-LEVEL GLOBAL INFERENCE\n');
fprintf('================================\n');
disp(CrossDrugGlobal);

fprintf('\nROTATION NULL (independent random directions, fixed norms)\n');
fprintf('============================================================\n');
disp(RotNullSummary);

fprintf('\nDRUG-STRATIFIED SIGN TEST FOR K\n');
fprintf('===============================\n');
disp(DrugStratifiedSign);

fprintf('\nROBUSTNESS SUMMARY\n');
fprintf('==================\n');
disp(RobustSummary);

fprintf('\nMATCHED-UNIT DRUG COMPARISONS (exploratory)\n');
fprintf('===============================================\n');
showPairs = PairSummary(PairSummary.metric=="contraction_K" & PairSummary.n_matched>=4,:);
disp(showPairs(:,{'drug_A','drug_B','n_matched','median_B_minus_A', 'sign_consistency','exact_sign_p'}));

fprintf('\nINDEPENDENT 500-HVG FEATURE-SPACE VALIDATION\n');
fprintf('=============================================\n');
disp(IndValidation);

fprintf('\nPUBLISHED NEFTEL META-MODULE COVERAGE\n');
fprintf('=====================================\n');
disp(NeftelCoverage);

fprintf('\nFULL NEFTEL MODULE SENSITIVITY\n');
fprintf('==============================\n');
disp(NeftelFullValidation);
disp(NeftelAssignmentSummary);

fprintf('\nNEFTEL CROSS-FITTED HELD-OUT GEOMETRY\n');
fprintf('======================================\n');
disp(CrossfitValidation);
fprintf('Median assignment stability across repeated splits:\n');
fprintf('  A vs B agreement: %.3f\n',median(CrossfitAssignment.agreement_A_B,'omitnan'));
fprintf('  A vs B kappa:     %.3f\n',median(CrossfitAssignment.kappa_A_B,'omitnan'));
fprintf('  A vs compact:     %.3f\n',median(CrossfitAssignment.agreement_compact_A,'omitnan'));
fprintf('  B vs compact:     %.3f\n',median(CrossfitAssignment.agreement_compact_B,'omitnan'));

fprintf('\nBIOLOGICAL UNIT / COHORT AUDIT\n');
fprintf('======================\n');
disp(PatientAudit(:,{'patient','pw_prefixes','n_cells_all','n_malignant','n_vehicle_malignant', 'n_nonvehicle_treatments','n_geometry_drugs','contributes_geometry'}));

fprintf('\nInterpretation:\n');
fprintf('  C near +1 : malignant states respond in similar directions.\n');
fprintf('  C near  0 : state responses are directionally heterogeneous.\n');
fprintf('  C < 0     : some state responses oppose one another.\n');
fprintf('  K < 0     : state geometry contracts after treatment.\n');
fprintf('  K > 0     : state geometry expands after treatment.\n');
fprintf('\nResults saved to:\n%s\n',outDir);
fprintf('All publication tables and sensitivity outputs are saved in the results folder.\n\n');
end

function meta = inspectAndDetectLoom(file, markerGenes)
infoM = h5info(file,'/matrix');
shape = infoM.Dataspace.Size;

rowInfo = h5info(file,'/row_attrs');
colInfo = h5info(file,'/col_attrs');
rowPaths = string.empty;
for i=1:numel(rowInfo.Datasets)
    rowPaths(end+1) = "/row_attrs/" + string(rowInfo.Datasets(i).Name);
end
colPaths = string.empty;
for i=1:numel(colInfo.Datasets)
    colPaths(end+1) = "/col_attrs/" + string(colInfo.Datasets(i).Name);
end

bestGeneScore = -inf; genePath = ""; nGenes = NaN;
for i=1:numel(rowPaths)
    try
        nGuess = max(shape);
        s = readTextDataset(file,rowPaths(i),[]);
        if isempty(s), continue; end
        hits = nnz(ismember(upper(strtrim(s)), upper(string(markerGenes))));
        nameBonus = 0;
        nm = lower(rowPaths(i));
        if contains(nm,'gene') || contains(nm,'symbol'), nameBonus=5; end
        score = hits*10 + nameBonus;
        if score > bestGeneScore
            bestGeneScore=score; genePath=rowPaths(i); nGenes=numel(s);
        end
    catch
    end
end
if genePath=="" || bestGeneScore < 20
    error('Could not identify a gene-name field under /row_attrs.');
end

if shape(1)==nGenes
    nCells=shape(2); genesOnRows=true;
elseif shape(2)==nGenes
    nCells=shape(1); genesOnRows=false;
else
    error('Matrix dimensions do not match the detected gene list.');
end

bestPatient=-inf; bestTreat=-inf; bestCell=-inf;
patientPath=""; treatmentPath=""; celltypePath="";
for i=1:numel(colPaths)
    try
        s = readTextDataset(file,colPaths(i),nCells);
        if numel(s)~=nCells, continue; end
        nm=lower(colPaths(i)); vals=lower(s);

        patHits = nnz(~cellfun(@isempty,regexp(cellstr(vals),'pw\d+','once')));
        patScore = patHits/max(1,nCells)*20;
        if contains(nm,'patient') || contains(nm,'subject') || contains(nm,'donor'), patScore=patScore+10; end
        if patScore>bestPatient, bestPatient=patScore; patientPath=colPaths(i); end

        trHits = nnz(contains(vals,'panobin') | contains(vals,'etopos') | contains(vals,'vehicle') | contains(vals,'dmso') | contains(vals,'topotec') | contains(vals,'ispin') | contains(vals,'ana-12') | contains(vals,'ro492') | contains(vals,'givinostat') | contains(vals,'tazemet'));
        trScore = trHits/max(1,nCells)*20;
        if contains(nm,'drug') || contains(nm,'treat') || contains(nm,'condition'), trScore=trScore+10; end
        if trScore>bestTreat, bestTreat=trScore; treatmentPath=colPaths(i); end

        categoryHits = nnz(contains(vals,'transform')) + nnz(contains(vals,'myeloid')) + nnz(contains(vals,'oligodend')) + nnz(contains(vals,'malig'));
        cellScore = categoryHits/max(1,nCells)*20;
        if contains(nm,'cell') || contains(nm,'type') || contains(nm,'class'), cellScore=cellScore+10; end
        if numel(unique(vals))>=2 && numel(unique(vals))<100, cellScore=cellScore+2; end
        if cellScore>bestCell, bestCell=cellScore; celltypePath=colPaths(i); end
    catch
    end
end
if patientPath=="" || bestPatient<1, error('Could not detect patient metadata.'); end
if treatmentPath=="" || bestTreat<1, error('Could not detect treatment metadata.'); end
if celltypePath=="" || bestCell<1, error('Could not detect cell-type metadata.'); end

libPath="";
for i=1:numel(colPaths)
    nm=lower(colPaths(i));
    if contains(nm,'numi') || contains(nm,'n_umi') || contains(nm,'ncount') || contains(nm,'total_count') || contains(nm,'total_umi') || contains(nm,'molecules')
        try
            x=readNumericDataset(file,colPaths(i));
            if isvector(x) && numel(x)==nCells
                libPath=colPaths(i); break;
            end
        catch
        end
    end
end

meta.matrixPath='/matrix';
meta.genePath=char(genePath);
meta.patientPath=char(patientPath);
meta.treatmentPath=char(treatmentPath);
meta.celltypePath=char(celltypePath);
meta.librarySizePath=char(libPath);
meta.nGenes=nGenes;
meta.nCells=nCells;
meta.genesOnRows=genesOnRows;
meta.rowAttrPaths=rowPaths;
meta.colAttrPaths=colPaths;
end

function s = readTextDataset(file,path,expectedN)
x = h5read(file,char(path));
if isstring(x)
    s=x(:);
elseif iscell(x)
    s=string(x(:));
elseif ischar(x)
    if isempty(expectedN)
        if size(x,2)>=size(x,1), s=string(cellstr(x')); else, s=string(cellstr(x)); end
    elseif size(x,2)==expectedN
        s=string(cellstr(x'));
    elseif size(x,1)==expectedN
        s=string(cellstr(x));
    else
        s=string(x(:));
    end
elseif isnumeric(x) && ~isvector(x) && (isa(x,'uint8') || isa(x,'int8'))
    if ~isempty(expectedN) && size(x,2)==expectedN
        c=cell(expectedN,1);
        for j=1:expectedN
            q=char(x(:,j)'); q=q(q~=0); c{j}=q;
        end
        s=string(c);
    elseif ~isempty(expectedN) && size(x,1)==expectedN
        c=cell(expectedN,1);
        for j=1:expectedN
            q=char(x(j,:)); q=q(q~=0); c{j}=q;
        end
        s=string(c);
    else
        error('Unsupported encoded string matrix at %s.',path);
    end
elseif isnumeric(x) && isvector(x)
    s=string(x(:));
else
    error('Unsupported text dataset type at %s.',path);
end
s=strtrim(s(:));
end

function lib = getOrComputeLibrarySize(file,meta)
[filePath,fileBase,~]=fileparts(file);
cacheFile=fullfile(filePath,[fileBase '_library_size_cache.mat']);

if exist(cacheFile,'file')
    q=load(cacheFile,'lib','nCells_cache','nGenes_cache');
    if isfield(q,'lib') && numel(q.lib)==meta.nCells
        if (~isfield(q,'nCells_cache') || q.nCells_cache==meta.nCells) && (~isfield(q,'nGenes_cache') || q.nGenes_cache==meta.nGenes)
            lib=double(q.lib(:));
            fprintf('Loaded cached total library sizes: %s\n',cacheFile);
            return;
        end
    end
end

fprintf(['No per-cell library-size metadata is present.\n' 'Computing total counts across all %d genes for %d cells.\n' 'This is the slow step and is done only once; progress is shown below.\n'], meta.nGenes,meta.nCells);

lib=zeros(meta.nCells,1,'double');
blockCells=500;

if meta.genesOnRows
    for firstCell=1:blockCells:meta.nCells
        nThis=min(blockCells,meta.nCells-firstCell+1);
        A=h5read(file,meta.matrixPath,[1 firstCell],[meta.nGenes nThis]);
        lib(firstCell:firstCell+nThis-1)=double(sum(A,1))';
        if firstCell==1 || mod(firstCell-1,blockCells*25)==0 || firstCell+nThis-1==meta.nCells
            fprintf('  library-size scan: %d / %d cells (%.1f%%)\n', firstCell+nThis-1,meta.nCells,100*(firstCell+nThis-1)/meta.nCells);
        end
    end
else
    for firstCell=1:blockCells:meta.nCells
        nThis=min(blockCells,meta.nCells-firstCell+1);
        A=h5read(file,meta.matrixPath,[firstCell 1],[nThis meta.nGenes]);
        lib(firstCell:firstCell+nThis-1)=double(sum(A,2));
        if firstCell==1 || mod(firstCell-1,blockCells*25)==0 || firstCell+nThis-1==meta.nCells
            fprintf('  library-size scan: %d / %d cells (%.1f%%)\n', firstCell+nThis-1,meta.nCells,100*(firstCell+nThis-1)/meta.nCells);
        end
    end
end

nCells_cache=meta.nCells;
nGenes_cache=meta.nGenes;
save(cacheFile,'lib','nCells_cache','nGenes_cache');
fprintf('Saved library-size cache: %s\n',cacheFile);
end

function x = readNumericDataset(file,path)
x=h5read(file,char(path));
if iscell(x), x=cell2mat(x); end
x=x(:);
end

function u = geneListUpper(g)
u=upper(strtrim(string(g(:))));
u=regexprep(u,'\.[0-9]+$','');
end

function [found,idx,missing] = mapGenes(geneUpper,requested)
req=string(requested(:)); found=string.empty; idx=[]; missing=string.empty;
for i=1:numel(req)
    j=find(geneUpper==upper(req(i)),1,'first');
    if isempty(j)
        missing(end+1)=req(i);
    else
        found(end+1)=req(i);
        idx(end+1)=j;
    end
end
found=cellstr(found); missing=cellstr(missing);
end

function X = readSelectedGenesFromLoom(file,meta,gidx)
nG=numel(gidx); nC=meta.nCells;
X=zeros(nC,nG,'single');
fprintf('Reading %d marker genes from /matrix\n',nG);
for q=1:nG
    gi=gidx(q);
    if meta.genesOnRows
        v=h5read(file,meta.matrixPath,[gi 1],[1 nC]);
    else
        v=h5read(file,meta.matrixPath,[1 gi],[nC 1]);
    end
    X(:,q)=single(v(:));
end
end

function p = normalizePatient(raw)
r=upper(strtrim(string(raw(:)))); p=strings(size(r));
for i=1:numel(r)
    tok=regexp(char(r(i)),'PW\d+','match','once');
    if ~isempty(tok), p(i)=string(tok); end
end
end

function t = normalizeTreatment(raw)
r=lower(strtrim(string(raw(:)))); t=strings(size(r));
for i=1:numel(r)
    s=r(i);
    if contains(s,'vehicle') || contains(s,'dmso') || s=="control" || s=="ctrl" || s=="none"
        t(i)="vehicle";
    elseif contains(s,'panobin')
        t(i)="panobinostat";
    elseif contains(s,'etopos')
        t(i)="etoposide";
    elseif contains(s,'topotec')
        t(i)="topotecan";
    elseif contains(s,'ro492') || contains(s,'r0492')
        t(i)="RO4929097";
    elseif contains(s,'tazemet')
        t(i)="tazemetostat";
    elseif contains(s,'ispin') || contains(s,'ispen')
        t(i)="ispinesib";
    elseif contains(s,'ana-12') || contains(s,'ana12')
        t(i)="ANA-12";
    elseif contains(s,'givinostat')
        t(i)="givinostat";
    elseif contains(s,'vorinostat')
        t(i)="vorinostat";
    elseif contains(s,'temozolomide') || contains(s,'tmz')
        t(i)="temozolomide";
    else
        s=regexprep(s,'[0-9\.]+\s*(nm|um|µm|mm).*','');
        s=strtrim(s);
        if strlength(s)>0, t(i)=s; end
    end
end
end

function R = computeGeometryTable(patient,treatment,mask,state,margin,Z,states,cfg,confidentOnly)
rows={};
pats=unique(patient(mask));
for ip=1:numel(pats)
    p=pats(ip);
    base=mask & patient==p;
    if confidentOnly, base=base & margin>=cfg.confidentMargin; end
    ctrl=base & treatment=="vehicle";
    if nnz(ctrl)<cfg.minCellsPatientControl, continue; end
    drugs=unique(treatment(base & treatment~="vehicle" & treatment~=""));
    for id=1:numel(drugs)
        d=drugs(id);
        tr=base & treatment==d;
        if nnz(tr)<cfg.minCellsPatientControl/2, continue; end
        C0=[]; C1=[]; n0=[]; n1=[]; used=string.empty;
        p0=[]; p1=[];
        for k=1:4
            a=ctrl & state==k; b=tr & state==k;
            if nnz(a)>=cfg.minCellsPerState && nnz(b)>=cfg.minCellsPerState
                C0(end+1,:)=mean(Z(a,:),1,'omitnan');
                C1(end+1,:)=mean(Z(b,:),1,'omitnan');
                n0(end+1)=nnz(a); n1(end+1)=nnz(b);
                used(end+1)=string(states{k});
            end
            p0(k)=nnz(a)/max(1,nnz(ctrl));
            p1(k)=nnz(b)/max(1,nnz(tr));
        end
        ns=size(C0,1);
        if ns<cfg.minStatesForMetric, continue; end
        Delta=C1-C0;
        coherence=meanPairwiseCosine(Delta,cfg.eps);
        D0=meanPairwiseDistance(C0);
        D1=meanPairwiseDistance(C1);
        contraction=log2((D1+cfg.eps)/(D0+cfg.eps));
        magnitude=mean(sqrt(sum(Delta.^2,2)),'omitnan');
        tv=0.5*sum(abs(p1-p0));
        rows(end+1,:)={p,d,nnz(ctrl),nnz(tr),ns,strjoin(used,','), coherence,contraction,magnitude,D0,D1,tv,min(n0),min(n1)};
    end
end
if isempty(rows)
    R=table(); return;
end
R=cell2table(rows,'VariableNames',{'patient','drug','n_control','n_treated', 'n_states','states_used','coherence_C','contraction_K','magnitude_M', 'distance_control','distance_treated','state_shift_TV','min_state_control','min_state_treated'});
R.patient=string(R.patient); R.drug=string(R.drug);
numvars={'n_control','n_treated','n_states','coherence_C','contraction_K','magnitude_M', 'distance_control','distance_treated','state_shift_TV','min_state_control','min_state_treated'};
for i=1:numel(numvars)
    if iscell(R.(numvars{i})), R.(numvars{i})=cell2mat(R.(numvars{i})); end
end
end

function c = meanPairwiseCosine(D,epsv)
v=[];
for i=1:size(D,1)-1
    for j=i+1:size(D,1)
        den=norm(D(i,:))*norm(D(j,:));
        if den>epsv, v(end+1)=dot(D(i,:),D(j,:))/den; end
    end
end
if isempty(v), c=NaN; else, c=mean(v); end
end

function d = meanPairwiseDistance(C)
v=[];
for i=1:size(C,1)-1
    for j=i+1:size(C,1)
        v(end+1)=norm(C(i,:)-C(j,:));
    end
end
if isempty(v), d=NaN; else, d=mean(v); end
end

function S = summarizeAcrossPatients(R,cfg)
drugs=unique(R.drug); rows={};
for i=1:numel(drugs)
    d=drugs(i); A=R(R.drug==d,:); n=height(A);
    if n<cfg.minPatientsReplication, continue; end
    C=A.coherence_C; K=A.contraction_K; M=A.magnitude_M;
    [cLo,cHi]=bootstrapMedianCI(C,cfg.nBootstrap);
    [kLo,kHi]=bootstrapMedianCI(K,cfg.nBootstrap);
    kSign=mean(sign(K)==sign(median(K,'omitnan')));
    cSign=mean(sign(C)==sign(median(C,'omitnan')));
    pK=exactSignP(K);
    rows(end+1,:)={d,n,median(C,'omitnan'),cLo,cHi,cSign, median(K,'omitnan'),kLo,kHi,kSign,pK,median(M,'omitnan')};
end
if isempty(rows)
    S=table(); return;
end
S=cell2table(rows,'VariableNames',{'drug','n_patients','median_C','C_CI_low','C_CI_high', 'C_sign_consistency','median_K','K_CI_low','K_CI_high','K_sign_consistency', 'K_exact_sign_p','median_M'});
S.drug=string(S.drug);
numvars=S.Properties.VariableNames(2:end);
for i=1:numel(numvars)
    if iscell(S.(numvars{i})), S.(numvars{i})=cell2mat(S.(numvars{i})); end
end
end

function [lo,hi]=bootstrapMedianCI(x,B)
x=x(isfinite(x)); n=numel(x);
if n<2, lo=NaN; hi=NaN; return; end
b=zeros(B,1);
for k=1:B
    b(k)=median(x(randi(n,n,1)));
end
q=prctile(b,[2.5 97.5]); lo=q(1); hi=q(2);
end

function p=exactSignP(x)
x=x(isfinite(x) & x~=0); n=numel(x);
if n==0, p=NaN; return; end
k=min(nnz(x>0),nnz(x<0));
prob=0;
for j=0:k
    prob=prob+nchoosek(n,j)*(0.5^n);
end
p=min(1,2*prob);
end

function PairSummary = summarizeMatchedDrugPairs(R,cfg)
drugs=unique(R.drug);
metrics=["contraction_K","coherence_C","magnitude_M"];
rows={};

for ia=1:numel(drugs)-1
    for ib=ia+1:numel(drugs)
        dA=drugs(ia); dB=drugs(ib);
        A=R(R.drug==dA,:);
        B=R(R.drug==dB,:);
        [commonP,idxA,idxB]=intersect(A.patient,B.patient,'stable');
        n=numel(commonP);
        if n<3, continue; end

        for im=1:numel(metrics)
            metric=metrics(im);
            x=A.(metric)(idxA);
            y=B.(metric)(idxB);
            delta=y-x; % B minus A
            [lo,hi]=bootstrapMedianCI(delta,cfg.nBootstrap);
            nz=delta(delta~=0 & isfinite(delta));
            if isempty(nz)
                consistency=NaN;
            else
                consistency=max(mean(nz>0),mean(nz<0));
            end
            rows(end+1,:)={dA,dB,metric,n,median(delta,'omitnan'), lo,hi,nnz(delta>0),nnz(delta<0),consistency,exactSignP(delta), strjoin(commonP,';')};
        end
    end
end

if isempty(rows)
    PairSummary=table(); return;
end

PairSummary=cell2table(rows,'VariableNames', {'drug_A','drug_B','metric','n_matched','median_B_minus_A', 'delta_CI_low','delta_CI_high','n_B_gt_A','n_B_lt_A', 'sign_consistency','exact_sign_p','patients'});
PairSummary.drug_A=string(PairSummary.drug_A);
PairSummary.drug_B=string(PairSummary.drug_B);
PairSummary.metric=string(PairSummary.metric);
PairSummary.patients=string(PairSummary.patients);

numvars={'n_matched','median_B_minus_A','delta_CI_low','delta_CI_high', 'n_B_gt_A','n_B_lt_A','sign_consistency','exact_sign_p'};
for i=1:numel(numvars)
    if iscell(PairSummary.(numvars{i}))
        PairSummary.(numvars{i})=cell2mat(PairSummary.(numvars{i}));
    end
end
end

function S = summarizeRobustness(Robust)
metrics=["coherence_C","contraction_K","magnitude_M"];
rows={};
for i=1:numel(metrics)
    m=metrics(i);
    x=Robust.(m+"_R");
    y=Robust.(m+"_Rconf");
    ok=isfinite(x) & isfinite(y);
    x=x(ok); y=y(ok);
    r=pearsonSimple(x,y);
    if isempty(x)
        signAgree=NaN; medAbs=NaN; medSigned=NaN;
    else
        signAgree=mean(sign(x)==sign(y));
        medAbs=median(abs(y-x),'omitnan');
        medSigned=median(y-x,'omitnan');
    end
    rows(end+1,:)={m,numel(x),r,signAgree,medAbs,medSigned};
end
S=cell2table(rows,'VariableNames', {'metric','n_pairs','pearson_r','sign_agreement', 'median_abs_change','median_signed_change'});
S.metric=string(S.metric);
numvars={'n_pairs','pearson_r','sign_agreement','median_abs_change','median_signed_change'};
for i=1:numel(numvars)
    if iscell(S.(numvars{i}))
        S.(numvars{i})=cell2mat(S.(numvars{i}));
    end
end
end

function S = summarizeMetricRelationships(R)
pairs = {
    'coherence_C','contraction_K','C vs K';
    'state_shift_TV','contraction_K','state shift vs K';
    'magnitude_M','contraction_K','magnitude vs K'
    };
rows={};

for i=1:size(pairs,1)
    xname=string(pairs{i,1});
    yname=string(pairs{i,2});
    label=string(pairs{i,3});
    x=R.(xname); y=R.(yname);
    rawR=pearsonSimple(x,y);

    xw=x; yw=y;
    pats=unique(R.patient);
    for ip=1:numel(pats)
        q=R.patient==pats(ip);
        xw(q)=x(q)-mean(x(q),'omitnan');
        yw(q)=y(q)-mean(y(q),'omitnan');
    end
    withinR=pearsonSimple(xw,yw);
    rows(end+1,:)={label,xname,yname,height(R),rawR,withinR};
end

S=cell2table(rows,'VariableNames', {'relationship','x_metric','y_metric','n_observations', 'raw_pearson_r','within_patient_centered_r'});
S.relationship=string(S.relationship);
S.x_metric=string(S.x_metric);
S.y_metric=string(S.y_metric);
numvars={'n_observations','raw_pearson_r','within_patient_centered_r'};
for i=1:numel(numvars)
    if iscell(S.(numvars{i}))
        S.(numvars{i})=cell2mat(S.(numvars{i}));
    end
end
end

function T = summarizeMagnitudeAdjustedK(R)

counts=groupsummary(table(R.patient),'Var1');
multiPatients=counts.Var1(counts.GroupCount>=2);
use=ismember(R.patient,multiPatients);
A=R(use,:);

Mw=A.magnitude_M;
Kw=A.contraction_K;
pats=unique(A.patient);
for ip=1:numel(pats)
    q=A.patient==pats(ip);
    Mw(q)=Mw(q)-mean(Mw(q),'omitnan');
    Kw(q)=Kw(q)-mean(Kw(q),'omitnan');
end

den=sum(Mw.^2,'omitnan');
if den<=0
    beta=NaN;
else
    beta=sum(Mw.*Kw,'omitnan')/den;
end
resid=Kw-beta.*Mw;

T=table(A.patient,A.drug,A.magnitude_M,A.contraction_K,Mw,Kw,resid, 'VariableNames',{'patient','drug','magnitude_M','contraction_K', 'magnitude_within_patient','K_within_patient','K_residual_after_M'});
T.Properties.UserData.beta=beta;
end

function r=pearsonSimple(x,y)
ok=isfinite(x) & isfinite(y);
x=x(ok); y=y(ok);
if numel(x)<2 || std(x)==0 || std(y)==0
    r=NaN; return;
end
C=corrcoef(x,y);
r=C(1,2);
end

function T = buildPatientAudit(patient,sequencingRaw,tissueRaw,treatment,malignant,analysisMask,R)
pats=unique(patient(patient~=""));
rows={};
for i=1:numel(pats)
    p=pats(i);
    q=patient==p;
    qm=q & malignant;
    qa=q & analysisMask;

    seq=unique(strtrim(string(sequencingRaw(q))));
    seq=seq(seq~="");
    pw=unique(normalizePatient(seq));
    pw=pw(pw~="");
    tis=unique(strtrim(string(tissueRaw(q))));
    tis=tis(tis~="");
    tr=unique(treatment(qm));
    tr=tr(tr~="");
    nonveh=tr(tr~="vehicle");

    rr=R(R.patient==p,:);
    rows(end+1,:)={p,nnz(q),nnz(qm),nnz(qm & treatment=="vehicle"), nnz(qa),numel(nonveh),height(rr),height(rr)>0, strjoin(pw,';'),strjoin(seq,';'),strjoin(tis,';'),strjoin(tr,';')};
end
T=cell2table(rows,'VariableNames', {'patient','n_cells_all','n_malignant','n_vehicle_malignant', 'n_analysis_cells','n_nonvehicle_treatments','n_geometry_drugs', 'contributes_geometry','pw_prefixes','sequencing_names','tissue_ids','treatments'});
T.patient=string(T.patient);
T.pw_prefixes=string(T.pw_prefixes);
T.sequencing_names=string(T.sequencing_names);
T.tissue_ids=string(T.tissue_ids);
T.treatments=string(T.treatments);
numvars={'n_cells_all','n_malignant','n_vehicle_malignant', 'n_analysis_cells','n_nonvehicle_treatments','n_geometry_drugs'};
for i=1:numel(numvars)
    if iscell(T.(numvars{i})), T.(numvars{i})=cell2mat(T.(numvars{i})); end
end
if iscell(T.contributes_geometry), T.contributes_geometry=cell2mat(T.contributes_geometry); end
end

function [hvgGenes,hvgIndex,Selection] = selectIndependentHVGs( file,meta,genes,lib,patient,treatment,analysisMask,stateMarkers,cfg)

[fileDir,base,~]=fileparts(file);
cacheFile=fullfile(fileDir,sprintf('%s_independent_HVG%d_cache.mat', base,cfg.nIndependentHVG));

if exist(cacheFile,'file')
    try
        C=load(cacheFile);
        if isfield(C,'cacheVersion') && strcmp(C.cacheVersion,cfg.hvgCacheVersion) && isfield(C,'hvgIndex') && numel(C.hvgIndex)==cfg.nIndependentHVG && isfield(C,'nGenes') && C.nGenes==meta.nGenes
            hvgIndex=C.hvgIndex;
            hvgGenes=C.hvgGenes;
            Selection=C.Selection;
            fprintf('Loaded cached independent %d-HVG panel: %s\n', cfg.nIndependentHVG,cacheFile);
            return;
        end
    catch
    end
end

fprintf('\nSelecting %d independent HVGs from VEHICLE malignant cells only\n', cfg.nIndependentHVG);
fprintf('State markers, mitochondrial and ribosomal genes are excluded.\n');

nG=meta.nGenes;
scoreCacheFile=fullfile(fileDir,sprintf('%s_independent_HVG_sampling_cache.mat',base));
loadedSamplingCache=false;

if exist(scoreCacheFile,'file')
    try
        Q=load(scoreCacheFile);
        if isfield(Q,'cacheVersion') && strcmp(Q.cacheVersion,cfg.hvgCacheVersion) && isfield(Q,'nGenes') && Q.nGenes==meta.nGenes && isfield(Q,'varSum') && isfield(Q,'varN') && isfield(Q,'detectTotal')
            varSum=Q.varSum;
            varN=Q.varN;
            detectTotal=Q.detectTotal;
            sampleTotal=Q.sampleTotal;
            sampledPatients=Q.sampledPatients;
            loadedSamplingCache=true;
            fprintf('Loaded cached vehicle-cell HVG sampling statistics: %s\n',scoreCacheFile);
        end
    catch
    end
end

if ~loadedSamplingCache
    varSum=zeros(1,nG);
    varN=zeros(1,nG);
    detectTotal=zeros(1,nG);
    sampleTotal=0;
    sampledPatients=0;

    allPats=unique(patient(analysisMask));
    hasDrug=false(size(allPats));
    for ii=1:numel(allPats)
        q=analysisMask & patient==allPats(ii);
        hasDrug(ii)=any(treatment(q)~="vehicle" & treatment(q)~="");
    end
    pats=allPats(hasDrug);

    for ip=1:numel(pats)
    p=pats(ip);
    idx=find(analysisMask & patient==p & treatment=="vehicle");
    if numel(idx)<cfg.minCellsPatientControl, continue; end

    breaks=[0; find(diff(idx)>1); numel(idx)];
    runStart=zeros(numel(breaks)-1,1);
    runLen=zeros(numel(breaks)-1,1);
    for r=1:numel(runLen)
        a=breaks(r)+1; b=breaks(r+1);
        runStart(r)=idx(a);
        runLen(r)=b-a+1;
    end
    [~,ord]=sort(runLen,'descend');
    ord=ord(1:min(cfg.hvgMaxRunsPerPatient,numel(ord)));

    sumP=zeros(1,nG);
    sumSqP=zeros(1,nG);
    detectP=zeros(1,nG);
    nP=0;

    for rr=1:numel(ord)
        r=ord(rr);
        take=min(cfg.hvgCellsPerRun,runLen(r));
        offset=floor((runLen(r)-take)/2);
        first=runStart(r)+offset;
        cellIdx=(first:first+take-1)';

        raw=readCellBlockAllGenes(file,meta,first,take);
        raw=single(raw);
        denom=single(lib(cellIdx));
        denom(denom<=0)=1;
        Y=log1p(1e4 .* raw ./ denom);

        sumP=sumP + sum(double(Y),1);
        sumSqP=sumSqP + sum(double(Y).^2,1);
        detectP=detectP + sum(raw>0,1);
        nP=nP+take;
        clear raw Y
    end

    if nP>=2
        mu=sumP/nP;
        v=(sumSqP - nP.*mu.^2)/max(1,nP-1);
        v(v<0 & v>-1e-10)=0;
        good=isfinite(v);
        varSum(good)=varSum(good)+v(good);
        varN(good)=varN(good)+1;
        detectTotal=detectTotal+detectP;
        sampleTotal=sampleTotal+nP;
        sampledPatients=sampledPatients+1;
    end
        fprintf('  HVG sampling: %s, %d vehicle cells (cumulative %d)\n', char(p),nP,sampleTotal);
    end

    cacheVersion=cfg.hvgCacheVersion;
    nGenes=meta.nGenes;
    save(scoreCacheFile,'cacheVersion','nGenes','varSum','varN','detectTotal', 'sampleTotal','sampledPatients');
    fprintf('Saved vehicle-cell HVG sampling cache: %s\n',scoreCacheFile);
end

score=(varSum ./ max(1,varN))';
detFrac=(detectTotal ./ max(1,sampleTotal))';
gUpper=geneListUpper(genes);
gUpper=gUpper(:);
markerUpper=upper(string(stateMarkers(:)));

exclude=ismember(gUpper,markerUpper) | startsWith(gUpper,"MT-") | startsWith(gUpper,"RPL") | startsWith(gUpper,"RPS") | gUpper=="MALAT1";

eligible=(~exclude) & (detFrac>=cfg.hvgMinDetection) & isfinite(score) & (score>0);

if numel(score)~=meta.nGenes || numel(detFrac)~=meta.nGenes || numel(eligible)~=meta.nGenes
    error('Internal HVG vector-shape error: expected %d genes.',meta.nGenes);
end

cand=find(eligible);
[~,ord]=sort(score(cand),'descend');
if numel(ord)<cfg.nIndependentHVG
    error('Only %d genes passed independent-HVG filters; %d requested.', numel(ord),cfg.nIndependentHVG);
end
hvgIndex=cand(ord(1:cfg.nIndependentHVG));
hvgGenes=string(genes(hvgIndex));
rank=(1:cfg.nIndependentHVG)';
selScore=score(hvgIndex);
selDetect=detFrac(hvgIndex);
Selection=table(rank,hvgGenes(:),selScore(:),selDetect(:), 'VariableNames',{'rank','gene','vehicle_within_patient_variance','detection_fraction'});

cacheVersion=cfg.hvgCacheVersion;
nGenes=meta.nGenes;
save(cacheFile,'cacheVersion','nGenes','hvgIndex','hvgGenes','Selection','sampleTotal','sampledPatients');
fprintf('Saved independent-HVG cache: %s\n',cacheFile);
fprintf('Selected %d HVGs from %d vehicle cells across %d TissueID units.\n\n', numel(hvgIndex),sampleTotal,sampledPatients);
end

function X=readCellBlockAllGenes(file,meta,firstCell,nCells)
if meta.genesOnRows
    X=h5read(file,meta.matrixPath,[1 firstCell],[meta.nGenes nCells]);
    X=X';
else
    X=h5read(file,meta.matrixPath,[firstCell 1],[nCells meta.nGenes]);
end
end

function v=readOneGeneFromLoom(file,meta,gidx)
if meta.genesOnRows
    v=h5read(file,meta.matrixPath,[gidx 1],[1 meta.nCells]);
else
    v=h5read(file,meta.matrixPath,[1 gidx],[meta.nCells 1]);
end
v=double(v(:));
end

function Rind = computeIndependentHVGGeometry( file,meta,hvgIndex,lib,patient,treatment,analysisMask,state,R,states,cfg)

[fileDir,base,~]=fileparts(file);
cacheFile=fullfile(fileDir,sprintf('%s_independent_HVG%d_geometry_cache.mat', base,numel(hvgIndex)));

if exist(cacheFile,'file')
    try
        C=load(cacheFile);
        if isfield(C,'cacheVersion') && strcmp(C.cacheVersion,cfg.hvgCacheVersion) && isfield(C,'hvgIndexSaved') && isequal(C.hvgIndexSaved(:),hvgIndex(:)) && isfield(C,'Rind')
            Rind=C.Rind;
            fprintf('Loaded cached independent-HVG geometry: %s\n',cacheFile);
            return;
        end
    catch
    end
end

fprintf('Computing C/K in independent %d-HVG expression space\n',numel(hvgIndex));
fprintf('State labels use the 34-marker assignments.\n');

pats=unique(R.patient,'stable');
nP=numel(pats);
nR=height(R);

ctrlAll=cell(nP,1);
ctrlState=cell(nP,4);
for ip=1:nP
    p=pats(ip);
    ctrlAll{ip}=find(analysisMask & patient==p & treatment=="vehicle");
    for k=1:4
        ctrlState{ip,k}=find(analysisMask & patient==p & treatment=="vehicle" & state==k);
    end
end

rowP=zeros(nR,1);
valid=false(nR,4);
trState=cell(nR,4);
for r=1:nR
    rowP(r)=find(pats==R.patient(r),1,'first');
    used=split(string(R.states_used(r)),',');
    for k=1:4
        valid(r,k)=any(used==string(states{k}));
        if valid(r,k)
            trState{r,k}=find(analysisMask & patient==R.patient(r) & treatment==R.drug(r) & state==k);
        end
    end
end

normSq=zeros(nR,4);
dotSum=zeros(nR,4,4);
d0Sq=zeros(nR,4,4);
d1Sq=zeros(nR,4,4);

for g=1:numel(hvgIndex)
    y=readOneGeneFromLoom(file,meta,hvgIndex(g));
    y=log1p(1e4 .* y ./ lib);

    mu=zeros(nP,1);
    sd=zeros(nP,1);
    for ip=1:nP
        z=y(ctrlAll{ip});
        mu(ip)=mean(z,'omitnan');
        sd(ip)=std(z,0,'omitnan');
        if ~isfinite(sd(ip)) || sd(ip)<1e-6, sd(ip)=1; end
    end

    for r=1:nR
        ip=rowP(r);
        used=find(valid(r,:));
        c0=nan(1,4); c1=nan(1,4);
        for kk=1:numel(used)
            k=used(kk);
            a=ctrlState{ip,k};
            b=trState{r,k};
            c0(k)=(mean(y(a),'omitnan')-mu(ip))/sd(ip);
            c1(k)=(mean(y(b),'omitnan')-mu(ip))/sd(ip);
            delta=c1(k)-c0(k);
            normSq(r,k)=normSq(r,k)+delta.^2;
        end
        for ii=1:numel(used)-1
            i=used(ii);
            di=c1(i)-c0(i);
            for jj=ii+1:numel(used)
                j=used(jj);
                dj=c1(j)-c0(j);
                dotSum(r,i,j)=dotSum(r,i,j)+di*dj;
                d0Sq(r,i,j)=d0Sq(r,i,j)+(c0(i)-c0(j)).^2;
                d1Sq(r,i,j)=d1Sq(r,i,j)+(c1(i)-c1(j)).^2;
            end
        end
    end

    if mod(g,25)==0 || g==numel(hvgIndex)
        fprintf('  independent-space geometry: %d / %d genes (%.1f%%)\n', g,numel(hvgIndex),100*g/numel(hvgIndex));
    end
end

C=zeros(nR,1); K=zeros(nR,1); M=zeros(nR,1);
D0=zeros(nR,1); D1=zeros(nR,1);
for r=1:nR
    used=find(valid(r,:));
    cosv=[]; d0v=[]; d1v=[];
    for ii=1:numel(used)-1
        i=used(ii);
        for jj=ii+1:numel(used)
            j=used(jj);
            den=sqrt(normSq(r,i)*normSq(r,j));
            if den>cfg.eps
                cosv(end+1)=dotSum(r,i,j)/den;
            end
            d0v(end+1)=sqrt(max(0,d0Sq(r,i,j)));
            d1v(end+1)=sqrt(max(0,d1Sq(r,i,j)));
        end
    end
    C(r)=mean(cosv,'omitnan');
    D0(r)=mean(d0v,'omitnan');
    D1(r)=mean(d1v,'omitnan');
    K(r)=log2((D1(r)+cfg.eps)/(D0(r)+cfg.eps));
    M(r)=mean(sqrt(normSq(r,used)),'omitnan');
end

Rind=R;
Rind.coherence_C=C;
Rind.contraction_K=K;
Rind.magnitude_M=M;
Rind.distance_control=D0;
Rind.distance_treated=D1;

cacheVersion=cfg.hvgCacheVersion;
hvgIndexSaved=hvgIndex;
save(cacheFile,'cacheVersion','hvgIndexSaved','Rind');
fprintf('Saved independent-HVG geometry cache: %s\n\n',cacheFile);
end

function S=comparePrimaryIndependentGeometry(R,Rind)
metrics=["coherence_C","contraction_K","magnitude_M"];
rows={};
for i=1:numel(metrics)
    m=metrics(i);
    x=R.(m); y=Rind.(m);
    ok=isfinite(x) & isfinite(y);
    x=x(ok); y=y(ok);
    if isempty(x)
        r=NaN; signAgree=NaN; medAbs=NaN; medSigned=NaN;
    else
        r=pearsonSimple(x,y);
        signAgree=mean(sign(x)==sign(y));
        medAbs=median(abs(y-x),'omitnan');
        medSigned=median(y-x,'omitnan');
    end
    rows(end+1,:)={m,numel(x),r,signAgree,medAbs,medSigned};
end
S=cell2table(rows,'VariableNames', {'metric','n_pairs','pearson_r','sign_agreement', 'median_abs_difference','median_signed_difference'});
S.metric=string(S.metric);
numvars=S.Properties.VariableNames(2:end);
for i=1:numel(numvars)
    if iscell(S.(numvars{i})), S.(numvars{i})=cell2mat(S.(numvars{i})); end
end
end

function M=getNeftelModules2019()
M.AC = { 'CST3','S100B','SLC1A3','HEPN1','HOPX','MT3','SPARCL1', 'MLC1','GFAP','FABP7','BCAN','PON2','METTL7B','SPARC', 'GATM','RAMP1','PMP2','AQP4','DBI','EDNRB','PTPRZ1', 'CLU','PMP22','ATP1A2','S100A16','HEY1','PCDHGC3','TTYH1', 'NDRG2','PRCP','ATP1B2','AGT','PLTP','GPM6B','F3', 'RAB31','PPAP2B','ANXA5','TSPAN7'
    };

M.OPC = { 'BCAN','PLP1','GPR17','FIBIN','LHFPL3','OLIG1','PSAT1', 'SCRG1','OMG','APOD','SIRT2','TNR','THY1','PHYHIPL', 'SOX2-OT','NKAIN4','LPPR1','PTPRZ1','VCAN','DBI','PMP2', 'CNP','TNS3','LIMA1','CA10','PCDHGC3','CNTN1','SCD5', 'P2RX7','CADM2','TTYH1','FGF12','TMEM206','NEU4','FXYD6', 'RNF13','RTKN','GPM6B','LMF1','ALCAM','PGRMC1','HRASLS', 'BCAS1','RAB31','PLLP','FABP5','NLGN3','SERINC5','EPB41L2', 'GPR37L1'
    };

M.NPC1 = { 'DLL3','DLL1','SOX4','TUBB3','HES6','TAGLN3','NEU4', 'MARCKSL1','CD24','STMN1','TCF12','BEX1','OLIG1','MAP2', 'FXYD6','PTPRS','MLLT11','NPPA','BCAN','MEST','ASCL1', 'BTG2','DCX','NXPH1','HN1','PFN2','SCG3','MYT1', 'CHD7','GPR56','TUBA1A','PCBP4','ETV1','SHD','TNR', 'AMOTL2','DBN1','HIP1','ABAT','ELAVL4','LMF1','GRIK2', 'SERINC5','TSPAN13','ELMO1','GLCCI1','SEZ6L','LRRN1','SEZ6', 'SOX11'
    };

M.NPC2 = { 'STMN2','CD24','RND3','HMP19','TUBB3','MIAT','DCX', 'NSG1','ELAVL4','MLLT11','DLX6-AS1','SOX11','NREP','FNBP1L', 'TAGLN3','STMN4','DLX5','SOX4','MAP1B','RBFOX2','IGFBPL1', 'STMN1','HN1','TMEM161B-AS1','DPYSL3','SEPT3','PKIA','ATP1B1', 'DYNC1I1','CD200','SNAP25','PAK3','NDRG4','KIF5A','UCHL1', 'ENO2','KIF5C','DDAH2','TUBB2A','LBH','LOC150568','TCF4', 'GNG3','NFIB','DPYSL5','CRABP1','DBN1','NFIX','CEP170', 'BLCAP'
    };

M.MES1 = { 'CHI3L1','ANXA2','ANXA1','CD44','VIM','MT2A','C1S', 'NAMPT','EFEMP1','C1R','SOD2','IFITM3','TIMP1','SPP1', 'A2M','S100A11','MT1X','S100A10','FN1','LGALS1','S100A16', 'CLIC1','MGST1','RCAN1','TAGLN2','NPC2','SERPING1','C8orf4', 'EMP1','APOE','CTSB','C3','LGALS3','MT1E','EMP3', 'SERPINA3','ACTN1','PRDX6','IGFBP7','SERPINE1','PLP2','MGP', 'CLIC4','GFPT2','GSN','NNMT','TUBA1C','GJA1','TNFRSF1A', 'WWTR1'
    };

M.MES2 = { 'HILPDA','ADM','DDIT3','NDRG1','HERPUD1','DNAJB9','TRIB3', 'ENO2','AKAP12','SQSTM1','MT1X','ATF3','NAMPT','NRN1', 'SLC2A1','BNIP3','LGALS3','INSIG2','IGFBP3','PPP1R15A','VIM', 'PLOD2','GBE1','SLC2A3','FTL','WARS','ERO1L','XPOT', 'HSPA5','GDF15','ANXA2','EPAS1','LDHA','P4HA1','SERTAD1', 'PFKP','PGK1','EGLN3','SLC6A6','CA9','BNIP3L','RPL21', 'TRAM1','UFM1','ASNS','GOLT1B','ANGPTL4','SLC39A14','CDKN1A', 'HSPA9'
    };
end

function [Zgene,geneNames,moduleCols,Coverage] = loadNeftelModuleExpression( file,meta,genes,lib,patient,treatment,analysisMask,M,cfg)

[fileDir,base,~]=fileparts(file);
cacheFile=fullfile(fileDir,[base '_neftel_module_expression_cache.mat']);
moduleNames={'AC','OPC','NPC1','NPC2','MES1','MES2'};

allGenes={};
for i=1:numel(moduleNames)
    allGenes=[allGenes M.(moduleNames{i})];
end
allGenes=unique(allGenes,'stable');

geneUpper=geneListUpper(genes);
[foundGenes,geneIndex,missingAll]=mapGenes(geneUpper,allGenes);
geneNames=string(foundGenes(:));
nA=nnz(analysisMask);
analysisIndex=find(analysisMask);

moduleCols=cell(numel(moduleNames),1);
covRows={};
for i=1:numel(moduleNames)
    published=string(M.(moduleNames{i}));
    present=published(ismember(upper(published),upper(geneNames)));
    missing=published(~ismember(upper(published),upper(geneNames)));
    moduleCols{i}=find(ismember(upper(geneNames),upper(present)));
    covRows(end+1,:)={string(moduleNames{i}),numel(published),numel(present), numel(present)/numel(published),strjoin(missing,';')};
end
Coverage=cell2table(covRows,'VariableNames', {'module','n_published','n_found','fraction_found','missing_genes'});
Coverage.module=string(Coverage.module);
Coverage.missing_genes=string(Coverage.missing_genes);
for v={'n_published','n_found','fraction_found'}
    if iscell(Coverage.(v{1})), Coverage.(v{1})=cell2mat(Coverage.(v{1})); end
end

if any(Coverage.fraction_found<0.80)
    warning('One or more Neftel modules have <80%% gene coverage in the loom.');
end
for i=1:numel(moduleNames)
    if numel(moduleCols{i})<20
        warning('%s has only %d detected genes.',moduleNames{i},numel(moduleCols{i}));
    end
end

if exist(cacheFile,'file')
    try
        C=load(cacheFile);
        if isfield(C,'cacheVersion') && strcmp(C.cacheVersion,cfg.neftelCacheVersion) && isfield(C,'geneNamesSaved') && isequal(C.geneNamesSaved(:),geneNames(:)) && isfield(C,'analysisIndexSaved') && isequal(C.analysisIndexSaved(:),analysisIndex(:))
            Zgene=C.Zgene;
            fprintf('Loaded cached Neftel module expression: %s\n',cacheFile);
            return;
        end
    catch
    end
end

fprintf('\nReading and standardizing %d unique Neftel module genes\n',numel(geneIndex));
fprintf('Reference for every gene is the patient-matched vehicle malignant population.\n');

patientA=patient(analysisIndex);
treatmentA=treatment(analysisIndex);
pats=unique(patientA);
Zgene=nan(nA,numel(geneIndex),'single');

for g=1:numel(geneIndex)
    raw=readOneGeneFromLoom(file,meta,geneIndex(g));
    y=log1p(1e4 .* raw ./ lib);
    yA=y(analysisIndex);
    z=nan(nA,1,'single');

    for ip=1:numel(pats)
        q=patientA==pats(ip);
        qc=q & treatmentA=="vehicle";
        mu=mean(yA(qc),'omitnan');
        sd=std(yA(qc),0,'omitnan');
        if ~isfinite(sd) || sd<1e-6, sd=1; end
        z(q)=single((yA(q)-mu)./sd);
    end
    Zgene(:,g)=z;

    if mod(g,25)==0 || g==numel(geneIndex)
        fprintf('  Neftel module genes: %d / %d (%.1f%%)\n', g,numel(geneIndex),100*g/numel(geneIndex));
    end
end

cacheVersion=cfg.neftelCacheVersion;
geneNamesSaved=geneNames;
analysisIndexSaved=analysisIndex;
try
    save(cacheFile,'cacheVersion','geneNamesSaved','analysisIndexSaved','Zgene','-v7.3');
catch
    save(cacheFile,'cacheVersion','geneNamesSaved','analysisIndexSaved','Zgene');
end
fprintf('Saved Neftel module-expression cache: %s\n\n',cacheFile);
end

function [state,margin,S4]=assignNeftelFourStates(Z,moduleCols,allowed)
if nargin<3 || isempty(allowed)
    allowed=true(size(Z,2),1);
end
allowed=logical(allowed(:));

scores=nan(size(Z,1),6);
for i=1:6
    cols=moduleCols{i};
    cols=cols(allowed(cols));
    if isempty(cols)
        error('No genes available for Neftel module %d in this split.',i);
    end
    scores(:,i)=mean(Z(:,cols),2,'omitnan');
end

S4=[scores(:,1),scores(:,2),max(scores(:,3:4),[],2),max(scores(:,5:6),[],2)];
[~,state]=max(S4,[],2);
tmp=sort(S4,2,'descend');
margin=tmp(:,1)-tmp(:,2);
end

function R=computeGeometryLocal(patient,treatment,state,Z,states,cfg)
rows={};
pats=unique(patient);
for ip=1:numel(pats)
    p=pats(ip);
    base=patient==p;
    ctrl=base & treatment=="vehicle";
    if nnz(ctrl)<cfg.minCellsPatientControl, continue; end
    drugs=unique(treatment(base & treatment~="vehicle" & treatment~=""));
    for id=1:numel(drugs)
        d=drugs(id);
        tr=base & treatment==d;
        if nnz(tr)<cfg.minCellsPatientControl/2, continue; end

        C0=[]; C1=[]; n0=[]; n1=[]; used=string.empty;
        p0=zeros(1,4); p1=zeros(1,4);
        for k=1:4
            a=ctrl & state==k;
            b=tr & state==k;
            if nnz(a)>=cfg.minCellsPerState && nnz(b)>=cfg.minCellsPerState
                C0(end+1,:)=mean(Z(a,:),1,'omitnan');
                C1(end+1,:)=mean(Z(b,:),1,'omitnan');
                n0(end+1)=nnz(a);
                n1(end+1)=nnz(b);
                used(end+1)=string(states{k});
            end
            p0(k)=nnz(a)/max(1,nnz(ctrl));
            p1(k)=nnz(b)/max(1,nnz(tr));
        end

        ns=size(C0,1);
        if ns<cfg.minStatesForMetric, continue; end
        Delta=C1-C0;
        coherence=meanPairwiseCosine(Delta,cfg.eps);
        D0=meanPairwiseDistance(C0);
        D1=meanPairwiseDistance(C1);
        K=log2((D1+cfg.eps)/(D0+cfg.eps));
        magnitude=mean(sqrt(sum(Delta.^2,2)),'omitnan');
        tv=0.5*sum(abs(p1-p0));

        rows(end+1,:)={p,d,nnz(ctrl),nnz(tr),ns,strjoin(used,','), coherence,K,magnitude,D0,D1,tv,min(n0),min(n1)};
    end
end

if isempty(rows)
    R=table(); return;
end
R=cell2table(rows,'VariableNames',{'patient','drug','n_control','n_treated', 'n_states','states_used','coherence_C','contraction_K','magnitude_M', 'distance_control','distance_treated','state_shift_TV', 'min_state_control','min_state_treated'});
R.patient=string(R.patient); R.drug=string(R.drug);
numvars=R.Properties.VariableNames(3:end);
for i=1:numel(numvars)
    if iscell(R.(numvars{i})) && numvars{i}~="states_used"
        R.(numvars{i})=cell2mat(R.(numvars{i}));
    end
end
end

function V=compareGeometryByKey(R,T)
metrics=["coherence_C","contraction_K","magnitude_M"];
keyR=R.patient+"|"+R.drug;
keyT=T.patient+"|"+T.drug;
[~,ia,ib]=intersect(keyR,keyT,'stable');
rows={};

for i=1:numel(metrics)
    m=metrics(i);
    x=R.(m)(ia); y=T.(m)(ib);
    ok=isfinite(x) & isfinite(y);
    x=x(ok); y=y(ok);
    rows(end+1,:)={m,numel(x),pearsonSimple(x,y), mean(sign(x)==sign(y)),median(abs(y-x),'omitnan'), median(y-x,'omitnan')};
end
V=cell2table(rows,'VariableNames',{'metric','n_pairs','pearson_r', 'sign_agreement','median_abs_difference','median_signed_difference'});
V.metric=string(V.metric);
for i=2:width(V)
    if iscell(V.(i)), V.(i)=cell2mat(V.(i)); end
end
end

function [Summary,Confusion]=compareStateAssignments(a,b,patient,states)
ok=a>=1 & a<=4 & b>=1 & b<=4;
a=a(ok); b=b(ok); patient=patient(ok);
po=mean(a==b);
pa=zeros(4,1); pb=zeros(4,1);
for k=1:4
    pa(k)=mean(a==k);
    pb(k)=mean(b==k);
end
pe=sum(pa.*pb);
if pe<1
    kappa=(po-pe)/(1-pe);
else
    kappa=NaN;
end

pats=unique(patient);
patientAgreement=nan(numel(pats),1);
for i=1:numel(pats)
    q=patient==pats(i);
    patientAgreement(i)=mean(a(q)==b(q));
end

Summary=table(nnz(ok),po,kappa,median(patientAgreement,'omitnan'), min(patientAgreement),max(patientAgreement), 'VariableNames',{'n_cells','overall_agreement','cohen_kappa', 'median_patient_agreement','min_patient_agreement','max_patient_agreement'});

rows={};
for i=1:4
    denom=max(1,nnz(a==i));
    for j=1:4
        n=nnz(a==i & b==j);
        rows(end+1,:)={string(states{i}),string(states{j}),n,n/denom};
    end
end
Confusion=cell2table(rows,'VariableNames', {'compact_state','neftel_state','n_cells','fraction_within_compact_state'});
Confusion.compact_state=string(Confusion.compact_state);
Confusion.neftel_state=string(Confusion.neftel_state);
if iscell(Confusion.n_cells), Confusion.n_cells=cell2mat(Confusion.n_cells); end
if iscell(Confusion.fraction_within_compact_state)
    Confusion.fraction_within_compact_state=cell2mat(Confusion.fraction_within_compact_state);
end
end

function [CrossAll,Summary,Validation,Assignment] = runNeftelCrossfit( Z,moduleCols,patient,treatment,compactState,Rprimary,states,cfg)

nGenes=size(Z,2);
nExpected=2*cfg.neftelCrossfitSplits;
allTables=cell(nExpected,1);
assignRows={};
tt=0;

for s=1:cfg.neftelCrossfitSplits
    [A,B,seedUsed]=balancedNeftelSplit(nGenes,moduleCols, cfg.neftelCrossfitSeed+s-1,cfg.neftelMinGenesPerHalf);

    [stateA,~]=assignNeftelFourStates(Z,moduleCols,A);
    [stateB,~]=assignNeftelFourStates(Z,moduleCols,B);

    [agr,kappa]=stateAgreementSimple(stateA,stateB);
    agrCA=mean(stateA==compactState);
    agrCB=mean(stateB==compactState);

    countsA=zeros(1,6); countsB=zeros(1,6);
    for m=1:6
        countsA(m)=nnz(A(moduleCols{m}));
        countsB(m)=nnz(B(moduleCols{m}));
    end
    assignRows(end+1,:)={s,seedUsed,agr,kappa,agrCA,agrCB, min(countsA),max(countsA),min(countsB),max(countsB)};

    T=computeGeometryLocal(patient,treatment,stateA,Z(:,B),states,cfg);
    tt=tt+1;
    if ~isempty(T)
        T.split_id=repmat(s,height(T),1);
        T.fold=repmat("A_assign_B_geometry",height(T),1);
        allTables{tt}=T;
    end

    T=computeGeometryLocal(patient,treatment,stateB,Z(:,A),states,cfg);
    tt=tt+1;
    if ~isempty(T)
        T.split_id=repmat(s,height(T),1);
        T.fold=repmat("B_assign_A_geometry",height(T),1);
        allTables{tt}=T;
    end

    fprintf('  Neftel cross-fit split %d/%d complete; A/B state agreement %.3f\n', s,cfg.neftelCrossfitSplits,agr);
end

allTables=allTables(~cellfun(@isempty,allTables));
if isempty(allTables)
    CrossAll=table(); Summary=table(); Validation=table(); Assignment=table();
    return;
end
CrossAll=vertcat(allTables{:});

Assignment=cell2table(assignRows,'VariableNames', {'split_id','seed_used','agreement_A_B','kappa_A_B', 'agreement_compact_A','agreement_compact_B', 'min_module_genes_A','max_module_genes_A', 'min_module_genes_B','max_module_genes_B'});
for i=1:width(Assignment)
    if iscell(Assignment.(i)), Assignment.(i)=cell2mat(Assignment.(i)); end
end

keyP=Rprimary.patient+"|"+Rprimary.drug;
rows={};
for i=1:height(Rprimary)
    p=Rprimary.patient(i); d=Rprimary.drug(i);
    q=CrossAll.patient==p & CrossAll.drug==d;
    C=CrossAll.coherence_C(q);
    K=CrossAll.contraction_K(q);
    M=CrossAll.magnitude_M(q);
    n=nnz(q);

    [cq25,cq75]=quartilesLocal(C);
    [kq25,kq75]=quartilesLocal(K);
    [mq25,mq75]=quartilesLocal(M);

    if n>0
        signAgree=mean(sign(K)==sign(Rprimary.contraction_K(i)));
        positiveFraction=mean(K>0);
    else
        signAgree=NaN; positiveFraction=NaN;
    end
    eligible=n>=ceil(cfg.neftelMinValidFoldFraction*nExpected);

    rows(end+1,:)={p,d,n,nExpected,eligible, median(C,'omitnan'),cq25,cq75, median(K,'omitnan'),kq25,kq75,signAgree,positiveFraction, median(M,'omitnan'),mq25,mq75};
end

Summary=cell2table(rows,'VariableNames', {'patient','drug','n_valid_folds','n_expected_folds','eligible_validation', 'coherence_C','C_q25','C_q75','contraction_K','K_q25','K_q75', 'K_sign_agreement_primary','K_positive_fold_fraction', 'magnitude_M','M_q25','M_q75'});
Summary.patient=string(Summary.patient); Summary.drug=string(Summary.drug);
for i=3:width(Summary)
    if iscell(Summary.(i)), Summary.(i)=cell2mat(Summary.(i)); end
end

Validation=compareCrossfitSummary(Rprimary,Summary);
end

function [A,B,seedUsed]=balancedNeftelSplit(nGenes,moduleCols,seed0,minPerModule)
saved=rng;
cleanup=onCleanup(@()rng(saved));
for attempt=0:999
    seedUsed=seed0+1000*attempt;
    rng(seedUsed,'twister');
    perm=randperm(nGenes);
    A=false(nGenes,1);
    A(perm(1:floor(nGenes/2)))=true;
    B=~A;

    ok=true;
    for m=1:numel(moduleCols)
        if nnz(A(moduleCols{m}))<minPerModule || nnz(B(moduleCols{m}))<minPerModule
            ok=false; break;
        end
    end
    if ok, return; end
end
error('Could not generate a balanced globally disjoint Neftel gene split.');
end

function [agreement,kappa]=stateAgreementSimple(a,b)
ok=isfinite(a) & isfinite(b);
a=a(ok); b=b(ok);
agreement=mean(a==b);
pa=zeros(4,1); pb=zeros(4,1);
for k=1:4
    pa(k)=mean(a==k); pb(k)=mean(b==k);
end
pe=sum(pa.*pb);
if pe<1, kappa=(agreement-pe)/(1-pe); else, kappa=NaN; end
end

function [q25,q75]=quartilesLocal(x)
x=sort(x(isfinite(x)));
if isempty(x)
    q25=NaN; q75=NaN; return;
end
q25=interpQuantile(x,0.25);
q75=interpQuantile(x,0.75);
end

function q=interpQuantile(x,p)
n=numel(x);
if n==1, q=x(1); return; end
pos=1+(n-1)*p;
lo=floor(pos); hi=ceil(pos);
if lo==hi
    q=x(lo);
else
    q=x(lo)+(pos-lo)*(x(hi)-x(lo));
end
end

function V=compareCrossfitSummary(R,S)
S=S(S.eligible_validation,:);
metrics=["coherence_C","contraction_K","magnitude_M"];
keyR=R.patient+"|"+R.drug;
keyS=S.patient+"|"+S.drug;
[~,ia,ib]=intersect(keyR,keyS,'stable');
rows={};
for i=1:numel(metrics)
    m=metrics(i);
    x=R.(m)(ia); y=S.(m)(ib);
    ok=isfinite(x) & isfinite(y);
    x=x(ok); y=y(ok);
    rows(end+1,:)={m,numel(x),pearsonSimple(x,y), mean(sign(x)==sign(y)),median(abs(y-x),'omitnan'), median(y-x,'omitnan')};
end
V=cell2table(rows,'VariableNames',{'metric','n_pairs','pearson_r', 'sign_agreement','median_abs_difference','median_signed_difference'});
V.metric=string(V.metric);
for i=2:width(V)
    if iscell(V.(i)), V.(i)=cell2mat(V.(i)); end
end
end

function writeNeftelModules(path,M)
moduleNames={'AC','OPC','NPC1','NPC2','MES1','MES2'};
rows={};
for i=1:numel(moduleNames)
    g=string(M.(moduleNames{i}));
    for j=1:numel(g)
        rows(end+1,:)={string(moduleNames{i}),j,g(j), "Neftel et al., Cell 2019, Table S2", "10.1016/j.cell.2019.06.024"};
    end
end
T=cell2table(rows,'VariableNames', {'module','rank_in_published_module','gene','source','doi'});
T.module=string(T.module); T.gene=string(T.gene);
T.source=string(T.source); T.doi=string(T.doi);
if iscell(T.rank_in_published_module)
    T.rank_in_published_module=cell2mat(T.rank_in_published_module);
end
writetable(T,path);
end

function s=prettyDrug(s)
s=string(s);
s=replace(s,"rsl3_ferrostatin","RSL3 + ferrostatin");
s=replace(s,"rsl3","RSL3");
s=replace(s,"panobinostat","panobinostat");
s=replace(s,"etoposide","etoposide");
s=replace(s,"ispinesib","ispinesib");
s=replace(s,"givinostat","givinostat");
s=replace(s,"topotecan","topotecan");
end

function writeKeyResults(path,R,Disc,S,PairSummary,RobustSummary,RelationshipSummary,MagnitudeAdjusted,Rind,IndSummary,IndValidation,PatientAudit,NeftelFullValidation,NeftelAssignmentSummary,CrossfitSummary,CrossfitValidation,CrossfitAssignment,NullSplit,NullPerm,NullSummary,CrossDrug,CrossDrugSummary,RotNull,RotNullSummary,DrugStratifiedSign,discoveryPatient)
fid=fopen(path,'w');
fprintf(fid,'GBM DRUG-STATE GEOMETRY: KEY RESULTS\n');
fprintf(fid,'====================================\n\n');
fprintf(fid,'Patient-drug observations retained: %d\n',height(R));
fprintf(fid,'Patients contributing retained drug comparisons: %d\n',numel(unique(R.patient)));
fprintf(fid,'Replicated drugs (>=3 TissueID units): %d\n\n',height(S));

if any(Disc.drug=="panobinostat")
    x=Disc(Disc.drug=="panobinostat",:);
    fprintf(fid,'PW030 panobinostat: C = %.4f, K = %.4f, M = %.4f, distance ratio = %.3f\n', x.coherence_C,x.contraction_K,x.magnitude_M,2.^x.contraction_K);
end

if any(S.drug=="ANA-12")
    x=S(S.drug=="ANA-12",:);
    fprintf(fid,'ANA-12 across %d patients: median K = %.4f, bootstrap CI [%.4f, %.4f], sign consistency = %.3f, sign P = %.4f\n', x.n_patients,x.median_K,x.K_CI_low,x.K_CI_high,x.K_sign_consistency,x.K_exact_sign_p);
end
if any(S.drug=="panobinostat")
    x=S(S.drug=="panobinostat",:);
    fprintf(fid,'Panobinostat across %d patients: median K = %.4f, bootstrap CI [%.4f, %.4f], sign consistency = %.3f, sign P = %.4f\n', x.n_patients,x.median_K,x.K_CI_low,x.K_CI_high,x.K_sign_consistency,x.K_exact_sign_p);
end
if any(S.drug=="rsl3")
    x=S(S.drug=="rsl3",:);
    fprintf(fid,'RSL3 across %d patients: median C = %.4f, median K = %.4f, K sign consistency = %.3f\n', x.n_patients,x.median_C,x.median_K,x.K_sign_consistency);
end

fprintf(fid,'\nNull calibration of C, K and M:\n');
for i=1:height(NullSummary)
    fprintf(fid,'  %s, %s: n=%d, median observed=%.4f, median null=%.4f, median excess=%+.4f, %d/%d with empirical P<0.05\n', char(NullSummary.null_model(i)),char(NullSummary.metric(i)), NullSummary.n_comparisons(i),NullSummary.median_observed(i), NullSummary.median_null(i),NullSummary.median_excess(i), NullSummary.n_significant_0p05(i),NullSummary.n_comparisons(i));
end
q=NullSplit.patient==discoveryPatient & NullSplit.drug=="panobinostat";
if any(q)
    fprintf(fid,'  PW030 panobinostat vehicle-split null: C=%.4f versus null median %.4f [%.4f, %.4f], P=%.4f\n', NullSplit.coherence_C(q),NullSplit.null_C_median(q), NullSplit.null_C_lo(q),NullSplit.null_C_hi(q),NullSplit.p_C(q));
end
fprintf(fid,'  Vehicle-split draws reuse the unit vehicle pool; label-permutation draws reuse the observed cells with condition labels shuffled within state.\n');

if ~isempty(CrossDrugSummary)
    fprintf(fid,'\nWithin-unit cross-drug coherence:\n');
    fprintf(fid,'  %d drug pairs in %d units; median within-drug C = %.4f, median between-drug C (different states) = %.4f, median difference = %+.4f\n', CrossDrugSummary.n_drug_pairs,CrossDrugSummary.n_units, CrossDrugSummary.median_within_drug_C,CrossDrugSummary.median_between_drug_C, CrossDrugSummary.median_within_minus_between);
    fprintf(fid,'  Median between-drug coherence for the same state across two drugs = %.4f\n', CrossDrugSummary.median_between_drug_same_state_C);
    fprintf(fid,'  Between-drug coherence close to within-drug coherence would indicate a shared response axis that is not drug specific.\n');
end

fprintf(fid,'\nRotation null (independent random directions per state, observed norms fixed):\n');
fprintf(fid,'  n=%d, median observed C=%.4f, median null C=%.4f, median excess=%+.4f, %d/%d with empirical P<0.05\n', RotNullSummary.n_comparisons,RotNullSummary.median_observed,RotNullSummary.median_null, RotNullSummary.median_excess,RotNullSummary.n_significant_0p05,RotNullSummary.n_comparisons);

fprintf(fid,'\nDrug-stratified sign test for K (one vote per drug, from its median-unit sign):\n');
fprintf(fid,'  %d of %d drugs vote positive\n', nnz(DrugStratifiedSign.vote>0),nnz(DrugStratifiedSign.vote~=0));

fprintf(fid,'\nExploratory matched-patient comparisons (K):\n');
targets = {
    "etoposide","panobinostat";
    "ispinesib","panobinostat"
    };
for i=1:2
    q=PairSummary.metric=="contraction_K" & PairSummary.drug_A==targets{i,1} & PairSummary.drug_B==targets{i,2};
    if any(q)
        x=PairSummary(q,:);
        fprintf(fid,'  %s minus %s: n=%d, median delta K=%.4f, %d/%d positive, exact sign P=%.4f\n', char(prettyDrug(x.drug_B)),char(prettyDrug(x.drug_A)), x.n_matched,x.median_B_minus_A,x.n_B_gt_A,x.n_matched,x.exact_sign_p);
    end
end

fprintf(fid,'\nRobustness to high-confidence state assignment:\n');
for i=1:height(RobustSummary)
    fprintf(fid,'  %s: n=%d, Pearson r=%.4f, sign agreement=%.3f, median absolute change=%.4f, median signed change=%+.4f\n', char(RobustSummary.metric(i)),RobustSummary.n_pairs(i), RobustSummary.pearson_r(i),RobustSummary.sign_agreement(i), RobustSummary.median_abs_change(i),RobustSummary.median_signed_change(i));
end

fprintf(fid,'\nDescriptive metric relationships:\n');
for i=1:height(RelationshipSummary)
    fprintf(fid,'  %s: raw r=%.4f, patient-centered r=%.4f\n', char(RelationshipSummary.relationship(i)), RelationshipSummary.raw_pearson_r(i), RelationshipSummary.within_patient_centered_r(i));
end

fprintf(fid,'\nMagnitude-adjusted K sensitivity (patients with >=2 drugs):\n');
fprintf(fid,'  within-patient K~M slope = %.4f\n',MagnitudeAdjusted.Properties.UserData.beta);
if any(MagnitudeAdjusted.patient==discoveryPatient & MagnitudeAdjusted.drug=="panobinostat")
    q=MagnitudeAdjusted.patient==discoveryPatient & MagnitudeAdjusted.drug=="panobinostat";
    fprintf(fid,'  PW030 panobinostat residual K after M adjustment = %.4f\n', MagnitudeAdjusted.K_residual_after_M(q));
end

fprintf(fid,'\nIndependent 500-HVG feature-space validation:\n');
for i=1:height(IndValidation)
    fprintf(fid,'  %s: n=%d, Pearson r=%.4f, sign agreement=%.3f, median absolute difference=%.4f\n', char(IndValidation.metric(i)),IndValidation.n_pairs(i), IndValidation.pearson_r(i),IndValidation.sign_agreement(i), IndValidation.median_abs_difference(i));
end
if any(Rind.patient==discoveryPatient & Rind.drug=="panobinostat")
    q=Rind.patient==discoveryPatient & Rind.drug=="panobinostat";
    fprintf(fid,'  PW030 panobinostat independent-space: C=%.4f, K=%.4f, M=%.4f\n', Rind.coherence_C(q),Rind.contraction_K(q),Rind.magnitude_M(q));
end

fprintf(fid,'\nPublished Neftel full-module sensitivity:\n');
for i=1:height(NeftelFullValidation)
    fprintf(fid,'  %s: n=%d, Pearson r=%.4f, sign agreement=%.3f, median absolute difference=%.4f\n', char(NeftelFullValidation.metric(i)),NeftelFullValidation.n_pairs(i), NeftelFullValidation.pearson_r(i),NeftelFullValidation.sign_agreement(i), NeftelFullValidation.median_abs_difference(i));
end
fprintf(fid,'  Compact versus full-module state assignment agreement = %.4f; kappa = %.4f\n', NeftelAssignmentSummary.overall_agreement(1),NeftelAssignmentSummary.cohen_kappa(1));

fprintf(fid,'\nCross-fitted held-out Neftel state-program geometry:\n');
for i=1:height(CrossfitValidation)
    fprintf(fid,'  %s: n=%d, Pearson r=%.4f, sign agreement=%.3f, median absolute difference=%.4f\n', char(CrossfitValidation.metric(i)),CrossfitValidation.n_pairs(i), CrossfitValidation.pearson_r(i),CrossfitValidation.sign_agreement(i), CrossfitValidation.median_abs_difference(i));
end
fprintf(fid,'  Median A/B state-assignment agreement = %.4f; median kappa = %.4f\n', median(CrossfitAssignment.agreement_A_B,'omitnan'), median(CrossfitAssignment.kappa_A_B,'omitnan'));
q=CrossfitSummary.patient==discoveryPatient & CrossfitSummary.drug=="panobinostat";
if any(q)
    fprintf(fid,'  PW030 panobinostat cross-fitted median: C=%.4f, K=%.4f; K sign agreement with primary across folds=%.3f\n', CrossfitSummary.coherence_C(q),CrossfitSummary.contraction_K(q), CrossfitSummary.K_sign_agreement_primary(q));
end

fprintf(fid,'\nCohort audit:\n');
fprintf(fid,'  Biological TissueID units in loom: %d\n',height(PatientAudit));
fprintf(fid,'  TissueID units contributing >=1 drug geometry comparison: %d\n', nnz(PatientAudit.contributes_geometry));
fprintf(fid,'  See patient_cohort_audit.csv for TissueID-to-PW/sample mapping and treatment coverage.\n');

fprintf(fid,'\nInterpretation notes:\n');
fprintf(fid,'  Matched-pair and magnitude-adjusted analyses are exploratory sensitivity analyses.\n');
fprintf(fid,'  The 500-HVG analysis tests broad transcriptomic feature-space dependence.\n');
fprintf(fid,'  The Neftel cross-fit keeps state-assignment genes and geometry genes non-overlapping within each fold.\n');
fclose(fid);
end

function writeMarkerFile(path,marker,missing)
fid=fopen(path,'w');
fields=fieldnames(marker);
for i=1:numel(fields)
    fprintf(fid,'%s: %s\n',fields{i},strjoin(marker.(fields{i}),', '));
end
fprintf(fid,'\nMissing from input: %s\n',strjoin(missing,', '));
fclose(fid);
end

function writeDetectionFile(path,meta)
fid=fopen(path,'w');
fprintf(fid,'matrix: %s\n',meta.matrixPath);
fprintf(fid,'gene: %s\n',meta.genePath);
fprintf(fid,'patient: %s\n',meta.patientPath);
fprintf(fid,'treatment: %s\n',meta.treatmentPath);
fprintf(fid,'celltype: %s\n',meta.celltypePath);
fprintf(fid,'library_size: %s\n',meta.librarySizePath);
fprintf(fid,'nGenes: %d\n',meta.nGenes);
fprintf(fid,'nCells: %d\n',meta.nCells);
fprintf(fid,'genesOnRows: %d\n',meta.genesOnRows);
fclose(fid);
end

function Null = computeVehicleSplitNull(patient,treatment,mask,state,Z,R,states,cfg)
rows={};
for i=1:height(R)
    p=R.patient(i); d=R.drug(i);
    ctrl=mask & patient==p & treatment=="vehicle";
    tr=mask & patient==p & treatment==d;
    used=strsplit(R.states_used{i},',');
    Zk={}; aSize=[]; bSize=[]; Ssum=[];
    for u=1:numel(used)
        k=find(strcmp(states,strtrim(used{u})));
        if isempty(k), continue; end
        v=find(ctrl & state==k);
        nCtrl=numel(v); nTr=nnz(tr & state==k);
        b=round(nCtrl*nTr/max(1,nCtrl+nTr));
        a=nCtrl-b;
        if a>=cfg.minCellsPerState && b>=cfg.minCellsPerState
            Zk{end+1}=Z(v,:);
            aSize(end+1)=a; bSize(end+1)=b;
            Ssum(end+1,:)=sum(Z(v,:),1,'omitnan');
        end
    end
    ns=numel(Zk);
    if ns<cfg.minStatesForMetric
        rows(end+1,:)=nullRowStats(p,d,ns,R.coherence_C(i),R.contraction_K(i), R.magnitude_M(i),[],[],[],cfg.nNullDraws);
        continue;
    end
    cNull=nan(cfg.nNullDraws,1); kNull=nan(cfg.nNullDraws,1); mNull=nan(cfg.nNullDraws,1);
    for it=1:cfg.nNullDraws
        C0=zeros(ns,size(Z,2)); C1=zeros(ns,size(Z,2));
        for u=1:ns
            o=randperm(size(Zk{u},1));
            sb=sum(Zk{u}(o(1:bSize(u)),:),1,'omitnan');
            C1(u,:)=sb/bSize(u);
            C0(u,:)=(Ssum(u,:)-sb)/aSize(u);
        end
        Delta=C1-C0;
        cNull(it)=meanPairwiseCosine(Delta,cfg.eps);
        kNull(it)=log2((meanPairwiseDistance(C1)+cfg.eps)/(meanPairwiseDistance(C0)+cfg.eps));
        mNull(it)=mean(sqrt(sum(Delta.^2,2)),'omitnan');
    end
    rows(end+1,:)=nullRowStats(p,d,ns,R.coherence_C(i),R.contraction_K(i), R.magnitude_M(i),cNull,kNull,mNull,cfg.nNullDraws);
    fprintf('  vehicle-split null: %s / %s done\n',char(p),char(d));
end
Null=finishNullTable(rows);
end

function Null = computeLabelPermutationNull(patient,treatment,mask,state,Z,R,states,cfg)
rows={};
for i=1:height(R)
    p=R.patient(i); d=R.drug(i);
    ctrl=mask & patient==p & treatment=="vehicle";
    tr=mask & patient==p & treatment==d;
    used=strsplit(R.states_used{i},',');
    Zk={}; aSize=[]; bSize=[]; Ssum=[];
    for u=1:numel(used)
        k=find(strcmp(states,strtrim(used{u})));
        if isempty(k), continue; end
        a=find(ctrl & state==k); b=find(tr & state==k);
        if numel(a)>=cfg.minCellsPerState && numel(b)>=cfg.minCellsPerState
            pool=[a;b];
            Zk{end+1}=Z(pool,:);
            aSize(end+1)=numel(a); bSize(end+1)=numel(b);
            Ssum(end+1,:)=sum(Z(pool,:),1,'omitnan');
        end
    end
    ns=numel(Zk);
    if ns<cfg.minStatesForMetric
        rows(end+1,:)=nullRowStats(p,d,ns,R.coherence_C(i),R.contraction_K(i), R.magnitude_M(i),[],[],[],cfg.nNullDraws);
        continue;
    end
    cNull=nan(cfg.nNullDraws,1); kNull=nan(cfg.nNullDraws,1); mNull=nan(cfg.nNullDraws,1);
    for it=1:cfg.nNullDraws
        C0=zeros(ns,size(Z,2)); C1=zeros(ns,size(Z,2));
        for u=1:ns
            o=randperm(size(Zk{u},1));
            sb=sum(Zk{u}(o(aSize(u)+1:end),:),1,'omitnan');
            C1(u,:)=sb/bSize(u);
            C0(u,:)=(Ssum(u,:)-sb)/aSize(u);
        end
        Delta=C1-C0;
        cNull(it)=meanPairwiseCosine(Delta,cfg.eps);
        kNull(it)=log2((meanPairwiseDistance(C1)+cfg.eps)/(meanPairwiseDistance(C0)+cfg.eps));
        mNull(it)=mean(sqrt(sum(Delta.^2,2)),'omitnan');
    end
    rows(end+1,:)=nullRowStats(p,d,ns,R.coherence_C(i),R.contraction_K(i), R.magnitude_M(i),cNull,kNull,mNull,cfg.nNullDraws);
    fprintf('  label-permutation null: %s / %s done\n',char(p),char(d));
end
Null=finishNullTable(rows);
end

function out = nullRowStats(p,d,ns,cObs,kObs,mObs,cNull,kNull,mNull,B)
if isempty(cNull)
    out={p,d,ns,cObs,NaN,NaN,NaN,NaN,kObs,NaN,NaN,NaN,NaN,mObs,NaN,NaN,NaN,NaN};
    return;
end
cq=prctile(cNull,[2.5 50 97.5]);
kq=prctile(kNull,[2.5 50 97.5]);
mq=prctile(mNull,[2.5 50 97.5]);
pC=(1+nnz(cNull>=cObs))/(B+1);
pK=(1+nnz(abs(kNull)>=abs(kObs)))/(B+1);
pM=(1+nnz(mNull>=mObs))/(B+1);
out={p,d,ns,cObs,cq(2),cq(1),cq(3),pC,kObs,kq(2),kq(1),kq(3),pK, mObs,mq(2),mq(1),mq(3),pM};
end

function Null = finishNullTable(rows)
if isempty(rows)
    Null=table(); return;
end
Null=cell2table(rows,'VariableNames',{'patient','drug','n_states', 'coherence_C','null_C_median','null_C_lo','null_C_hi','p_C', 'contraction_K','null_K_median','null_K_lo','null_K_hi','p_K', 'magnitude_M','null_M_median','null_M_lo','null_M_hi','p_M'});
Null.patient=string(Null.patient); Null.drug=string(Null.drug);
numvars=Null.Properties.VariableNames(3:end);
for i=1:numel(numvars)
    if iscell(Null.(numvars{i})), Null.(numvars{i})=cell2mat(Null.(numvars{i})); end
end
end

function S = summarizeNullCalibration(NullSplit,NullPerm)
sets={'vehicle split',NullSplit;'label permutation',NullPerm};
metrics={'coherence_C','null_C_median','p_C';
    'contraction_K','null_K_median','p_K';
    'magnitude_M','null_M_median','p_M'};
rows={};
for s=1:size(sets,1)
    T=sets{s,2};
    for m=1:size(metrics,1)
        obs=T.(metrics{m,1}); nul=T.(metrics{m,2}); pv=T.(metrics{m,3});
        ok=isfinite(pv) & isfinite(obs) & isfinite(nul);
        if nnz(ok)==0
            rows(end+1,:)={sets{s,1},metrics{m,1},0,NaN,NaN,NaN,0,NaN};
            continue;
        end
        rows(end+1,:)={sets{s,1},metrics{m,1},nnz(ok), median(obs(ok)),median(nul(ok)),median(obs(ok)-nul(ok)), nnz(pv(ok)<0.05),mean(pv(ok)<0.05)};
    end
end
S=cell2table(rows,'VariableNames',{'null_model','metric','n_comparisons', 'median_observed','median_null','median_excess','n_significant_0p05', 'fraction_significant_0p05'});
S.null_model=string(S.null_model); S.metric=string(S.metric);
numvars=S.Properties.VariableNames(3:end);
for i=1:numel(numvars)
    if iscell(S.(numvars{i})), S.(numvars{i})=cell2mat(S.(numvars{i})); end
end
end

function X = computeCrossDrugCoherence(patient,treatment,mask,state,Z,R,states,cfg)
rows={};
pats=unique(R.patient);
for ip=1:numel(pats)
    p=pats(ip);
    A=R(R.patient==p,:);
    if height(A)<2, continue; end
    ctrl=mask & patient==p & treatment=="vehicle";
    Dl={}; lab=string.empty; st={};
    for id=1:height(A)
        d=A.drug(id);
        tr=mask & patient==p & treatment==d;
        used=strsplit(A.states_used{id},',');
        Dd=[]; su=string.empty;
        for u=1:numel(used)
            k=find(strcmp(states,strtrim(used{u})));
            if isempty(k), continue; end
            a=ctrl & state==k; b=tr & state==k;
            Dd(end+1,:)=mean(Z(b,:),1,'omitnan')-mean(Z(a,:),1,'omitnan');
            su(end+1)=string(states{k});
        end
        Dl{end+1}=Dd; lab(end+1)=d; st{end+1}=su;
    end
    for i=1:numel(Dl)-1
        for j=i+1:numel(Dl)
            [shared,ia,ib]=intersect(st{i},st{j},'stable');
            if numel(shared)<cfg.minStatesForMetric, continue; end
            Di=Dl{i}(ia,:); Dj=Dl{j}(ib,:);
            diffState=[]; sameState=[];
            for a=1:numel(shared)
                for b=1:numel(shared)
                    den=norm(Di(a,:))*norm(Dj(b,:));
                    if den<=cfg.eps, continue; end
                    c=dot(Di(a,:),Dj(b,:))/den;
                    if a==b, sameState(end+1)=c; else, diffState(end+1)=c; end
                end
            end
            if isempty(diffState), continue; end
            rows(end+1,:)={p,lab(i),lab(j),numel(shared), meanPairwiseCosine(Di,cfg.eps),meanPairwiseCosine(Dj,cfg.eps), mean(diffState),mean(sameState)};
        end
    end
end
if isempty(rows)
    X=table(); return;
end
X=cell2table(rows,'VariableNames',{'patient','drug_A','drug_B','n_shared_states', 'C_within_A','C_within_B','C_between_different_states','C_between_same_state'});
X.patient=string(X.patient); X.drug_A=string(X.drug_A); X.drug_B=string(X.drug_B);
numvars=X.Properties.VariableNames(4:end);
for i=1:numel(numvars)
    if iscell(X.(numvars{i})), X.(numvars{i})=cell2mat(X.(numvars{i})); end
end
end

function S = summarizeCrossDrugCoherence(X)
if isempty(X)
    S=table(); return;
end
within=0.5*(X.C_within_A+X.C_within_B);
between=X.C_between_different_states;
same=X.C_between_same_state;
S=table(height(X),numel(unique(X.patient)),median(within,'omitnan'), median(between,'omitnan'),median(same,'omitnan'), median(within-between,'omitnan'),mean(between<within), 'VariableNames',{'n_drug_pairs','n_units','median_within_drug_C', 'median_between_drug_C','median_between_drug_same_state_C', 'median_within_minus_between','fraction_within_exceeds_between'});
end


function [Unit,Global] = summarizeCrossDrugByUnit(X,cfg)
% Summarize cross-drug specificity using TissueID, not drug-pair, as the
% independent biological unit. This prevents the 52 drug pairs from being
% treated as 52 independent replicates when they are nested within units.

if isempty(X)
    Unit=table(); Global=table(); return;
end

pats=unique(X.patient);
rows={};

for ip=1:numel(pats)
    q=X.patient==pats(ip);

    within=0.5*(X.C_within_A(q)+X.C_within_B(q));
    betweenDiff=X.C_between_different_states(q);
    betweenSame=X.C_between_same_state(q);

    rows(end+1,:)={pats(ip),nnz(q), ...
        median(within,'omitnan'), ...
        median(betweenDiff,'omitnan'), ...
        median(betweenSame,'omitnan'), ...
        median(within-betweenDiff,'omitnan'), ...
        median(within-betweenSame,'omitnan')};
end

Unit=cell2table(rows,'VariableNames', ...
    {'patient','n_drug_pairs', ...
     'median_within_drug_C', ...
     'median_between_different_state_C', ...
     'median_between_same_state_C', ...
     'delta_within_minus_between_different', ...
     'delta_within_minus_between_same'});

Unit.patient=string(Unit.patient);
for i=2:width(Unit)
    vn=Unit.Properties.VariableNames{i};
    if iscell(Unit.(vn))
        Unit.(vn)=cell2mat(Unit.(vn));
    end
end

comparisons={ ...
    'delta_within_minus_between_different','within minus between-drug, different-state'; ...
    'delta_within_minus_between_same','within minus between-drug, same-state'};

grows={};

for i=1:size(comparisons,1)
    x=Unit.(comparisons{i,1});
    x=x(isfinite(x));

    [lo,hi]=bootstrapMedianCI(x,cfg.nBootstrap);

    nPositive=nnz(x>0);
    nNegative=nnz(x<0);

    % One-sided exact sign test for the prespecified direction:
    % within-drug coherence > between-drug coherence.
    pOne=exactSignPOneSidedPositive(x);

    grows(end+1,:)={string(comparisons{i,2}),numel(x), ...
        median(x,'omitnan'),lo,hi,nPositive,nNegative,pOne};
end

Global=cell2table(grows,'VariableNames', ...
    {'comparison','n_units','median_delta','CI_low','CI_high', ...
     'n_positive','n_negative','sign_p_one_sided'});

Global.comparison=string(Global.comparison);
for i=2:width(Global)
    vn=Global.Properties.VariableNames{i};
    if iscell(Global.(vn))
        Global.(vn)=cell2mat(Global.(vn));
    end
end
end

function p=exactSignPOneSidedPositive(x)
x=x(isfinite(x) & x~=0);
n=numel(x);

if n==0
    p=NaN;
    return;
end

k=nnz(x>0);

% P[X >= k] for X ~ Binomial(n,0.5), implemented without Statistics Toolbox.
p=0;
for j=k:n
    p=p+nchoosek(n,j)*(0.5^n);
end

p=min(1,p);
end

function RotNull = computeRotationNull(patient,treatment,mask,state,Z,R,states,cfg)
rows={};
G=size(Z,2);
for i=1:height(R)
    p=R.patient(i); d=R.drug(i);
    ctrl=mask & patient==p & treatment=="vehicle";
    tr=mask & patient==p & treatment==d;
    used=strsplit(R.states_used{i},',');
    normsK=[];
    for u=1:numel(used)
        k=find(strcmp(states,strtrim(used{u})));
        if isempty(k), continue; end
        a=ctrl & state==k; b=tr & state==k;
        if nnz(a)>=cfg.minCellsPerState && nnz(b)>=cfg.minCellsPerState
            delta=mean(Z(b,:),1,'omitnan')-mean(Z(a,:),1,'omitnan');
            normsK(end+1)=norm(delta);
        end
    end
    ns=numel(normsK);
    if ns<cfg.minStatesForMetric
        rows(end+1,:)={p,d,ns,R.coherence_C(i),NaN,NaN,NaN,NaN};
        continue;
    end
    cNull=nan(cfg.nNullDraws,1);
    for it=1:cfg.nNullDraws
        Delta=zeros(ns,G);
        for u=1:ns
            v=randn(1,G);
            v=v/norm(v);
            Delta(u,:)=v*normsK(u);
        end
        cNull(it)=meanPairwiseCosine(Delta,cfg.eps);
    end
    cq=prctile(cNull,[2.5 50 97.5]);
    pC=(1+nnz(cNull>=R.coherence_C(i)))/(cfg.nNullDraws+1);
    rows(end+1,:)={p,d,ns,R.coherence_C(i),cq(2),cq(1),cq(3),pC};
    fprintf('  rotation null: %s / %s done\n',char(p),char(d));
end
RotNull=cell2table(rows,'VariableNames',{'patient','drug','n_states', 'coherence_C','null_C_median','null_C_lo','null_C_hi','p_C'});
RotNull.patient=string(RotNull.patient); RotNull.drug=string(RotNull.drug);
numvars=RotNull.Properties.VariableNames(3:end);
for i=1:numel(numvars)
    if iscell(RotNull.(numvars{i})), RotNull.(numvars{i})=cell2mat(RotNull.(numvars{i})); end
end
end

function S = summarizeRotationNull(RotNull)
ok=isfinite(RotNull.p_C) & isfinite(RotNull.coherence_C) & isfinite(RotNull.null_C_median);
S=table(nnz(ok),median(RotNull.coherence_C(ok)),median(RotNull.null_C_median(ok)), median(RotNull.coherence_C(ok)-RotNull.null_C_median(ok)), nnz(RotNull.p_C(ok)<0.05),mean(RotNull.p_C(ok)<0.05), 'VariableNames',{'n_comparisons','median_observed','median_null', 'median_excess','n_significant_0p05','fraction_significant_0p05'});
end

function DrugSign = computeDrugStratifiedSignTest(R)
drugs=unique(R.drug);
rows={};
for i=1:numel(drugs)
    d=drugs(i);
    k=R.contraction_K(R.drug==d);
    k=k(k~=0);
    n=numel(k);
    if n==0
        continue;
    end
    nPos=nnz(k>0);
    vote=sign(nPos-n/2);
    rows(end+1,:)={d,n,nPos,vote};
end
DrugSign=cell2table(rows,'VariableNames',{'drug','n_units','n_positive','vote'});
DrugSign.drug=string(DrugSign.drug);
v=DrugSign.vote;
v=v(v~=0);
nD=numel(v);
nPosD=nnz(v>0);
pDrug=exactSignP([ones(nPosD,1);-ones(nD-nPosD,1)]);
fprintf('  drug-stratified sign test: %d of %d drugs vote K>0 (median unit), P=%.4f\n', nPosD,nD,pDrug);
end
