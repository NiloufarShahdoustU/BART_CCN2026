% the the final version, mix of 5 and 2
% Author: Nill

clear;
clc;
close all;
warning('off','all');

inputFolderName_IEDdata = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_23_IED_nonIED_ITs_2_chunks\';
fileList = dir(fullfile(inputFolderName_IEDdata, '*.LFPIED.mat'));
PatientsNum = length(fileList);

PostResponse = 1;
PreOutcome   = 2;
alpha = 0.05;

epoch_names = {'Post-Response','Pre-Outcome'};
epoch_colors = [
    0.85 0.33 0.10
    0.60 0.40 0.80
];

nonIEDtrials_bhvr_measure_mean = nan(PatientsNum,1);
IEDtrials_bhvr_measure_mean    = nan(PatientsNum,2);

for pt = 1:PatientsNum
    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1};
    disp("patient: " + ptID);

    IEDdata = [inputFolderName_IEDdata '\' ptID '.LFPIED.mat'];
    load(IEDdata);

    nChans = length(LFPIED.selectedChans)-1;


    RTs = LFPIED.RTs;
    ITs = LFPIED.ITs;
    RTsThreshold = 10;
    OutlierIndices = RTs >= RTsThreshold;

    isControl = LFPIED.isControl;      % 1 = control, 0 = not control
    isControl = isControl(:)';        
    keepIdx   = (~OutlierIndices) & (~isControl);

    RTs = RTs(keepIdx);
    ITs = ITs(keepIdx);

    IEDtrialsPostResponse = LFPIED.IEDtrialsPostResponse(:, keepIdx);
    IEDtrialsPreOutcome   = LFPIED.IEDtrialsPreOutcome(:, keepIdx);


    allTimePoints = IEDtrialsPostResponse + IEDtrialsPreOutcome;
    nonIEDIndices = find(all(allTimePoints == 0, 1));

    nonIEDtrials_bhvr_measure = ITs(nonIEDIndices);
    nonIEDtrials_bhvr_measure_mean(pt) = mean(nonIEDtrials_bhvr_measure);

    IEDTrials_bhvr_measure_MeanPerChan_PostResponse = nan(nChans,1);
    IEDTrials_bhvr_measure_MeanPerChan_PreOutcome   = nan(nChans,1);

    IEDtrialsPostResponse_only = IEDtrialsPostResponse .* ~IEDtrialsPreOutcome;
    IEDtrialsPreOutcome_only   = IEDtrialsPreOutcome   .* ~IEDtrialsPostResponse;

    pVal_PostResponse = nan(1,nChans);
    pVal_PreOutcome   = nan(1,nChans);

    NumberofPermutations = 10000;

    for chz = 1:nChans

        IEDTrials_bhvr_measure = IEDtrialsPostResponse_only(chz,:);
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure .* ITs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);

        if (numel(IEDTrials_bhvr_measure) > 0)
            pVal_PostResponse(chz) = permutationTest(IEDTrials_bhvr_measure, nonIEDtrials_bhvr_measure, NumberofPermutations);
            if pVal_PostResponse(chz) < alpha
                IEDTrials_bhvr_measure_MeanPerChan_PostResponse(chz) = nanmean(IEDTrials_bhvr_measure);
            end
        else
            pVal_PostResponse(chz) = NaN;
        end
        clear IEDTrials_bhvr_measure

        IEDTrials_bhvr_measure = IEDtrialsPreOutcome_only(chz,:);
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure .* ITs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);

        if (numel(IEDTrials_bhvr_measure) > 0)
            pVal_PreOutcome(chz) = permutationTest(IEDTrials_bhvr_measure, nonIEDtrials_bhvr_measure, NumberofPermutations);
            if pVal_PreOutcome(chz) < alpha
                IEDTrials_bhvr_measure_MeanPerChan_PreOutcome(chz) = nanmean(IEDTrials_bhvr_measure);
            end
        else
            pVal_PreOutcome(chz) = NaN;
        end
        clear IEDTrials_bhvr_measure

    end

    IEDtrials_bhvr_measure_mean(pt,PostResponse) = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PostResponse);
    IEDtrials_bhvr_measure_mean(pt,PreOutcome)   = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PreOutcome);

end

%%

name = 'IEDtrials_bhvr_measure_mean_ITs_PostResponse_PreOutcome';
save([outputFolderName name '.mat'],'IEDtrials_bhvr_measure_mean');

IEDtrials_bhvr_measure_mean_PostResponse = IEDtrials_bhvr_measure_mean(:,PostResponse);
IEDtrials_bhvr_measure_mean_PostResponse = IEDtrials_bhvr_measure_mean_PostResponse(~isnan(IEDtrials_bhvr_measure_mean_PostResponse));
IEDtrials_bhvr_measure_mean_PostResponse = IEDtrials_bhvr_measure_mean_PostResponse(IEDtrials_bhvr_measure_mean_PostResponse ~= 0);

IEDtrials_bhvr_measure_mean_PreOutcome = IEDtrials_bhvr_measure_mean(:,PreOutcome);
IEDtrials_bhvr_measure_mean_PreOutcome = IEDtrials_bhvr_measure_mean_PreOutcome(~isnan(IEDtrials_bhvr_measure_mean_PreOutcome));
IEDtrials_bhvr_measure_mean_PreOutcome = IEDtrials_bhvr_measure_mean_PreOutcome(IEDtrials_bhvr_measure_mean_PreOutcome ~= 0);

vec1 = nonIEDtrials_bhvr_measure_mean;
vec2 = IEDtrials_bhvr_measure_mean_PostResponse;
vec3 = IEDtrials_bhvr_measure_mean_PreOutcome;

allVecs = [vec1; vec2; vec3];
group = [ones(length(vec1),1); 2*ones(length(vec2),1); 3*ones(length(vec3),1)];

figSize = 4;

figure( ...
    'Units','inches', ...
    'Position',[1 1 figSize figSize], ...   
    'Visible','on');

boxplot(allVecs, group, ...
    'Labels', {'non-IED','Post-Response','Pre-Outcome'}, ...
    'Color', 'k', ...
    'Symbol', '');

xlabel('non-IED trials vs IED trials', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('normalized mean IT', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

set(findobj(gca, 'Type', 'line'), 'LineWidth', 0.5);

hold on;
jitterAmount = 0.1;
markerSize = 18;

x1 = 1 + (rand(size(vec1)) - 0.5) * jitterAmount;
scatter(x1, vec1, markerSize, ...
    'filled', ...
    'MarkerFaceColor', [0.5 0.5 0.5], ...
    'MarkerFaceAlpha', 0.4, ...
    'MarkerEdgeColor', 'none');

x2 = 2 + (rand(size(vec2)) - 0.5) * jitterAmount;
scatter(x2, vec2, markerSize, ...
    'filled', ...
    'MarkerFaceColor', epoch_colors(1,:), ...
    'MarkerFaceAlpha', 0.4, ...
    'MarkerEdgeColor', 'none');

x3 = 3 + (rand(size(vec3)) - 0.5) * jitterAmount;
scatter(x3, vec3, markerSize, ...
    'filled', ...
    'MarkerFaceColor', epoch_colors(2,:), ...
    'MarkerFaceAlpha', 0.4, ...
    'MarkerEdgeColor', 'none');

hold off;

set(gca, 'box', 'off', 'tickdir', 'out');
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 figSize figSize], ...
         'PaperSize', [figSize figSize]);
filename = 'ITs_boxplot_PostResponse_PreOutcome';
saveas(gcf, fullfile(outputFolderName, filename), 'pdf');

pValuePostResponse_ITs_pt = ranksum(vec1,vec2);
pValuePreOutcome_ITs_pt   = ranksum(vec1,vec3);

outputFileName = fullfile(outputFolderName, 'ITs_pValues_PostResponse_PreOutcome.txt');
fileID = fopen(outputFileName, 'w');

fprintf(fileID, 'P-values from Rank Sum Tests:\n');
fprintf(fileID, 'P-value (non-IED vs Post-Response): %.4f\n', pValuePostResponse_ITs_pt);
fprintf(fileID, 'P-value (non-IED vs Pre-Outcome): %.4f\n', pValuePreOutcome_ITs_pt);

fclose(fileID);
