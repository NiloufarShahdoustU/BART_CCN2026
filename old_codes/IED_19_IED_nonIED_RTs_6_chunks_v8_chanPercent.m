 % the the final version, mix of 5 and 2

% Author: Nill



clear;
clc;
close all;
warning('off','all');


%%


inputFolderName_IEDdata = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks_IED_11_v9';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_19_IEDnonIED_bhvr_analysis\v8\';

fileList = dir(fullfile(inputFolderName_IEDdata, '*.LFPIED.mat'));
PatientsNum = length(fileList);



PreOnset1 = 1;
PreOnset2 = 2;
PostOnset = 3;
PreResponse = 4;
PostResponse = 5;
PreOutcome = 6;
alpha = 0.05;
%%

% for 4 time periods for RTs

% these means are over trials 
nonIEDtrials_bhvr_measure_mean = nan(PatientsNum,1); 
IEDtrials_bhvr_measure_mean = nan(PatientsNum, 4);


AnatoicalLocsNums = 150; % arbitrary number
AnatomicalLocsPatientsPreOnset1 = zeros(AnatoicalLocsNums, PatientsNum);
AnatomicalLocsPatientsPreOnset2 = zeros(AnatoicalLocsNums, PatientsNum);
AnatomicalLocsPatientsPostOnset = zeros(AnatoicalLocsNums, PatientsNum);
AnatomicalLocsPatientsPreResponse = zeros(AnatoicalLocsNums, PatientsNum);
AnatomicalLocsPatientsAll = zeros(AnatoicalLocsNums, PatientsNum);

nan_cell_array = repmat({'nan'},AnatoicalLocsNums, 1);
AnatomicalLocsVecPreOnset1 = string(nan_cell_array);
AnatomicalLocsVecPreOnset2 = string(nan_cell_array);
AnatomicalLocsVecPostOnset = string(nan_cell_array);
AnatomicalLocsVecPreResponse = string(nan_cell_array);
AnatomicalLocsVecAll = string(nan_cell_array);


for pt = 1:PatientsNum
% for pt = 2:3
 fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1}; 
    disp("patient: " + ptID);

    %IED data read
    IEDdata = [inputFolderName_IEDdata '\' ptID '.LFPIED.mat'];
    load(IEDdata);

    nChans = length(LFPIED.selectedChans)-1; % removing the last selected chan

    anatomicalLocs = LFPIED.anatomicalLocs;

    selectedChans = LFPIED.selectedChans;
    selectedChans = selectedChans(1:end-1);
    RTs = LFPIED.RTs;
    RTsThreshold = 10;
    OutlierIndices = RTs >= RTsThreshold;
    RTs = RTs(~OutlierIndices);
    nTrials = length(RTs);

    IEDtrialsPreOnset1 = LFPIED.IEDtrialsPreOnset1(:, ~OutlierIndices);
    IEDtrialsPreOnset2 = LFPIED.IEDtrialsPreOnset2(:, ~OutlierIndices);
    IEDtrialsPostOnset = LFPIED.IEDtrialsPostOnset(:, ~OutlierIndices);
    IEDtrialsPreResponse = LFPIED.IEDtrialsPreResponse(:, ~OutlierIndices);
    % IEDtrialsPostResponse = LFPIED.IEDtrialsPostResponse(:, ~OutlierIndices);
    % IEDtrialsPreOutcome = LFPIED.IEDtrialsPreOutcome(:, ~OutlierIndices);

        % first I need to find the nonIED trials, I mean the trials that did
    % not have ANY IEDs in ANY timepoints before Response!!! in ANY chans! for this I am going
    % to sum all the IEDtrials a pairwise sum, then if all the rows of a
    % colum is 0, that colum (trial) is an nonIED trial! Pay attention I am
    % not taking into account the 
% 
%     allTimePoints = IEDtrialsPreOnset1 + IEDtrialsPreOnset2 + IEDtrialsPostOnset + IEDtrialsPreResponse;
%     nonIEDIndices = find(all(allTimePoints == 0, 1));

        allTimePoints = IEDtrialsPreOnset1 + IEDtrialsPreOnset2 + IEDtrialsPostOnset + IEDtrialsPreResponse;
        nonIEDIndices = find(all(allTimePoints == 0, 1));

   

        % now that I've found nonIED trials I need to have their behavioral
    % mearue vector: 
    % for the sake of each patient
    nonIEDtrials_bhvr_measure = RTs(nonIEDIndices);
    % for the sake of overall patients
    % I am doing standardized RTs:
    % the calculated means for each patient are influenced by the number of trials and the 
    % specific variability of that patient%s RTs. To address this, you can consider z-scoring
    % or standardizing the RTs for each patient before taking the mean. This ensures that differences
    % in baseline RTs and variability across patients are accounted for.
   


    % Step 2: Use the z-scored RTs to compute the mean for non-IED trials
    nonIEDtrials_bhvr_measure_mean(pt) = mean(nonIEDtrials_bhvr_measure); 



    IEDTrials_bhvr_measure_MeanPerChan_PreOnset1 = nan(nChans,1);
    IEDTrials_bhvr_measure_MeanPerChan_PreOnset2 = nan(nChans,1);
    IEDTrials_bhvr_measure_MeanPerChan_PostOnset = nan(nChans,1);
    IEDTrials_bhvr_measure_MeanPerChan_PreResponse = nan(nChans,1);
    % IEDTrials_bhvr_measure_MeanPerChan_PostResponse = nan(nChans,1);
    % IEDTrials_bhvr_measure_MeanPerChan_PreOutcome = nan(nChans,1);



    % now I need to take those IEDtrials that ONLY happened in 1 of the 4
    % time periods and not in the other 3:
    IEDtrialsPreOnset1_only = IEDtrialsPreOnset1.*~IEDtrialsPreOnset2.*~IEDtrialsPostOnset.*~IEDtrialsPreResponse;
    IEDtrialsPreOnset2_only = IEDtrialsPreOnset2.*~IEDtrialsPreOnset1.*~IEDtrialsPostOnset.*~IEDtrialsPreResponse;
    IEDtrialsPostOnset_only = IEDtrialsPostOnset.*~IEDtrialsPreOnset1.*~IEDtrialsPreOnset2.*~IEDtrialsPreResponse;
    IEDtrialsPreResponse_only = IEDtrialsPreResponse.*~IEDtrialsPreOnset1.*~IEDtrialsPreOnset2.*~IEDtrialsPostOnset;

    

    pVal_PreOnset1 = nan(1,nChans);
    pVal_PreOnset2 = nan(1,nChans);
    pVal_PostOnset = nan(1,nChans);
    pVal_PreResponse = nan(1,nChans);
    % pVal_PostResponse = nan(1,nChans);
    % pVal_PreOutcome = nan(1,nChans);

    NumberofPermutations = 10000;


    % stacked behavior measure for ied and non ied trials:
    

    

    for chz = 1:nChans
  
    
        %PreOnset1
        IEDTrials_bhvr_measure = IEDtrialsPreOnset1_only(chz,:); 
         % we are only taking the RTs that are not outliers
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*RTs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);


        if (size(IEDTrials_bhvr_measure)>0)
            pVal_PreOnset1(chz) = permutationTest(IEDTrials_bhvr_measure, nonIEDtrials_bhvr_measure, NumberofPermutations);
            if pVal_PreOnset1(chz)<alpha
                IEDTrials_bhvr_measure_MeanPerChan_PreOnset1(chz) = nanmean(IEDTrials_bhvr_measure);
            end
        else
            pVal_PreOnset1(chz) = NaN;
        end
        clear IEDTrials_bhvr_measure 

        %///////////////////////////////////////////////////////////////////////////////////////////////////////////

         %PreOnset2

  

        IEDTrials_bhvr_measure = IEDtrialsPreOnset2_only(chz,:); 
         % we are only taking the RTs that are not outliers
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*RTs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);


        if (size(IEDTrials_bhvr_measure)>0)
            pVal_PreOnset2(chz) = permutationTest(IEDTrials_bhvr_measure, nonIEDtrials_bhvr_measure, NumberofPermutations);
            if pVal_PreOnset2(chz)<alpha
                IEDTrials_bhvr_measure_MeanPerChan_PreOnset2(chz) = nanmean(IEDTrials_bhvr_measure);
            end
            else
            pVal_PreOnset2(chz) = NaN;
        end
        clear IEDTrials_bhvr_measure 


        %///////////////////////////////////////////////////////////////////////////////////////////////////////////

            %PostOnset


        IEDTrials_bhvr_measure = IEDtrialsPostOnset_only(chz,:); 
         % we are only taking the RTs that are not outliers
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*RTs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);



        if (size(IEDTrials_bhvr_measure)>0)
            pVal_PostOnset(chz) = permutationTest(IEDTrials_bhvr_measure, nonIEDtrials_bhvr_measure, NumberofPermutations);
            if pVal_PostOnset(chz)< alpha
                IEDTrials_bhvr_measure_MeanPerChan_PostOnset(chz) = nanmean(IEDTrials_bhvr_measure);
            end
        else
            pVal_PostOnset(chz) = NaN;
        end
        clear IEDTrials_bhvr_measure 

        %///////////////////////////////////////////////////////////////////////////////////////////////////////////

            %PreResponse


        IEDTrials_bhvr_measure = IEDtrialsPreResponse_only(chz,:); 
         % we are only taking the RTs that are not outliers
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*RTs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);


        if (size(IEDTrials_bhvr_measure)>0)
            pVal_PreResponse(chz) = permutationTest(IEDTrials_bhvr_measure, nonIEDtrials_bhvr_measure, NumberofPermutations);
            if pVal_PreResponse(chz)< alpha
                IEDTrials_bhvr_measure_MeanPerChan_PreResponse(chz) = nanmean(IEDTrials_bhvr_measure);
            end
            else
            pVal_PreResponse(chz) = NaN;
        end
        clear IEDTrials_bhvr_measure 

        %///////////////////////////////////////////////////////////////////////////////////////////////////////////


   
  
    end % end for chan



     IEDtrials_bhvr_measure_mean(pt,PreOnset1) = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PreOnset1);
    IEDtrials_bhvr_measure_mean(pt,PreOnset2) = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PreOnset2);
     IEDtrials_bhvr_measure_mean(pt,PostOnset) = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PostOnset);
    IEDtrials_bhvr_measure_mean(pt,PreResponse) = nanmean(IEDTrials_bhvr_measure_MeanPerChan_PreResponse);

%////////////////////////////////////////////////////////////////////////////////
    
    % for the sake of channel percentage

    % PreOnset1
    ChanIndices = find(~isnan(IEDTrials_bhvr_measure_MeanPerChan_PreOnset1));
    timePeriodAnatomicalLoc = anatomicalLocs(selectedChans(ChanIndices));
    for location=1:length(timePeriodAnatomicalLoc)
        element = timePeriodAnatomicalLoc(location);
        tempIndexInLocs = ismember(AnatomicalLocsVecPreOnset1, element);
        FoundIndexInLocs = find(tempIndexInLocs);
        if ~isempty(FoundIndexInLocs)
            AnatomicalLocsPatientsPreOnset1(FoundIndexInLocs,pt) = AnatomicalLocsPatientsPreOnset1(FoundIndexInLocs,pt)+1;
        else
            nan_index = find(AnatomicalLocsVecPreOnset1 == "nan", 1);
            AnatomicalLocsVecPreOnset1(nan_index) = element;
            AnatomicalLocsPatientsPreOnset1(nan_index,pt) = AnatomicalLocsPatientsPreOnset1(nan_index,pt)+1;
        end
        clear FoundIndexInLocs nan_index tempIndexInLocs element
    end
    clear ChanIndices timePeriodAnatomicalLoc



        % PreOnset2
    ChanIndices = find(~isnan(IEDTrials_bhvr_measure_MeanPerChan_PreOnset2));
    timePeriodAnatomicalLoc = anatomicalLocs(selectedChans(ChanIndices));
    for location=1:length(timePeriodAnatomicalLoc)
        element = timePeriodAnatomicalLoc(location);
        tempIndexInLocs = ismember(AnatomicalLocsVecPreOnset2, element);
        FoundIndexInLocs = find(tempIndexInLocs);
        if ~isempty(FoundIndexInLocs)
            AnatomicalLocsPatientsPreOnset2(FoundIndexInLocs,pt) = AnatomicalLocsPatientsPreOnset2(FoundIndexInLocs,pt)+1;
        else
            nan_index = find(AnatomicalLocsVecPreOnset2 == "nan", 1);
            AnatomicalLocsVecPreOnset2(nan_index) = element;
            AnatomicalLocsPatientsPreOnset2(nan_index,pt) = AnatomicalLocsPatientsPreOnset2(nan_index,pt)+1;
        end
        clear FoundIndexInLocs nan_index tempIndexInLocs element
    end
    clear ChanIndices timePeriodAnatomicalLoc


            % PostOnset
    ChanIndices = find(~isnan(IEDTrials_bhvr_measure_MeanPerChan_PostOnset));
    timePeriodAnatomicalLoc = anatomicalLocs(selectedChans(ChanIndices));
    for location=1:length(timePeriodAnatomicalLoc)
        element = timePeriodAnatomicalLoc(location);
        tempIndexInLocs = ismember(AnatomicalLocsVecPostOnset, element);
        FoundIndexInLocs = find(tempIndexInLocs);
        if ~isempty(FoundIndexInLocs)
            AnatomicalLocsPatientsPostOnset(FoundIndexInLocs,pt) = AnatomicalLocsPatientsPostOnset(FoundIndexInLocs,pt)+1;
        else
            nan_index = find(AnatomicalLocsVecPostOnset == "nan", 1);
            AnatomicalLocsVecPostOnset(nan_index) = element;
            AnatomicalLocsPatientsPostOnset(nan_index,pt) = AnatomicalLocsPatientsPostOnset(nan_index,pt)+1;
        end
        clear FoundIndexInLocs nan_index tempIndexInLocs element
    end
    clear ChanIndices timePeriodAnatomicalLoc



                % PreResponse
    ChanIndices = find(~isnan(IEDTrials_bhvr_measure_MeanPerChan_PreResponse));
    timePeriodAnatomicalLoc = anatomicalLocs(selectedChans(ChanIndices));
    for location=1:length(timePeriodAnatomicalLoc)
        element = timePeriodAnatomicalLoc(location);
        tempIndexInLocs = ismember(AnatomicalLocsVecPreResponse, element);
        FoundIndexInLocs = find(tempIndexInLocs);
        if ~isempty(FoundIndexInLocs)
            AnatomicalLocsPatientsPreResponse(FoundIndexInLocs,pt) = AnatomicalLocsPatientsPreResponse(FoundIndexInLocs,pt)+1;
        else
            nan_index = find(AnatomicalLocsVecPreResponse == "nan", 1);
            AnatomicalLocsVecPreResponse(nan_index) = element;
            AnatomicalLocsPatientsPreResponse(nan_index,pt) = AnatomicalLocsPatientsPreResponse(nan_index,pt)+1;
        end
        clear FoundIndexInLocs nan_index tempIndexInLocs element
    end
    clear ChanIndices timePeriodAnatomicalLoc



    clear IEDtrialsPreOnset1  IEDtrialsPreOnset2 IEDtrialsPostOnset IEDtrialsPreResponse 
    clear IEDtrialsPreOnset1_only  IEDtrialsPreOnset2_only IEDtrialsPostOnset_only IEDtrialsPreResponse_only 
end



for pt = 1:PatientsNum

    AnatomicalLocsAll = anatomicalLocs(selectedChans);
    for all=1:length(AnatomicalLocsAll)
        element = AnatomicalLocsAll(all);
        tempIndexInLocs = ismember(AnatomicalLocsVecAll, element);
        FoundIndexInLocs = find(tempIndexInLocs);
        if ~isempty(FoundIndexInLocs)
            AnatomicalLocsPatientsAll(FoundIndexInLocs,pt) = AnatomicalLocsPatientsAll(FoundIndexInLocs,pt)+1;
        else
            nan_index = find(AnatomicalLocsVecAll == "nan", 1);
            AnatomicalLocsVecAll(nan_index) = element;
            AnatomicalLocsPatientsAll(nan_index,pt) = AnatomicalLocsPatientsAll(nan_index,pt)+1;
        end
       clear FoundIndexInLocs nan_index tempIndexInLocs element
    end
    clear AnatomicalLocsAll

end

%%

name = 'IEDtrials_bhvr_measure_mean_RTs';
save([outputFolderName name '.mat'],'IEDtrials_bhvr_measure_mean');
%% preprocessing

IEDtrials_bhvr_measure_mean_PreOnset1 = IEDtrials_bhvr_measure_mean(:,PreOnset1);
IEDtrials_bhvr_measure_mean_PreOnset1 = IEDtrials_bhvr_measure_mean_PreOnset1(~isnan(IEDtrials_bhvr_measure_mean_PreOnset1));


IEDtrials_bhvr_measure_mean_PreOnset2 = IEDtrials_bhvr_measure_mean(:,PreOnset2);
IEDtrials_bhvr_measure_mean_PreOnset2 = IEDtrials_bhvr_measure_mean_PreOnset2(~isnan(IEDtrials_bhvr_measure_mean_PreOnset2));


IEDtrials_bhvr_measure_mean_PostOnset = IEDtrials_bhvr_measure_mean(:,PostOnset);
IEDtrials_bhvr_measure_mean_PostOnset = IEDtrials_bhvr_measure_mean_PostOnset(~isnan(IEDtrials_bhvr_measure_mean_PostOnset));

IEDtrials_bhvr_measure_mean_PreResponse = IEDtrials_bhvr_measure_mean(:,PreResponse);
IEDtrials_bhvr_measure_mean_PreResponse = IEDtrials_bhvr_measure_mean_PreResponse(~isnan(IEDtrials_bhvr_measure_mean_PreResponse));





%% visualization box plot


% for RTs only the first 4 time periods, are important: 

vec1 = nonIEDtrials_bhvr_measure_mean;
vec2 = IEDtrials_bhvr_measure_mean_PreOnset1;
vec3 = IEDtrials_bhvr_measure_mean_PreOnset2;
vec4 = IEDtrials_bhvr_measure_mean_PostOnset;
vec5 = IEDtrials_bhvr_measure_mean_PreResponse;


% Concatenate all vectors into a single column vector and create group identifiers
allVecs = [vec1; vec2; vec3; vec4; vec5];
group = [ones(length(vec1),1); 2*ones(length(vec2),1); 3*ones(length(vec3),1); 4*ones(length(vec4),1); 5*ones(length(vec5),1) ];

figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.2, 0.4], 'Visible', 'on'); 

% Plot the boxplot with black color
boxplotHandle = boxplot(allVecs, group, 'Labels', {'non-IED', 'PreOnset1', 'PreOnset2', 'PostOnset', 'PreResponse'}, 'Color', 'k', 'Symbol', '');
xlabel('non-IED trials vs IED trials', 'FontSize', 14, 'FontWeight', 'bold');
% ylim([0 5]);
ylabel('mean RT', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

% Make box plot lines thicker
set(findobj(gca, 'Type', 'line'), 'LineWidth', 0.5);

% Overlay jittered data points in gray with smaller marker size
hold on;
jitterAmount = 0.1;
markerSize = 10; % Smaller marker size
scatter(group + (rand(size(group)) - 0.5) * jitterAmount, allVecs, markerSize, ...
        'MarkerEdgeColor', [0.5 0.5 0.5], 'jitter', 'on', 'jitterAmount', jitterAmount);


% ylimValues = ylim; % Get current y-axis limits
% line([1.5 1.5], ylimValues, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); % Draw dashed line

% Final adjustments and saving the plot
hold off;
set(gca, 'box', 'off', 'tickdir', 'out');
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
filename = 'RTs_boxplot';
saveas(gcf, fullfile(outputFolderName, filename), 'pdf');


%% p_val for bhv measures accross patients

pValuePreOnset1_acc_pt = ranksum(vec1,vec2);
pValuePreOnset2_acc_pt = ranksum(vec1,vec3);
pValuePostOnset_acc_pt = ranksum(vec1,vec4);
pValuePreResponse_acc_pt = ranksum(vec1,vec5);

% Specify the output file name
outputFileName = fullfile(outputFolderName, 'RTs_pValues.txt');

% Open the file for writing
fileID = fopen(outputFileName, 'w');

% Write p-values to the file
fprintf(fileID, 'P-values from Rank Sum Tests:\n');
fprintf(fileID, 'P-value (non-IED vs PreOnset1): %.4f\n', pValuePreOnset1_acc_pt);
fprintf(fileID, 'P-value (non-IED vs PreOnset2): %.4f\n', pValuePreOnset2_acc_pt);
fprintf(fileID, 'P-value (non-IED vs PostOnset): %.4f\n', pValuePostOnset_acc_pt);
fprintf(fileID, 'P-value (non-IED vs PreResponse): %.4f\n', pValuePreResponse_acc_pt);

% Close the file
fclose(fileID);


%% visualization channel percentage

%cleaning AnatomicalLocsVec

% Replace elements starting with 'NaC' or containing 'Left Cerebral White Matter'
startsWithNaC = startsWith(AnatomicalLocsVecPreOnset1, "NaC");
containsLCWM = contains(AnatomicalLocsVecPreOnset1, "Left Cerebral White Matter");

% Combine both conditions and replace with 'nan'
AnatomicalLocsVecPreOnset1(startsWithNaC | containsLCWM) = "nan";
clear startsWithNaC containsLCWM

startsWithNaC = startsWith(AnatomicalLocsVecPreOnset2, "NaC");
containsLCWM = contains(AnatomicalLocsVecPreOnset2, "Left Cerebral White Matter");
AnatomicalLocsVecPreOnset2(startsWithNaC | containsLCWM) = "nan";
clear startsWithNaC containsLCWM

startsWithNaC = startsWith(AnatomicalLocsVecPostOnset, "NaC");
containsLCWM = contains(AnatomicalLocsVecPostOnset, "Left Cerebral White Matter");
AnatomicalLocsVecPostOnset(startsWithNaC | containsLCWM) = "nan";
clear startsWithNaC containsLCWM

startsWithNaC = startsWith(AnatomicalLocsVecPreResponse, "NaC");
containsLCWM = contains(AnatomicalLocsVecPreResponse, "Left Cerebral White Matter");
AnatomicalLocsVecPreResponse(startsWithNaC | containsLCWM) = "nan";
clear startsWithNaC containsLCWM

startsWithNaC = startsWith(AnatomicalLocsVecAll, "NaC");
containsLCWM = contains(AnatomicalLocsVecAll, "Left Cerebral White Matter");
AnatomicalLocsVecAll(startsWithNaC | containsLCWM) = "nan";
clear startsWithNaC containsLCWM






% cleaning AnatomicalLocsPatients based on AnatomicalLocsVec

missingIndices = find(AnatomicalLocsVecPreOnset1 == "nan");
AnatomicalLocsPatientsPreOnset1(missingIndices, :) = [];
AnatomicalLocsVecPreOnset1(missingIndices) = [];
clear missingIndices 




missingIndices = find(AnatomicalLocsVecPreOnset2 == "nan");
AnatomicalLocsPatientsPreOnset2(missingIndices, :) = [];
AnatomicalLocsVecPreOnset2(missingIndices) = [];
clear missingIndices 


missingIndices = find(AnatomicalLocsVecPostOnset == "nan");
AnatomicalLocsPatientsPostOnset(missingIndices, :) = [];
AnatomicalLocsVecPostOnset(missingIndices) = [];
clear missingIndices 


missingIndices = find(AnatomicalLocsVecPreResponse == "nan");
AnatomicalLocsPatientsPreResponse(missingIndices, :) = [];
AnatomicalLocsVecPreResponse(missingIndices) = [];
clear missingIndices 


missingIndicesAll = find(AnatomicalLocsVecAll == "nan");
AnatomicalLocsPatientsAll(missingIndicesAll, :) = [];
AnatomicalLocsVecAll(missingIndicesAll) = [];






ChanNumsPreOnset1 = sum(AnatomicalLocsPatientsPreOnset1, 2);
ChanNumsPreOnset2 = sum(AnatomicalLocsPatientsPreOnset2, 2);
ChanNumsPostOnset = sum(AnatomicalLocsPatientsPostOnset, 2);
ChanNumsPreResponse = sum(AnatomicalLocsPatientsPreResponse, 2);

% this part is for the purpose of visulization. I'm going to show all the
% number of channels within each area of the brain. So I'm gonna sum all
% the values in each row of AnatomicalLocsPatientsAll.
ChanNumsAll = sum(AnatomicalLocsPatientsAll, 2);



% finding percentage
% in this part I need to find out the percentage of the number of channels
% in an specific area in the analysis, among all channels in that brain
% area. I want to use it in the bar plot, I want to sort the bar plot based
% on the percentage and not the number of channels

ChanNumsPercentPreOnset1 = nan( length(ChanNumsPreOnset1),1);
for i = 1:length(ChanNumsPreOnset1)
    ChanNumInAll = find(AnatomicalLocsVecAll == AnatomicalLocsVecPreOnset1(i)); 
    ratioNum = ChanNumsAll(ChanNumInAll);
    ratioNumPercent = floor((ChanNumsPreOnset1(i)/ratioNum)*100);
    disp(ratioNumPercent);
    ChanNumsPercentPreOnset1(i)=ratioNumPercent;
end
clear ChanNumInAll ratioNum ratioNumPercent



ChanNumsPercentPreOnset2 = nan( length(ChanNumsPreOnset2),1);
for i = 1:length(ChanNumsPreOnset2)
    ChanNumInAll = find(AnatomicalLocsVecAll == AnatomicalLocsVecPreOnset2(i)); 
    ratioNum = ChanNumsAll(ChanNumInAll);
    ratioNumPercent = floor((ChanNumsPreOnset2(i)/ratioNum)*100);
    ChanNumsPercentPreOnset2(i)=ratioNumPercent;
end
clear ChanNumInAll ratioNum ratioNumPercent


ChanNumsPercentPostOnset = nan( length(ChanNumsPostOnset),1);
for i = 1:length(ChanNumsPostOnset)
    ChanNumInAll = find(AnatomicalLocsVecAll == AnatomicalLocsVecPostOnset(i)); 
    ratioNum = ChanNumsAll(ChanNumInAll);
    ratioNumPercent = floor((ChanNumsPostOnset(i)/ratioNum)*100);
    ChanNumsPercentPostOnset(i)=ratioNumPercent;
end
clear ChanNumInAll ratioNum ratioNumPercent



ChanNumsPercentPreResponse = nan( length(ChanNumsPreResponse),1);
for i = 1:length(ChanNumsPreResponse)
    ChanNumInAll = find(AnatomicalLocsVecAll == AnatomicalLocsVecPreResponse(i)); 
    ratioNum = ChanNumsAll(ChanNumInAll);
    ratioNumPercent = floor((ChanNumsPreResponse(i)/ratioNum)*100);
    ChanNumsPercentPreResponse(i)=ratioNumPercent;
end
clear ChanNumInAll ratioNum ratioNumPercent


%% vis percentage

figure('Units', 'normalized', 'Position', [0.1, 0, 0.4, 2]);

position1 = [0.22, 0.99, 0.6, 0.1];
position2 = [0.22, 0.75, 0.6, 0.1];
position3 = [0.22, 0.5, 0.6, 0.1]; 
position4 = [0.22, 0.2, 0.6, 0.1]; 


threshold = 1; % Set to the minimum value that should be visible

% PreOnset1 plot
subplot('Position', position1);
visibleIndices = ChanNumsPreOnset1 > threshold & ChanNumsPercentPreOnset1 > threshold;

values = ChanNumsPercentPreOnset1(visibleIndices);
[sortedValues, sortOrder] = sort(values, 'descend');
sortedLabels = AnatomicalLocsVecPreOnset1(visibleIndices);
sortedLabels = sortedLabels(sortOrder);
sortedChans = ChanNumsPreOnset1(visibleIndices);
sortedChans = sortedChans(sortOrder);
bar(sortedValues, 0.5, 'FaceColor', [0.329, 0.329, 0.306]);
xticks(1:length(sortedValues));
xticklabels(sortedLabels);
set(gca, 'XTickLabel', get(gca, 'XTickLabel'),  'FontSize', 10, 'FontWeight', 'bold');
% set(gca, 'YTickLabel', get(gca, 'YTickLabel'),  'FontSize', 10, 'FontWeight', 'bold');
% ylim([0 max(sortedValues)+5]);
text(max(xlim)-0.03*max(xlim), max(ylim), 'PreOnset1', 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
for i = 1:length(sortedValues)
    ChanNumInAll = find(AnatomicalLocsVecAll == sortedLabels(i)); 
    NumOfAllChans = ChanNumsAll(ChanNumInAll);
    labelText = sprintf('%d%% \n %d/%d', sortedValues(i), sortedChans(i), NumOfAllChans);
    text(i, sortedValues(i) + 1, labelText, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
box off;




% PreOnset2 plot
subplot('Position', position2);
visibleIndices = ChanNumsPreOnset2 > threshold & ChanNumsPercentPreOnset2 > threshold;

values = ChanNumsPercentPreOnset2(visibleIndices);
[sortedValues, sortOrder] = sort(values, 'descend');
sortedLabels = AnatomicalLocsVecPreOnset2(visibleIndices);
sortedLabels = sortedLabels(sortOrder);
sortedChans = ChanNumsPreOnset2(visibleIndices);
sortedChans = sortedChans(sortOrder);
bar(sortedValues, 0.5, 'FaceColor', [0.329, 0.329, 0.306]);
xticks(1:length(sortedValues));
xticklabels(sortedLabels);
set(gca, 'XTickLabel', get(gca, 'XTickLabel'),  'FontSize', 10, 'FontWeight', 'bold');
% set(gca, 'YTickLabel', get(gca, 'YTickLabel'),  'FontSize', 10, 'FontWeight', 'bold');
% ylim([0 max(sortedValues)+5]);
text(max(xlim)-0.03*max(xlim), max(ylim), 'PreOnset2', 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
for i = 1:length(sortedValues)
    ChanNumInAll = find(AnatomicalLocsVecAll == sortedLabels(i)); 
    NumOfAllChans = ChanNumsAll(ChanNumInAll);
    labelText = sprintf('%d%% \n %d/%d', sortedValues(i), sortedChans(i), NumOfAllChans);
    text(i, sortedValues(i) + 1, labelText, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
box off;




% PostOnset plot
subplot('Position', position3);
visibleIndices = ChanNumsPostOnset > threshold & ChanNumsPercentPostOnset > threshold;

values = ChanNumsPercentPostOnset(visibleIndices);
[sortedValues, sortOrder] = sort(values, 'descend');
sortedLabels = AnatomicalLocsVecPostOnset(visibleIndices);
sortedLabels = sortedLabels(sortOrder);
sortedChans = ChanNumsPostOnset(visibleIndices);
sortedChans = sortedChans(sortOrder);
bar(sortedValues, 0.5, 'FaceColor', [0.329, 0.329, 0.306]);
xticks(1:length(sortedValues));
xticklabels(sortedLabels);
set(gca, 'XTickLabel', get(gca, 'XTickLabel'),  'FontSize', 10, 'FontWeight', 'bold');
% set(gca, 'YTickLabel', get(gca, 'YTickLabel'),  'FontSize', 10, 'FontWeight', 'bold');
% ylim([0 max(sortedValues)+5]);
text(max(xlim)-0.03*max(xlim), max(ylim), 'PostOnset', 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
for i = 1:length(sortedValues)
    ChanNumInAll = find(AnatomicalLocsVecAll == sortedLabels(i)); 
    NumOfAllChans = ChanNumsAll(ChanNumInAll);
    labelText = sprintf('%d%% \n %d/%d', sortedValues(i), sortedChans(i), NumOfAllChans);
    text(i, sortedValues(i) + 1, labelText, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
box off;




% PreResponse plot
subplot('Position', position4);
visibleIndices = ChanNumsPreResponse > threshold & ChanNumsPercentPreResponse > threshold;

values = ChanNumsPercentPreResponse(visibleIndices);
[sortedValues, sortOrder] = sort(values, 'descend');
sortedLabels = AnatomicalLocsVecPreResponse(visibleIndices);
sortedLabels = sortedLabels(sortOrder);
sortedChans = ChanNumsPreResponse(visibleIndices);
sortedChans = sortedChans(sortOrder);
bar(sortedValues, 0.5, 'FaceColor', [0.329, 0.329, 0.306]);
xticks(1:length(sortedValues));
xticklabels(sortedLabels);
set(gca, 'XTickLabel', get(gca, 'XTickLabel'),  'FontSize', 10, 'FontWeight', 'bold');
% set(gca, 'YTickLabel', get(gca, 'YTickLabel'),  'FontSize', 10, 'FontWeight', 'bold');
% ylim([0 max(sortedValues)+5]);
text(max(xlim)-0.03*max(xlim), max(ylim), 'PreResponse', 'FontWeight', 'bold', 'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
for i = 1:length(sortedValues)
    ChanNumInAll = find(AnatomicalLocsVecAll == sortedLabels(i)); 
    NumOfAllChans = ChanNumsAll(ChanNumInAll);
    labelText = sprintf('%d%% \n %d/%d', sortedValues(i), sortedChans(i), NumOfAllChans);
    text(i, sortedValues(i) + 1, labelText, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end
box off;




% Title and saving
suptitle(" ");
annotation('textbox', [0.1, 0.9, 0.9, 0.1], 'String', 'percentage of channels with significantly different IED and non-IED RTs', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18);
annotation('textbox', [0.12, 0.45, 0.3, 0.06], 'String', 'channels across patients (%)', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 18, 'Rotation', 90);

set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', [screenposition(3:4)]);
filename = 'RTs_percentage';
saveas(gcf, fullfile(outputFolderName, filename), 'pdf');

