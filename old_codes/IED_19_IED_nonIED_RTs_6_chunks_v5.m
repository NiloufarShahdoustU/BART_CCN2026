 % explained in notes

% Author: Nill



clear;
clc;
close all;
warning('off','all');


%%

inputFolderName_IEDdata = '\\155.100.91.44\d\Data\Nill\BART\bad_chans_removed_IEDdata_LFPmat_6_chunks_IED_11_v7';
outputFolderName = '\\155.100.91.44\d\Code\Nill\BART\IED\IED_19_IEDnonIED_bhvr_analysis\v5\';
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



for pt = 1:PatientsNum
% for pt = 6:6
    fileNameParts = strsplit(fileList(pt).name, '.');
    ptID = fileNameParts{1}; 
    disp("patient: " + ptID);

    %IED data read
    IEDdata = [inputFolderName_IEDdata '\' ptID '.LFPIED.mat'];
    load(IEDdata);

    nChans = length(LFPIED.selectedChans)-1; % removing the last selected chan

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

    allTimePoints = IEDtrialsPreOnset1 + IEDtrialsPreOnset2 + IEDtrialsPostOnset + IEDtrialsPreResponse;
    nonIEDIndices = find(all(allTimePoints == 0, 1));
   

        % now that I've found nonIED trials I need to have their behavioral
    % mearue vector: 
    % for the sake of each patient
    NonIED_bhvr_msr = RTs(nonIEDIndices);

    % for the sake of overall patients
    % I am doing standardized RTs:
    % the calculated means for each patient are influenced by the number of trials and the 
    % specific variability of that patient%s RTs. To address this, you can consider z-scoring
    % or standardizing the RTs for each patient before taking the mean. This ensures that differences
    % in baseline RTs and variability across patients are accounted for.
   


    % Step 2: Use the z-scored RTs to compute the mean for non-IED trials
    nonIEDtrials_bhvr_measure_mean(pt) = mean(NonIED_bhvr_msr)/mean(RTs); 





    % now I need to take those IEDtrials that ONLY happened in 1 of the 4
    % time periods and not in the other 3:
    IEDtrialsPreOnset1_only = IEDtrialsPreOnset1.*~IEDtrialsPreOnset2.*~IEDtrialsPostOnset.*~IEDtrialsPreResponse;
    IEDtrialsPreOnset2_only = IEDtrialsPreOnset2.*~IEDtrialsPreOnset1.*~IEDtrialsPostOnset.*~IEDtrialsPreResponse;
    IEDtrialsPostOnset_only = IEDtrialsPostOnset.*~IEDtrialsPreOnset1.*~IEDtrialsPreOnset2.*~IEDtrialsPreResponse;
    IEDtrialsPreResponse_only = IEDtrialsPreResponse.*~IEDtrialsPreOnset1.*~IEDtrialsPreOnset2.*~IEDtrialsPostOnset;

    


    NumberofPermutations = 10000;


    % stacked behavior measure for ied and non ied trials:
    
    %PreOnset1
    IED_StackedChans_bhvr_msr_PreOnset1 = [];


    %PreOnset2
    IED_StackedChans_bhvr_msr_PreOnset2 = [];

    %PostOnset
    IED_StackedChans_bhvr_msr_PostOnset = [];

    %PreResponse
    IED_StackedChans_bhvr_msr_PreResponse = [];

    

    for chz = 1:nChans
  
    
        %PreOnset1
        IEDTrials_bhvr_measure = IEDtrialsPreOnset1_only(chz,:); 
         % we are only taking the RTs that are not outliers
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*RTs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);



        
        if(~isempty(IEDTrials_bhvr_measure))

             IED_StackedChans_bhvr_msr_PreOnset1 = [IED_StackedChans_bhvr_msr_PreOnset1, IEDTrials_bhvr_measure]; 
        end

        %///////////////////////////////////////////////////////////////////////////////////////////////////////////

         %PreOnset2

  

        IEDTrials_bhvr_measure = IEDtrialsPreOnset2_only(chz,:); 
         % we are only taking the RTs that are not outliers
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*RTs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);



        if(~isempty(IEDTrials_bhvr_measure))

            IED_StackedChans_bhvr_msr_PreOnset2 = [IED_StackedChans_bhvr_msr_PreOnset2, IEDTrials_bhvr_measure]; 
        end
        %///////////////////////////////////////////////////////////////////////////////////////////////////////////

            %PostOnset


        IEDTrials_bhvr_measure = IEDtrialsPostOnset_only(chz,:); 
         % we are only taking the RTs that are not outliers
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*RTs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);



        if(~isempty(IEDTrials_bhvr_measure))

            IED_StackedChans_bhvr_msr_PostOnset = [IED_StackedChans_bhvr_msr_PostOnset, IEDTrials_bhvr_measure]; 
        end
        %///////////////////////////////////////////////////////////////////////////////////////////////////////////

            %PreResponse


        IEDTrials_bhvr_measure = IEDtrialsPreResponse_only(chz,:); 
         % we are only taking the RTs that are not outliers
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure.*RTs;
        IEDTrials_bhvr_measure = IEDTrials_bhvr_measure(IEDTrials_bhvr_measure ~= 0);



        if(~isempty(IEDTrials_bhvr_measure))

            IED_StackedChans_bhvr_msr_PreResponse = [IED_StackedChans_bhvr_msr_PreResponse, IEDTrials_bhvr_measure];
        end
        %///////////////////////////////////////////////////////////////////////////////////////////////////////////


   
  
    end % end for chan


    % permutation test for RTs:
    pValPreOnset1 = permutationTest(IED_StackedChans_bhvr_msr_PreOnset1,NonIED_bhvr_msr, NumberofPermutations);
    pValPreOnset2 = permutationTest(IED_StackedChans_bhvr_msr_PreOnset2,NonIED_bhvr_msr, NumberofPermutations);
    pValPostOnset = permutationTest(IED_StackedChans_bhvr_msr_PostOnset,NonIED_bhvr_msr, NumberofPermutations);
    pValPreResponse = permutationTest(IED_StackedChans_bhvr_msr_PreResponse,NonIED_bhvr_msr, NumberofPermutations);
   

    if pValPreOnset1 < 0.05
        IEDtrials_bhvr_measure_mean(pt,PreOnset1) = mean(IED_StackedChans_bhvr_msr_PreOnset1)/mean(RTs);
    end

    if pValPreOnset2 < 0.05
        IEDtrials_bhvr_measure_mean(pt,PreOnset2) = mean(IED_StackedChans_bhvr_msr_PreOnset2)/mean(RTs);
    end


    if pValPostOnset < 0.05
        IEDtrials_bhvr_measure_mean(pt,PostOnset) = mean(IED_StackedChans_bhvr_msr_PostOnset)/mean(RTs);
    end


    if pValPreResponse < 0.05
        IEDtrials_bhvr_measure_mean(pt,PreResponse) = mean(IED_StackedChans_bhvr_msr_PreResponse)/mean(RTs);
    end




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




%% visualization


% for RTs only the first 4 time periods, are important: 

vec1 = nonIEDtrials_bhvr_measure_mean;
vec2 = IEDtrials_bhvr_measure_mean_PreOnset1;
vec3 = IEDtrials_bhvr_measure_mean_PreOnset2;
vec4 = IEDtrials_bhvr_measure_mean_PostOnset;
vec5 = IEDtrials_bhvr_measure_mean_PreResponse;


% Concatenate all vectors into a single column vector and create group identifiers
allVecs = [vec1; vec2; vec3; vec4; vec5];
group = [ones(length(vec1),1); 2*ones(length(vec2),1); 3*ones(length(vec3),1); 4*ones(length(vec4),1); 5*ones(length(vec5),1) ];

figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.1, 0.25], 'Visible', 'on'); 

% Plot the boxplot with black color
boxplotHandle = boxplot(allVecs, group, 'Labels', {'non-IED', 'PreOnset1', 'PreOnset2', 'PostOnset', 'PreResponse'}, 'Color', 'k', 'Symbol', '');
xlabel('non-IED trials vs IED trials', 'FontSize', 14, 'FontWeight', 'bold');
% ylim([0 5]);
ylabel('normalized mean RT', 'FontSize', 14, 'FontWeight', 'bold');
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


%% debug
