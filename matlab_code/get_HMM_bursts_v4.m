%here, we will discretize the data we have collected and fit it to a hidden
%markov model, similar to Lammers/Galstyan.
function [outspatial, difftotvecbin2, bins] = get_HMM_bursts_v4(file, rep, bins, numbins)%, spatialAPbins)
%The first thing to do is to discretize the data into bins. There is a background
%fluorescence value for each experiment, which is the 0 point. But we only
%discretize the values past this threshold.

formatSpec = "G:\\My Drive\\t48_dynamics\\MANUSCRIPT\\code_for_sam_revision_hiddenmarkov_letsfinishthispaper\\%s\\trajectories_0%d.mat";
A1 = file;
A2 = rep;
str = sprintf(formatSpec,A1,A2);

load(str);

if strcmp(A1,'t48_control')
    out = 'control';
elseif strcmp(A1,'t48_dl1_best')
    out = 'dl1best';
elseif strcmp(A1,'t48_dl1_deletion')
    out = 'dl1del';
elseif strcmp(A1,'t48_dl2_best')
    out = 'dl2best';
elseif strcmp(A1,'t48_dl2_deletion')
    out = 'dl2del';
elseif strcmp(A1,'t48_zld_best')
    out = 'zldgood';
end
% out = regexp(A1,filesep,'split');
% out = out{1}(5:end);
namethresh = ['thresh_base_' out 'n' num2str(A2)]; %extract condition name to get threshold value automatically
out = strrep(out,'dl1','');
nameactive = ['activeM_' out 'n' num2str(A2)];
nameactive = strrep(nameactive,'dl1','');
namebin = ['loc_' out 'n' num2str(A2)];

threshold = load('G:\Shared drives\Lim_Lab\Sam\Affinity_Project\control_best_deletion_post_analysis.mat',namethresh);
active = load('G:\Shared drives\Lim_Lab\Sam\Affinity_Project\control_best_deletion_post_analysis.mat',nameactive);
spatialAPbins = load('G:\Shared drives\Lim_Lab\Sam\Affinity_Project\control_best_deletion_post_analysis.mat',namebin);

active = cell2mat(struct2cell(active));
threshold = 0.8*cell2mat(struct2cell(threshold));
spatialAPbins = cell2mat(struct2cell(spatialAPbins))*numbins; %16; changed to numbins


%we set a threshold of 2500. This is when we assume the promoter has
%reached steady state- once it turns on. 
%%
%now, we take the trajectories after they turn on for the first time, and
%we discretize their levels. 


% M = M(:,active);
for i = 1:length(active)
    M(:,i) = smooth(M(:,active(i)),4);
end
totvec = double((M(:)));
% totvec = totvec(totvec > 0.04); %this is hardcoded as the corresponding threshold to the 2500 uint16 value.

%%
%there are two options here- one is whether you want to read in bins from a
%previous run or want to generate them new for this run. This is useful
%because you might want to generate the bins for a control group and then
%use them as input for the next group. 

 %the number of bins is usually 8!
if isempty(bins) == 1 %signifying you dont start out with bins
totvec = double((M(:)));
bins = linspace(min(totvec), (max(totvec)), numbins); %makes the bins

totvecbin = zeros(size(totvec));

for i = 2:numbins
    for j = 1:length(totvec)
        if and(totvec(j) > bins(i-1), totvec(j) < bins(i))
            totvecbin(j) = i-1;
        end
    end
end

totvecbin = reshape(totvecbin,size(M));
end

if isempty(bins) == 0 %signifying you start out with bins
totvec = double((M(:)));

totvecbin = zeros(size(totvec));

for i = 2:length(bins)
    for j = 1:length(totvec)
        if and(totvec(j) > bins(i-1), totvec(j) < bins(i))
            totvecbin(j) = i-1;
        end
    end
end

totvecbin = reshape(totvecbin,size(M));
end

%%
%finally, we get a discretized set of states.

for i = 1:length(totvecbin(1, :))
    difftotvecbin(:, i) = diff(totvecbin(:, i));
end

difftotvecbintrans = difftotvecbin';

difftotvecbin2 = difftotvecbintrans - min(min(difftotvecbintrans)) +1 ; 

minemstate = min(min(difftotvecbintrans));
maxemstate = max(max(difftotvecbintrans));

  
trans = [0.5,0.5;
      0.50,0.50]; %hardcoded but it learns from these
  
negstates = abs(minemstate)+1;
posstates = abs(maxemstate);
  
emis1 = [1/(negstates).*ones(1, negstates) zeros(1, posstates)];
emis2 = [zeros(1, negstates) 1/(posstates).*ones(1, posstates)];

emis = [emis1;emis2];
%%
%ADDITION 2/24/2020: implementing spatial bins: there is probably an
%optimal way to do this. The way I settled upon is binning the embryo
%beforehand, using another script that is also in this folder.

%if isempty(spatialAPbins) == 0
    %what spatialAPbins is supposed to be is a discretized set of cells.
    outspatial = struct;
    
    % Sam's Method
    for i = 1:numbins %16
        tmp = find(spatialAPbins == i-1);
        tmp2 = difftotvecbin2(tmp,:);
        [estTRs,estEs] = hmmtrain(tmp2,trans,emis);
        outspatial(i).Pon = estTRs(1,2);
        outspatial(i).Poff = estTRs(2,1);
       
        outspatial(i).Emmis = estEs;
        outspatial(i).bintrajs = tmp2;
        tmp3 = tmp2-mode(tmp2);
        for q = 1:size(tmp3,1)
            for j = 1:size(tmp3,2)-1
                if tmp3(q,j) > 0 && tmp3(q,j+1) == 0
                    tmp3(q,j+1) = 1;
                elseif tmp3(q,j) < 0 && tmp3(q,j+1) <= 0
                    tmp3(q,j+1) = 0;
                    tmp3(q,j) = 0;
                elseif tmp3(q,j) < 0 && tmp3(q,j+1) > 0
                    tmp3(q,j) = 0;
                end
                if tmp3(q,j) > 1
                    tmp3(q,j) = 1;
                end
            end
            if tmp3(q,end) > 1
                tmp3(q,end) = 1;
            elseif tmp3(q,end) < 0
                tmp3(q,end) = 0;
            end
        end
        outspatial(i).promoter = tmp3;
        statesmat = [];
        for ii = 1:length(tmp2(:, 1))
        
        tempst = hmmviterbi(tmp2(ii, :),estTRs,estEs);
        statesmat = [statesmat; tempst];
        end
        outspatial(i).statesmat = statesmat;
        clear tmp tmp3
    end
end  
    
%   Sidu's method

%     for i = 1:max(spatialAPbins)
%         difftotvecbintemp = [];
%         for j = 1:length(difftotvecbin2)
%             if spatialAPbins(j) == i
%                 difftotvecbintemp = [difftotvecbintemp difftotvecbin2(j, :)];
%             end
%         end
%         [estTRs,estEs] = hmmtrain(difftotvecbintemp,trans,emis);
%         outspatial(i).Trans = estTRs;
%         outspatial(i).Emmis = estEs;
%     end
% end




% 
% %%
% if strcmp(version, 'bulk') == 1
%     [estTRs,estEs] = hmmtrain(difftotvecbin2,trans,emis);
% end
% 
% if strcmp(version, 'separate') == 1
%     estTRs = zeros(2, length(difftotvecbin2(:, 1)));
%     for i = 1:length(difftotvecbin2(:, 1))
%         if strcmp(state, 'all') == 1
%             [estTR,estE] = hmmtrain(difftotvecbin2(i, :),trans,emis);
%             estTRs(:, i) = [estTR(1, 2) estTR(2, 1)];
%         end
%         if strcmp(state, 'active') == 1
%             temp = difftotvecbin2(i, :);
%             idx = find(M>threshold); %supposed to determine when nuclei becomes active
%             sol = idx(1);
%             [estTR,estE] = hmmtrain(difftotvecbin2(i, sol:end),trans,emis);
%             estTRs(:, i) = [estTR(1, 2) estTR(2, 1)];
%         end
%         
%     
%     end
% end




