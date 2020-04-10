%Analysis with Normalized Time
clear;
load('G:\Shared drives\Lim_Lab\Sam\Affinity_Project\control_best_deletion_v3')

rep = 3;

active_Signal_c{1} = mean(M_controln1(:,activeM_controln1),1);
active_Signal_c{2} = mean(M_controln2(:,activeM_controln2),1);
active_Signal_c{3} = mean(M_controln4(:,activeM_controln4),1);
active_Signal_b{1} = mean(M_bestn1(:,activeM_bestn1),1);
active_Signal_b{2} = mean(M_bestn2(:,activeM_bestn2),1);
active_Signal_b{3} = mean(M_bestn4(:,activeM_bestn4),1);
active_Signal_d{1} = mean(M_deln1(:,activeM_deln1),1);
active_Signal_d{2} = mean(M_deln2(:,activeM_deln2),1);
active_Signal_d{3} = mean(M_deln4(:,activeM_deln4),1);
active_Signal_dl2b{1} = mean(M_dl2bestn1(:,activeM_dl2bestn1),1);
active_Signal_dl2b{2} = mean(M_dl2bestn2(:,activeM_dl2bestn2),1);
active_Signal_dl2b{3} = mean(M_dl2bestn5(:,activeM_dl2bestn5),1);
active_Signal_dl2d{1} = mean(M_dl2deln1(:,activeM_dl2deln1),1);
active_Signal_dl2d{2} = mean(M_dl2deln2(:,activeM_dl2deln2),1);
active_Signal_dl2d{3} = mean(M_dl2deln4(:,activeM_dl2deln4),1);
active_Signal_zldgood{1} = mean(M_zldgoodn2(:,activeM_zldgoodn2),1);
active_Signal_zldgood{2} = mean(M_zldgoodn3(:,activeM_zldgoodn3),1);
active_Signal_zldgood{3} = mean(M_zldgoodn4(:,activeM_zldgoodn4),1);
active_Signal_zldweak{1} = mean(M_zldweakn1(:,activeM_zldweakn1),1);
active_Signal_zldweak{2} = mean(M_zldweakn2(:,activeM_zldweakn2),1);
active_Signal_zldweak{3} = mean(M_zldweakn5(:,activeM_zldweakn5),1);

rep_Signal_c(1) = mean(active_Signal_c{1});
rep_Signal_c(2) = mean(active_Signal_c{2});
rep_Signal_c(3) = mean(active_Signal_c{3});
rep_Signal_b(1) = mean(active_Signal_b{1});
rep_Signal_b(2) = mean(active_Signal_b{2});
rep_Signal_b(3) = mean(active_Signal_b{3});
rep_Signal_d(1) = mean(active_Signal_d{1});
rep_Signal_d(2) = mean(active_Signal_d{2});
rep_Signal_d(3) = mean(active_Signal_d{3});
rep_Signal_dl2b(1) = mean(active_Signal_dl2b{1});
rep_Signal_dl2b(2) = mean(active_Signal_dl2b{2});
rep_Signal_dl2b(3) = mean(active_Signal_dl2b{3});
rep_Signal_dl2d(1) = mean(active_Signal_dl2d{1});
rep_Signal_dl2d(2) = mean(active_Signal_dl2d{2});
rep_Signal_dl2d(3) = mean(active_Signal_dl2d{3});
rep_Signal_zldgood(1) = mean(active_Signal_zldgood{1});
rep_Signal_zldgood(2) = mean(active_Signal_zldgood{2});
rep_Signal_zldgood(3) = mean(active_Signal_zldgood{3});
rep_Signal_zldweak(1) = mean(active_Signal_zldweak{1});
rep_Signal_zldweak(2) = mean(active_Signal_zldweak{2});
rep_Signal_zldweak(3) = mean(active_Signal_zldweak{3});

% stdrep_Signal_c = std(active_Signal_c);
% stdrep_Signal_b = std(active_Signal_b);
% stdrep_Signal_d = std(active_Signal_d);
% stdrep_Signal_dl2b = std(active_Signal_dl2b);
% stdrep_Signal_dl2d = std(active_Signal_dl2d);
% stdrep_Signal_zldgood = std(active_Signal_zldgood);
% stdrep_Signal_zldweak = std(active_Signal_zldweak);

figure;
for i = 1:rep
    subplot(rep,1,i)
    histogram(active_Signal_zldweak{i},25);
    xlim([0 4500])
    title(['Zld Weak Replicate ' num2str(i)])
end

figure;
tmp = active_Signal_zldweak;
c = {'Rep 1', 'Rep 2', 'Rep 3'};
g = [zeros(1,length(tmp{1})) ones(1,length(tmp{2})) 2*ones(1,length(tmp{3}))];
boxplot([tmp{1} tmp{2} tmp{3}],g,'Labels',c)
ylabel('Signal');
ylim([0 4500])
title('Nuclei Signal - Zld Weak')





norm_time = .1:.005:1;

time_controln1 = T_controln1/T_controln1(end);
time_controln2 = T_controln2/T_controln2(end);
time_controln4 = T_controln4/T_controln4(end);
time_controln3 = T_controln3/T_controln3(end);
time_controln5 = T_controln5/T_controln5(end);
time_controln6 = T_controln6/T_controln6(end);
time_bestn1 = T_bestn1/T_bestn1(end);
time_bestn2 = T_bestn2/T_bestn2(end);
time_bestn4 = T_bestn4/T_bestn4(end);
time_deln1 = T_deln1/T_deln1(end);
time_deln2 = T_deln2/T_deln2(end);
time_deln4 = T_deln4/T_deln4(end);
time_dl2bestn1 = T_dl2bestn1/T_dl2bestn1(end);
time_dl2bestn2 = T_dl2bestn2/T_dl2bestn2(end);
time_dl2bestn5 = T_dl2bestn5/T_dl2bestn5(end);
time_dl2deln1 = T_dl2deln1/T_dl2deln1(end);
time_dl2deln2 = T_dl2deln2/T_dl2deln2(end);
time_dl2deln4 = T_dl2deln4/T_dl2deln4(end);
time_zldgoodn2 = T_zldgoodn2/T_zldgoodn2(end);
time_zldgoodn3 = T_zldgoodn3/T_zldgoodn3(end);
time_zldgoodn4 = T_zldgoodn4/T_zldgoodn4(end);
time_zldweakn1 = T_zldweakn1/T_zldweakn1(end);
time_zldweakn2 = T_zldweakn2/T_zldweakn2(end);
time_zldweakn5 = T_zldweakn5/T_zldweakn5(end);
time_mcphetn1 = T_mcphetn1/T_mcphetn1(end);
time_mcphetn2 = T_mcphetn2/T_mcphetn2(end);
time_mcphetn3 = T_mcphetn3/T_mcphetn3(end);
time_dlhetcontroln1 = T_dlhetcontroln1/T_dlhetcontroln1(end);
time_dlhetcontroln3 = T_dlhetcontroln3/T_dlhetcontroln3(end);
time_dlhetcontroln4 = T_dlhetcontroln4/T_dlhetcontroln4(end);
time_dlhetdl2bestn2 = T_dlhetdl2bestn2/T_dlhetdl2bestn2(end);
time_dlhetdl2bestn3 = T_dlhetdl2bestn3/T_dlhetdl2bestn3(end);
time_dlhetdl2bestn5 = T_dlhetdl2bestn5/T_dlhetdl2bestn5(end);
time_twideln1 = T_twideln1/T_twideln1(end);
time_twideln2 = T_twideln2/T_twideln2(end);

%het
time_hetcontroln1 = T_hetcontroln1/T_hetcontroln1(end);
time_hetcontroln2 = T_hetcontroln2/T_hetcontroln2(end);
time_hetcontroln3 = T_hetcontroln3/T_hetcontroln3(end);
time_hetbestn1 = T_hetbestn1/T_hetbestn1(end);
time_hetbestn2 = T_hetbestn2/T_hetbestn2(end);
time_hetbestn3 = T_hetbestn3/T_hetbestn3(end);
time_hetdl2bestn1 = T_hetdl2bestn1/T_hetdl2bestn1(end);
time_hetdl2bestn2 = T_hetdl2bestn2/T_hetdl2bestn2(end);
time_hetdl2bestn3 = T_hetdl2bestn3/T_hetdl2bestn3(end);
time_hetdl2deln1 = T_hetdl2deln1/T_hetdl2deln1(end);
time_hetdl2deln2 = T_hetdl2deln2/T_hetdl2deln2(end);
time_hetdl2deln3 = T_hetdl2deln3/T_hetdl2deln3(end);
time_hetzldbestn1 = T_hetzldbestn1/T_hetzldbestn1(end);
time_hetzldbestn2 = T_hetzldbestn2/T_hetzldbestn2(end);
time_hetzldbestn3 = T_hetzldbestn3/T_hetzldbestn3(end);

%Revision
time_Rcontroln1 = T_Rcn1/T_Rcn1(end);
time_Rcontroln2 = T_Rcn2/T_Rcn2(end);
time_Rcontroln3 = T_Rcn3/T_Rcn3(end);

time_Rdzn1 = T_Rdzn1/T_Rdzn1(end);
time_Rdzn2 = T_Rdzn2/T_Rdzn2(end);
time_Rdzn3 = T_Rdzn3/T_Rdzn3(end);

time_Rdn1 = T_Rdn1/T_Rdn1(end);
time_Rdn2 = T_Rdn2/T_Rdn2(end);
time_Rdn3 = T_Rdn3/T_Rdn3(end);
%% Activation Over Time
%Number active nuclei over time
nuclei_by_T_c(1,:) = interp1(time_controln1,nuclei_by_T_controln1,norm_time);
nuclei_by_T_c(2,:) = interp1(time_controln2,nuclei_by_T_controln2,norm_time);
nuclei_by_T_c(3,:) = interp1(time_controln4,nuclei_by_T_controln4,norm_time);
nuclei_by_T_c(4,:) = interp1(time_controln3,nuclei_by_T_controln3,norm_time);
nuclei_by_T_c(5,:) = interp1(time_controln5,nuclei_by_T_controln5,norm_time);
nuclei_by_T_c(6,:) = interp1(time_controln6,nuclei_by_T_controln6,norm_time);

nuclei_by_T_b(1,:) = interp1(time_bestn1,nuclei_by_T_bestn1,norm_time);
nuclei_by_T_b(2,:) = interp1(time_bestn2,nuclei_by_T_bestn2,norm_time);
nuclei_by_T_b(3,:) = interp1(time_bestn4,nuclei_by_T_bestn4,norm_time);

nuclei_by_T_d(1,:) = interp1(time_deln1,nuclei_by_T_deln1,norm_time);
nuclei_by_T_d(2,:) = interp1(time_deln2,nuclei_by_T_deln2,norm_time);
nuclei_by_T_d(3,:) = interp1(time_deln4,nuclei_by_T_deln4,norm_time);

nuclei_by_T_dl2b(1,:) = interp1(time_dl2bestn1,nuclei_by_T_dl2bestn1,norm_time);
nuclei_by_T_dl2b(2,:) = interp1(time_dl2bestn2,nuclei_by_T_dl2bestn2,norm_time);
nuclei_by_T_dl2b(3,:) = interp1(time_dl2bestn5,nuclei_by_T_dl2bestn5,norm_time);

nuclei_by_T_dl2d(1,:) = interp1(time_dl2deln1,nuclei_by_T_dl2deln1,norm_time);
nuclei_by_T_dl2d(2,:) = interp1(time_dl2deln2,nuclei_by_T_dl2deln2,norm_time);
nuclei_by_T_dl2d(3,:) = interp1(time_dl2deln4,nuclei_by_T_dl2deln4,norm_time);

nuclei_by_T_zg(1,:) = interp1(time_zldgoodn2,nuclei_by_T_zldgoodn2,norm_time);
nuclei_by_T_zg(2,:) = interp1(time_zldgoodn3,nuclei_by_T_zldgoodn3,norm_time);
nuclei_by_T_zg(3,:) = interp1(time_zldgoodn4,nuclei_by_T_zldgoodn4,norm_time);

nuclei_by_T_zw(1,:) = interp1(time_zldweakn1,nuclei_by_T_zldweakn1,norm_time);
nuclei_by_T_zw(2,:) = interp1(time_zldweakn2,nuclei_by_T_zldweakn2,norm_time);
nuclei_by_T_zw(3,:) = interp1(time_zldweakn5,nuclei_by_T_zldweakn5,norm_time);

nuclei_by_T_mcphet(1,:) = interp1(time_mcphetn1,nuclei_by_T_mcphetn1,norm_time);
nuclei_by_T_mcphet(2,:) = interp1(time_mcphetn2,nuclei_by_T_mcphetn2,norm_time);
nuclei_by_T_mcphet(3,:) = interp1(time_mcphetn3,nuclei_by_T_mcphetn3,norm_time);

nuclei_by_T_dlhetcontrol(1,:) = interp1(time_dlhetcontroln1,nuclei_by_T_dlhetcontroln1,norm_time);
nuclei_by_T_dlhetcontrol(2,:) = interp1(time_dlhetcontroln3,nuclei_by_T_dlhetcontroln3,norm_time);
nuclei_by_T_dlhetcontrol(3,:) = interp1(time_dlhetcontroln4,nuclei_by_T_dlhetcontroln4,norm_time);

nuclei_by_T_dlhetdl2best(1,:) = interp1(time_dlhetdl2bestn2,nuclei_by_T_dlhetdl2bestn2,norm_time);
nuclei_by_T_dlhetdl2best(2,:) = interp1(time_dlhetdl2bestn3,nuclei_by_T_dlhetdl2bestn3,norm_time);
nuclei_by_T_dlhetdl2best(3,:) = interp1(time_dlhetdl2bestn5,nuclei_by_T_dlhetdl2bestn5,norm_time);


ave_nuc_T_c = mean(nuclei_by_T_c,1);
std_nuc_T_c = std(nuclei_by_T_c,1);

ave_nuc_T_b = mean(nuclei_by_T_b,1);
std_nuc_T_b = std(nuclei_by_T_b,1);

ave_nuc_T_d = mean(nuclei_by_T_d,1);
std_nuc_T_d = std(nuclei_by_T_d,1);

ave_nuc_T_dl2b = mean(nuclei_by_T_dl2b,1);
std_nuc_T_dl2b = std(nuclei_by_T_dl2b,1);

ave_nuc_T_dl2d = mean(nuclei_by_T_dl2d,1);
std_nuc_T_dl2d = std(nuclei_by_T_dl2d,1);

ave_nuc_T_zg = mean(nuclei_by_T_zg,1);
std_nuc_T_zg = std(nuclei_by_T_zg,1);

ave_nuc_T_zw = mean(nuclei_by_T_zw,1);
std_nuc_T_zw = std(nuclei_by_T_zw,1);

ave_nuc_T_mcphet = mean(nuclei_by_T_mcphet,1);
std_nuc_T_mcphet = std(nuclei_by_T_mcphet,1);

ave_nuc_T_dlhetcontrol = mean(nuclei_by_T_dlhetcontrol,1);
std_nuc_T_dlhetcontrol = std(nuclei_by_T_dlhetcontrol,1);

ave_nuc_T_dlhetdl2best = mean(nuclei_by_T_dlhetdl2best,1);
std_nuc_T_dlhetdl2best = std(nuclei_by_T_dlhetdl2best,1);

%het 
nuclei_by_T_hetc(1,:) = interp1(time_hetcontroln1,nuclei_by_T_hetcontroln1,norm_time);
nuclei_by_T_hetc(2,:) = interp1(time_hetcontroln2,nuclei_by_T_hetcontroln2,norm_time);
nuclei_by_T_hetc(3,:) = interp1(time_hetcontroln3,nuclei_by_T_hetcontroln3,norm_time);

nuclei_by_T_hetb(1,:) = interp1(time_hetbestn1,nuclei_by_T_hetbestn1,norm_time);
nuclei_by_T_hetb(2,:) = interp1(time_hetbestn2,nuclei_by_T_hetbestn2,norm_time);
nuclei_by_T_hetb(3,:) = interp1(time_hetbestn3,nuclei_by_T_hetbestn3,norm_time);

nuclei_by_T_hetdl2b(1,:) = interp1(time_hetdl2bestn1,nuclei_by_T_hetdl2bestn1,norm_time);
nuclei_by_T_hetdl2b(2,:) = interp1(time_hetdl2bestn2,nuclei_by_T_hetdl2bestn2,norm_time);
nuclei_by_T_hetdl2b(3,:) = interp1(time_hetdl2bestn3,nuclei_by_T_hetdl2bestn3,norm_time);

nuclei_by_T_hetdl2d(1,:) = interp1(time_hetdl2deln1,nuclei_by_T_hetdl2deln1,norm_time);
nuclei_by_T_hetdl2d(2,:) = interp1(time_hetdl2deln2,nuclei_by_T_hetdl2deln2,norm_time);
nuclei_by_T_hetdl2d(3,:) = interp1(time_hetdl2deln3,nuclei_by_T_hetdl2deln3,norm_time);

nuclei_by_T_hetzg(1,:) = interp1(time_hetzldbestn1,nuclei_by_T_hetzldbestn1,norm_time);
nuclei_by_T_hetzg(2,:) = interp1(time_hetzldbestn2,nuclei_by_T_hetzldbestn2,norm_time);
nuclei_by_T_hetzg(3,:) = interp1(time_hetzldbestn3,nuclei_by_T_hetzldbestn3,norm_time);

ave_nuc_T_hetc = mean(nuclei_by_T_hetc,1);
std_nuc_T_hetc = std(nuclei_by_T_hetc,1);

ave_nuc_T_hetb = mean(nuclei_by_T_hetb,1);
std_nuc_T_hetb = std(nuclei_by_T_hetb,1);

ave_nuc_T_hetdl2b = mean(nuclei_by_T_hetdl2b,1);
std_nuc_T_hetdl2b = std(nuclei_by_T_hetdl2b,1);

ave_nuc_T_hetdl2d = mean(nuclei_by_T_hetdl2d,1);
std_nuc_T_hetdl2d = std(nuclei_by_T_hetdl2d,1);

ave_nuc_T_hetzg = mean(nuclei_by_T_hetzg,1);
std_nuc_T_hetzg = std(nuclei_by_T_hetzg,1);

%Revision
nuclei_by_T_Rc(1,:) = interp1(time_Rcontroln1,nuclei_by_T_Rcn1,norm_time);
nuclei_by_T_Rc(2,:) = interp1(time_Rcontroln2,nuclei_by_T_Rcn2,norm_time);
nuclei_by_T_Rc(3,:) = interp1(time_Rcontroln3,nuclei_by_T_Rcn3,norm_time);

nuclei_by_T_Rdz(1,:) = interp1(time_Rdzn1,nuclei_by_T_Rdzn1,norm_time);
nuclei_by_T_Rdz(2,:) = interp1(time_Rdzn2,nuclei_by_T_Rdzn2,norm_time);
nuclei_by_T_Rdz(3,:) = interp1(time_Rdzn3,nuclei_by_T_Rdzn3,norm_time);

nuclei_by_T_Rd(1,:) = interp1(time_Rdn1,nuclei_by_T_Rdn1,norm_time);
nuclei_by_T_Rd(2,:) = interp1(time_Rdn2,nuclei_by_T_Rdn2,norm_time);
nuclei_by_T_Rd(3,:) = interp1(time_Rdn3,nuclei_by_T_Rdn3,norm_time);

ave_nuc_T_Rc = mean(nuclei_by_T_Rc,1);
std_nuc_T_Rc = std(nuclei_by_T_Rc,1);
ave_nuc_T_Rdz = mean(nuclei_by_T_Rdz,1);
std_nuc_T_Rdz = std(nuclei_by_T_Rdz,1);

%Normalized active nuclei over time
act_by_T_c(1,:) = interp1(T_controln1/T_controln1(end),nuclei_by_T_controln1./nuclei_by_T_controln1(end),norm_time);
act_by_T_c(2,:) = interp1(T_controln2/T_controln2(end),nuclei_by_T_controln2./nuclei_by_T_controln2(end),norm_time);
act_by_T_c(3,:) = interp1(T_controln4/T_controln4(end),nuclei_by_T_controln4./nuclei_by_T_controln4(end),norm_time);
act_by_T_c(4,:) = interp1(T_controln3/T_controln3(end),nuclei_by_T_controln3./nuclei_by_T_controln3(end),norm_time);
act_by_T_c(5,:) = interp1(T_controln5/T_controln5(end),nuclei_by_T_controln5./nuclei_by_T_controln5(end),norm_time);
act_by_T_c(6,:) = interp1(T_controln6/T_controln6(end),nuclei_by_T_controln6./nuclei_by_T_controln6(end),norm_time);

act_by_T_b(1,:) = interp1(T_bestn1/T_bestn1(end),nuclei_by_T_bestn1./nuclei_by_T_bestn1(end),norm_time);
act_by_T_b(2,:) = interp1(T_bestn2/T_bestn2(end),nuclei_by_T_bestn2./nuclei_by_T_bestn2(end),norm_time);
act_by_T_b(3,:) = interp1(T_bestn4/T_bestn4(end),nuclei_by_T_bestn4./nuclei_by_T_bestn4(end),norm_time);

act_by_T_d(1,:) = interp1(T_deln1/T_deln1(end),nuclei_by_T_deln1./nuclei_by_T_deln1(end),norm_time);
act_by_T_d(2,:) = interp1(T_deln2/T_deln2(end),nuclei_by_T_deln2./nuclei_by_T_deln2(end),norm_time);
act_by_T_d(3,:) = interp1(T_deln4/T_deln4(end),nuclei_by_T_deln4./nuclei_by_T_deln4(end),norm_time);

act_by_T_dl2b(1,:) = interp1(T_dl2bestn1/T_dl2bestn1(end),nuclei_by_T_dl2bestn1./nuclei_by_T_dl2bestn1(end),norm_time);
act_by_T_dl2b(2,:) = interp1(T_dl2bestn2/T_dl2bestn2(end),nuclei_by_T_dl2bestn2./nuclei_by_T_dl2bestn2(end),norm_time);
act_by_T_dl2b(3,:) = interp1(T_dl2bestn5/T_dl2bestn5(end),nuclei_by_T_dl2bestn5./nuclei_by_T_dl2bestn5(end),norm_time);

act_by_T_dl2d(1,:) = interp1(T_dl2deln1/T_dl2deln1(end),nuclei_by_T_dl2deln1./nuclei_by_T_dl2deln1(end),norm_time);
act_by_T_dl2d(2,:) = interp1(T_dl2deln2/T_dl2deln2(end),nuclei_by_T_dl2deln2./nuclei_by_T_dl2deln2(end),norm_time);
act_by_T_dl2d(3,:) = interp1(T_dl2deln4/T_dl2deln4(end),nuclei_by_T_dl2deln4./nuclei_by_T_dl2deln4(end),norm_time);

act_by_T_zg(1,:) = interp1(T_zldgoodn2/T_zldgoodn2(end),nuclei_by_T_zldgoodn2./nuclei_by_T_zldgoodn2(end),norm_time);
act_by_T_zg(2,:) = interp1(T_zldgoodn3/T_zldgoodn3(end),nuclei_by_T_zldgoodn3./nuclei_by_T_zldgoodn3(end),norm_time);
act_by_T_zg(3,:) = interp1(T_zldgoodn4/T_zldgoodn4(end),nuclei_by_T_zldgoodn4./nuclei_by_T_zldgoodn4(end),norm_time);

act_by_T_zw(1,:) = interp1(T_zldweakn1/T_zldweakn1(end),nuclei_by_T_zldweakn1./nuclei_by_T_zldweakn1(end),norm_time);
act_by_T_zw(2,:) = interp1(T_zldweakn2/T_zldweakn2(end),nuclei_by_T_zldweakn2./nuclei_by_T_zldweakn2(end),norm_time);
act_by_T_zw(3,:) = interp1(T_zldweakn5/T_zldweakn5(end),nuclei_by_T_zldweakn5./nuclei_by_T_zldweakn5(end),norm_time);

act_by_T_mcphet(1,:) = interp1(T_mcphetn1/T_mcphetn1(end),nuclei_by_T_mcphetn1./nuclei_by_T_mcphetn1(end),norm_time);
act_by_T_mcphet(2,:) = interp1(T_mcphetn2/T_mcphetn2(end),nuclei_by_T_mcphetn2./nuclei_by_T_mcphetn2(end),norm_time);
act_by_T_mcphet(3,:) = interp1(T_mcphetn3/T_mcphetn3(end),nuclei_by_T_mcphetn3./nuclei_by_T_mcphetn3(end),norm_time);

act_by_T_dlhetcontrol(1,:) = interp1(T_dlhetcontroln1/T_dlhetcontroln1(end),nuclei_by_T_dlhetcontroln1./nuclei_by_T_dlhetcontroln1(end),norm_time);
act_by_T_dlhetcontrol(2,:) = interp1(T_dlhetcontroln3/T_dlhetcontroln3(end),nuclei_by_T_dlhetcontroln3./nuclei_by_T_dlhetcontroln3(end),norm_time);
act_by_T_dlhetcontrol(3,:) = interp1(T_dlhetcontroln4/T_dlhetcontroln4(end),nuclei_by_T_dlhetcontroln4./nuclei_by_T_dlhetcontroln4(end),norm_time);

act_by_T_dlhetdl2best(1,:) = interp1(T_dlhetdl2bestn2/T_dlhetdl2bestn2(end),nuclei_by_T_dlhetdl2bestn2./nuclei_by_T_dlhetdl2bestn2(end),norm_time);
act_by_T_dlhetdl2best(2,:) = interp1(T_dlhetdl2bestn3/T_dlhetdl2bestn3(end),nuclei_by_T_dlhetdl2bestn3./nuclei_by_T_dlhetdl2bestn3(end),norm_time);
act_by_T_dlhetdl2best(3,:) = interp1(T_dlhetdl2bestn5/T_dlhetdl2bestn5(end),nuclei_by_T_dlhetdl2bestn5./nuclei_by_T_dlhetdl2bestn5(end),norm_time);

act_by_T_td(1,:) = interp1(T_twideln1/T_twideln1(end),nuclei_by_T_twideln1./nuclei_by_T_twideln1(end),norm_time);
act_by_T_td(2,:) = interp1(T_twideln2/T_twideln2(end),nuclei_by_T_twideln2./nuclei_by_T_twideln2(end),norm_time);

ave_act_T_c = mean(act_by_T_c,1);
std_act_T_c = std(act_by_T_c,1);

ave_act_T_b = mean(act_by_T_b,1);
std_act_T_b = std(act_by_T_b,1);

ave_act_T_d = mean(act_by_T_d,1);
std_act_T_d = std(act_by_T_d,1);

ave_act_T_dl2b = mean(act_by_T_dl2b,1);
std_act_T_dl2b = std(act_by_T_dl2b,1);

ave_act_T_dl2d = mean(act_by_T_dl2d,1);
std_act_T_dl2d = std(act_by_T_dl2d,1);

ave_act_T_zg = mean(act_by_T_zg,1);
std_act_T_zg = std(act_by_T_zg,1);

ave_act_T_td = mean(act_by_T_td,1);
std_act_T_td = std(act_by_T_td,1);

ave_act_T_zw = mean(act_by_T_zw,1);
std_act_T_zw = std(act_by_T_zw,1);

ave_act_T_mcphet = mean(act_by_T_mcphet,1);
std_act_T_mcphet = std(act_by_T_mcphet,1);

ave_act_T_dlhetcontrol = mean(act_by_T_dlhetcontrol,1);
std_act_T_dlhetcontrol = std(act_by_T_dlhetcontrol,1);

ave_act_T_dlhetdl2best = mean(act_by_T_dlhetdl2best,1);
std_act_T_dlhetdl2best = std(act_by_T_dlhetdl2best,1);


%het
act_by_T_hetc(1,:) = interp1(T_hetcontroln1/T_hetcontroln1(end),nuclei_by_T_hetcontroln1./nuclei_by_T_hetcontroln1(end),norm_time);
act_by_T_hetc(2,:) = interp1(T_hetcontroln2/T_hetcontroln2(end),nuclei_by_T_hetcontroln2./nuclei_by_T_hetcontroln2(end),norm_time);
act_by_T_hetc(3,:) = interp1(T_hetcontroln3/T_hetcontroln3(end),nuclei_by_T_hetcontroln3./nuclei_by_T_hetcontroln3(end),norm_time);

act_by_T_hetb(1,:) = interp1(T_hetbestn1/T_hetbestn1(end),nuclei_by_T_hetbestn1./nuclei_by_T_hetbestn1(end),norm_time);
act_by_T_hetb(2,:) = interp1(T_hetbestn2/T_hetbestn2(end),nuclei_by_T_hetbestn2./nuclei_by_T_hetbestn2(end),norm_time);
act_by_T_hetb(3,:) = interp1(T_hetbestn3/T_hetbestn3(end),nuclei_by_T_hetbestn3./nuclei_by_T_hetbestn3(end),norm_time);

act_by_T_hetdl2b(1,:) = interp1(T_hetdl2bestn1/T_hetdl2bestn1(end),nuclei_by_T_hetdl2bestn1./nuclei_by_T_hetdl2bestn1(end),norm_time);
act_by_T_hetdl2b(2,:) = interp1(T_hetdl2bestn2/T_hetdl2bestn2(end),nuclei_by_T_hetdl2bestn2./nuclei_by_T_hetdl2bestn2(end),norm_time);
act_by_T_hetdl2b(3,:) = interp1(T_hetdl2bestn3/T_hetdl2bestn3(end),nuclei_by_T_hetdl2bestn3./nuclei_by_T_hetdl2bestn3(end),norm_time);

act_by_T_hetdl2d(1,:) = interp1(T_hetdl2deln1/T_hetdl2deln1(end),nuclei_by_T_hetdl2deln1./nuclei_by_T_hetdl2deln1(end),norm_time);
act_by_T_hetdl2d(2,:) = interp1(T_hetdl2deln2/T_hetdl2deln2(end),nuclei_by_T_hetdl2deln2./nuclei_by_T_hetdl2deln2(end),norm_time);
act_by_T_hetdl2d(3,:) = interp1(T_hetdl2deln3/T_hetdl2deln3(end),nuclei_by_T_hetdl2deln3./nuclei_by_T_hetdl2deln3(end),norm_time);

act_by_T_hetzg(1,:) = interp1(T_hetzldbestn1/T_hetzldbestn1(end),nuclei_by_T_hetzldbestn1./nuclei_by_T_hetzldbestn1(end),norm_time);
act_by_T_hetzg(2,:) = interp1(T_hetzldbestn2/T_hetzldbestn2(end),nuclei_by_T_hetzldbestn2./nuclei_by_T_hetzldbestn2(end),norm_time);
act_by_T_hetzg(3,:) = interp1(T_hetzldbestn3/T_hetzldbestn3(end),nuclei_by_T_hetzldbestn3./nuclei_by_T_hetzldbestn3(end),norm_time);

ave_act_T_hetc = mean(act_by_T_hetc,1);
std_act_T_hetc = std(act_by_T_hetc,1);

ave_act_T_hetb = mean(act_by_T_hetb,1);
std_act_T_hetb = std(act_by_T_hetb,1);

ave_act_T_hetdl2b = mean(act_by_T_hetdl2b,1);
std_act_T_hetdl2b = std(act_by_T_hetdl2b,1);

ave_act_T_hetdl2d = mean(act_by_T_hetdl2d,1);
std_act_T_hetdl2d = std(act_by_T_hetdl2d,1);

ave_act_T_hetzg = mean(act_by_T_hetzg,1);
std_act_T_hetzg = std(act_by_T_hetzg,1);

%Revision
act_by_T_Rc(1,:) = interp1(T_Rcn1/T_Rcn1(end),nuclei_by_T_Rcn1./nuclei_by_T_Rcn1(end),norm_time);
act_by_T_Rc(2,:) = interp1(T_Rcn2/T_Rcn2(end),nuclei_by_T_Rcn2./nuclei_by_T_Rcn2(end),norm_time);
act_by_T_Rc(3,:) = interp1(T_Rcn3/T_Rcn3(end),nuclei_by_T_Rcn3./nuclei_by_T_Rcn3(end),norm_time);
act_by_T_Rdz(1,:) = interp1(T_Rdzn1/T_Rdzn1(end),nuclei_by_T_Rdzn1./nuclei_by_T_Rdzn1(end),norm_time);
act_by_T_Rdz(2,:) = interp1(T_Rdzn2/T_Rdzn2(end),nuclei_by_T_Rdzn2./nuclei_by_T_Rdzn2(end),norm_time);
act_by_T_Rdz(3,:) = interp1(T_Rdzn3/T_Rdzn3(end),nuclei_by_T_Rdzn3./nuclei_by_T_Rdzn3(end),norm_time);
act_by_T_Rd(1,:) = interp1(T_Rdn1/T_Rdn1(end),nuclei_by_T_Rdn1./nuclei_by_T_Rdn1(end),norm_time);
act_by_T_Rd(2,:) = interp1(T_Rdn2/T_Rdn2(end),nuclei_by_T_Rdn2./nuclei_by_T_Rdn2(end),norm_time);
act_by_T_Rd(3,:) = interp1(T_Rdn3/T_Rdn3(end),nuclei_by_T_Rdn3./nuclei_by_T_Rdn3(end),norm_time);

ave_act_T_Rc = mean(act_by_T_Rc,1);
std_act_T_Rc = std(act_by_T_Rc,1);
ave_act_T_Rdz = mean(act_by_T_Rdz,1);
std_act_T_Rdz = std(act_by_T_Rdz,1);
ave_act_T_Rd = mean(act_by_T_Rd,1);
std_act_T_Rd = std(act_by_T_Rd,1);


figure;
plot(norm_time,nuclei_by_T_c(1,:)); hold on
plot(norm_time,nuclei_by_T_c(2,:))
plot(norm_time,nuclei_by_T_c(3,:))

%{
%Poster Plots
figure;
shadedErrorBar(norm_time,ave_act_T_c,std_act_T_c,{'color',rgb('black'),'LineWidth',2},1); hold on
shadedErrorBar(norm_time,ave_act_T_b,std_act_T_b,{'color',rgb('red'),'LineWidth',2},1); hold on
shadedErrorBar(norm_time,ave_act_T_d,std_act_T_d,{'--','color',rgb('lightcoral'),'LineWidth',2},1); hold on
title('Active Nuclei vs Time')
ylabel('Percent Active Nuclei')
xlabel('t''')
ylim([0 1])
legend('control', 'Dl1 strong','Dl1 weak','Location', 'northwest')

shadedErrorBar(norm_time,ave_act_T_c,std_act_T_c,{'color',rgb('black'),'LineWidth',2},1); hold on
shadedErrorBar(norm_time,ave_act_T_dl2b,std_act_T_dl2b,{'color',rgb('blue'),'LineWidth',2},1); hold on 
shadedErrorBar(norm_time,ave_act_T_dl2d,std_act_T_dl2d,{'--','color',rgb('cornflowerblue'),'LineWidth',2},1);
title('Number of Active Nuclei')
ylabel('Percent Active Nuclei')
xlabel('t''')
ylim([0 1])
legend('control', 'Dl2 strong', 'Dl2 weak','Location', 'northwest')


shadedErrorBar(norm_time,ave_nuc_T_c,std_nuc_T_c,{'color',rgb('black'),'LineWidth',2},1); hold on
shadedErrorBar(norm_time,ave_act_T_zg,std_act_T_zg,{'color',rgb('green'),'LineWidth',2},1);
title('Active Nuclei vs Time')
ylabel('Percent Active Nuclei')
xlabel('t''')
ylim([0 1])
legend('control', 'Zld strong','Location', 'northwest')


plot(norm_time,ave_nuc_T_c,'color',rgb('blue'),'LineWidth',2); hold on
plot(norm_time,ave_nuc_T_zg,'color',rgb('green'),'LineWidth',2)
title('Active Nuclei vs Time')
ylabel('Number Active Nuclei')
xlabel('t''')
ylim([0 610])
legend('control','Zld good','Location', 'northwest')


%het
shadedErrorBar(norm_time,ave_act_T_hetc,std_act_T_hetc,{'color',rgb('black'),'LineWidth',2},1); hold on
shadedErrorBar(norm_time,ave_act_T_hetb,std_act_T_hetb,{'color',rgb('silver'),'LineWidth',2},1); hold on 
title('Number of Active Nuclei')
ylabel('Number Active Nuclei')
xlabel('t''')
ylim([0 1])
xlim([.1 1])
legend('het control', 'het Dl1 strong', 'Location', 'northwest')

shadedErrorBar(norm_time,ave_nuc_T_hetc,std_nuc_T_hetc,{'color',rgb('black'),'LineWidth',2},1); hold on
shadedErrorBar(norm_time,ave_nuc_T_hetb,std_nuc_T_hetb,{'color',rgb('silver'),'LineWidth',2},1); hold on 
title('Number of Active Nuclei')
ylabel('Number Active Nuclei')
xlabel('t''')
ylim([0 610])
xlim([.1 1])
legend('het control', 'het Dl1 strong', 'Location', 'northwest')


shadedErrorBar(norm_time,ave_nuc_T_hetc,std_nuc_T_hetc,{'color',rgb('black'),'LineWidth',2},1); hold on
shadedErrorBar(norm_time,ave_nuc_T_hetdl2b,std_nuc_T_hetdl2b,{'color',rgb('purple'),'LineWidth',2},1); hold on 
shadedErrorBar(norm_time,ave_nuc_T_hetdl2d,std_nuc_T_hetdl2d,{'--','color',rgb('plum'),'LineWidth',2},1);
title('Number of Active Nuclei')
ylabel('Number Active Nuclei')
xlabel('t''')
ylim([0 610])
xlim([.1 1])
legend('het control', 'het Dl2 strong', 'het Dl2 weak','Location', 'northwest')

shadedErrorBar(norm_time,ave_nuc_T_hetc,std_nuc_T_hetc,{'color',rgb('black'),'LineWidth',2},1); hold on
shadedErrorBar(norm_time,ave_nuc_T_hetzg,std_nuc_T_hetzg,{'color',rgb('teal'),'LineWidth',2},1); hold on 
title('Number of Active Nuclei')
ylabel('Number Active Nuclei')
xlabel('t''')
ylim([0 610])
xlim([.1 1])
legend('het control', 'het Zld strong', 'Location', 'northwest')
%}

%Number active nuclei over time
figure;
shadedErrorBar(norm_time,ave_nuc_T_c,std_nuc_T_c,{'color',rgb('black'),'LineWidth',2}); hold on
shadedErrorBar(norm_time,ave_nuc_T_b,std_nuc_T_b,{'color',rgb('red'),'LineWidth',2});hold on
shadedErrorBar(norm_time,ave_nuc_T_d,std_nuc_T_d,{'--','color',rgb('lightcoral'),'LineWidth',2});
shadedErrorBar(norm_time,ave_nuc_T_dl2b,std_nuc_T_dl2b,{'color',rgb('blue'),'LineWidth',2});
shadedErrorBar(norm_time,ave_nuc_T_dl2d,std_nuc_T_dl2d,{'--','color',rgb('cornflowerblue'),'LineWidth',2});
shadedErrorBar(norm_time,ave_nuc_T_zg,std_nuc_T_zg,{'color',rgb('green'),'LineWidth',2});
shadedErrorBar(norm_time,ave_nuc_T_zw,std_nuc_T_zw,{'color',rgb('brown')});
title('Active Nuclei vs Time')
ylabel('Percent Active Nuclei')
xlabel('t''')
ylim([0 610])
legend('Control', 'Dl1 strong','Dl1 weak', 'Dl2 strong', 'Dl2 weak','zld good','zld weak','Location', 'northwest')

figure;
plot(norm_time,ave_nuc_T_c,'color',rgb('blue')); hold on
plot(norm_time,ave_nuc_T_b,'color',rgb('black'))
plot(norm_time,ave_nuc_T_d,'color',rgb('red'))
plot(norm_time,ave_nuc_T_dl2b,'color',rgb('darkorange'))
plot(norm_time,ave_nuc_T_dl2d,'color',rgb('purple'))
plot(norm_time,ave_nuc_T_zg,'color',rgb('green'))
plot(norm_time,ave_nuc_T_zw,'color',rgb('brown'))
title('Number Active Nuclei vs Time')
ylabel('Number of Active Nuclei')
xlabel('t''')
legend('control', 'dl1 strong','dl1 weak', 'dl2 strong', 'dl2 weak','zld good','zld weak','Location', 'northwest')

%half dorsal
figure;
errorbar(norm_time,ave_nuc_T_mcphet,std_nuc_T_mcphet,'color',rgb('blue')); hold on
errorbar(norm_time,ave_nuc_T_dlhetcontrol,std_nuc_T_dlhetcontrol,'color',rgb('black'));
errorbar(norm_time,ave_nuc_T_dlhetdl2best,std_nuc_T_dlhetdl2best,'color',rgb('red'));
plot(norm_time,ave_nuc_T_mcphet,'color',rgb('blue')); hold on
plot(norm_time,ave_nuc_T_dlhetcontrol,'color',rgb('black'))
plot(norm_time,ave_nuc_T_dlhetdl2best,'color',rgb('red'))
title('Number Active Nuclei vs Time')
ylabel('Number of Active Nuclei')
xlabel('t''')
ylim([0 550])
legend('mcp control', 'dl het control','dl het dl2 best','Location', 'northwest')

%Normalized activation over time
figure;
shadedErrorBar(norm_time,ave_act_T_c,std_act_T_c,{'color',rgb('black'),'LineWidth',2}); hold on
shadedErrorBar(norm_time,ave_act_T_b,std_act_T_b,{'color',rgb('red'),'LineWidth',2}); hold on
shadedErrorBar(norm_time,ave_act_T_d,std_act_T_d,{'--','color',rgb('lightcoral'),'LineWidth',2});
shadedErrorBar(norm_time,ave_act_T_dl2b,std_act_T_dl2b,{'color',rgb('blue'),'LineWidth',2});
shadedErrorBar(norm_time,ave_act_T_dl2d,std_act_T_dl2d,{'--','color',rgb('cornflowerblue'),'LineWidth',2});
shadedErrorBar(norm_time,ave_act_T_zg,std_act_T_zg,{'color',rgb('green'),'LineWidth',2});
title('Active Nuclei vs Time')
ylabel('Percent Active Nuclei')
xlabel('t''')
ylim([0 1])
xlim([.1 1])
legend('control', 'dl1 strong','dl1 weak', 'dl2 strong', 'dl2 weak','zld good','zld weak','Location', 'northwest')


%half dorsal
figure;
errorbar(norm_time,ave_act_T_mcphet,std_act_T_mcphet,'color',rgb('blue')); hold on
errorbar(norm_time,ave_act_T_dlhetcontrol,std_act_T_dlhetcontrol,'color',rgb('black'));
errorbar(norm_time,ave_act_T_dlhetdl2best,std_act_T_dlhetdl2best,'color',rgb('red'));
plot(norm_time,ave_act_T_mcphet,'color',rgb('blue')); hold on
plot(norm_time,ave_act_T_dlhetcontrol,'color',rgb('black'))
plot(norm_time,ave_act_T_dlhetdl2best,'color',rgb('red'))
title('Active Nuclei vs Time')
ylabel('Percent Active Nuclei')
xlabel('t''')
ylim([0 1])
legend('mcp control', 'dl het control','dl het dl2 best','Location', 'northwest')

figure;
plot(norm_time,ave_act_T_c,'color',rgb('black')); hold on
plot(norm_time,ave_act_T_b,'color',rgb('red'))
plot(norm_time,ave_act_T_d,'color',rgb('lightcoral'))
plot(norm_time,ave_act_T_dl2b,'color',rgb('blue'))
plot(norm_time,ave_act_T_dl2d,'color',rgb('cornflowerblue'))
plot(norm_time,ave_act_T_zg,'color',rgb('green'))
plot(norm_time,ave_act_T_zw,'color',rgb('brown'))
title('Active Nuclei vs Time')
ylabel('Percent Active Nuclei')
xlabel('t''')
ylim([0 1])
legend('control', 'dl1 strong','dl1 weak', 'dl2 strong', 'dl2 weak','zld good','Location', 'northwest')


%% Instantaneous Activation
for i = 1:length(T_controln1)
    inst_c{1}(i) = length(find(M_controln1(i,activeM_controln1) > .8*thresh_base_controln1));
end
for i = 1:length(T_controln2)
    inst_c{2}(i) = length(find(M_controln2(i,activeM_controln2) > .8*thresh_base_controln2));
end
for i = 1:length(T_controln4)
    inst_c{3}(i) = length(find(M_controln4(i,activeM_controln4) > .8*thresh_base_controln4));
end
for i = 1:length(T_bestn1)
    inst_b{1}(i) = length(find(M_bestn1(i,activeM_bestn1) > .8*thresh_base_dl1bestn1));
end
for i = 1:length(T_bestn2)
    inst_b{2}(i) = length(find(M_bestn2(i,activeM_bestn2) > .8*thresh_base_dl1bestn2));
end
for i = 1:length(T_bestn4)
    inst_b{3}(i) = length(find(M_bestn4(i,activeM_bestn4) > .8*thresh_base_dl1bestn4));
end
for i = 1:length(T_deln1)
    inst_d{1}(i) = length(find(M_deln1(i,activeM_deln1) > .8*thresh_base_dl1deln1));
end
for i = 1:length(T_deln2)
    inst_d{2}(i) = length(find(M_deln2(i,activeM_deln2) > .8*thresh_base_dl1deln2));
end
for i = 1:length(T_deln4)
    inst_d{3}(i) = length(find(M_deln4(i,activeM_deln4) > .8*thresh_base_dl1deln4));
end
for i = 1:length(T_dl2bestn1)
    inst_dl2b{1}(i) = length(find(M_dl2bestn1(i,activeM_dl2bestn1) > .8*thresh_base_dl2bestn1));
end
for i = 1:length(T_dl2bestn2)
    inst_dl2b{2}(i) = length(find(M_dl2bestn2(i,activeM_dl2bestn2) > .8*thresh_base_dl2bestn2));
end
for i = 1:length(T_dl2bestn5)
    inst_dl2b{3}(i) = length(find(M_dl2bestn5(i,activeM_dl2bestn5) > .8*thresh_base_dl2bestn5));
end
for i = 1:length(T_dl2deln1)
    inst_dl2d{1}(i) = length(find(M_dl2deln1(i,activeM_dl2deln1) > .8*thresh_base_dl2deln1));
end
for i = 1:length(T_dl2deln2)
    inst_dl2d{2}(i) = length(find(M_dl2deln2(i,activeM_dl2deln2) > .8*thresh_base_dl2deln2));
end
for i = 1:length(T_dl2deln4)
    inst_dl2d{3}(i) = length(find(M_dl2deln4(i,activeM_dl2deln4) > .8*thresh_base_dl2deln4));
end
for i = 1:length(T_zldgoodn2)
    inst_zg{1}(i) = length(find(M_zldgoodn2(i,activeM_zldgoodn2) > .8*thresh_base_zldgoodn2));
end
for i = 1:length(T_zldgoodn3)
    inst_zg{2}(i) = length(find(M_zldgoodn3(i,activeM_zldgoodn3) > .8*thresh_base_zldgoodn3));
end
for i = 1:length(T_zldgoodn4)
    inst_zg{3}(i) = length(find(M_zldgoodn4(i,activeM_zldgoodn4) > .8*thresh_base_zldgoodn4));
end

inst_norm_c(1,:) = interp1(T_controln1/T_controln1(end),inst_c{1},norm_time);
inst_norm_c(2,:) = interp1(T_controln2/T_controln2(end),inst_c{2},norm_time);
inst_norm_c(3,:) = interp1(T_controln4/T_controln4(end),inst_c{3},norm_time);
inst_norm_b(1,:) = interp1(T_bestn1/T_bestn1(end),inst_b{1},norm_time);
inst_norm_b(2,:) = interp1(T_bestn2/T_bestn2(end),inst_b{2},norm_time);
inst_norm_b(3,:) = interp1(T_bestn4/T_bestn4(end),inst_b{3},norm_time);
inst_norm_d(1,:) = interp1(T_deln1/T_deln1(end),inst_d{1},norm_time);
inst_norm_d(2,:) = interp1(T_deln2/T_deln2(end),inst_d{2},norm_time);
inst_norm_d(3,:) = interp1(T_deln4/T_deln4(end),inst_d{3},norm_time);
inst_norm_dl2b(1,:) = interp1(T_dl2bestn1/T_dl2bestn1(end),inst_dl2b{1},norm_time);
inst_norm_dl2b(2,:) = interp1(T_dl2bestn2/T_dl2bestn2(end),inst_dl2b{2},norm_time);
inst_norm_dl2b(3,:) = interp1(T_dl2bestn5/T_dl2bestn5(end),inst_dl2b{3},norm_time);
inst_norm_dl2d(1,:) = interp1(T_dl2deln1/T_dl2deln1(end),inst_dl2d{1},norm_time);
inst_norm_dl2d(2,:) = interp1(T_dl2deln2/T_dl2deln2(end),inst_dl2d{2},norm_time);
inst_norm_dl2d(3,:) = interp1(T_dl2deln4/T_dl2deln4(end),inst_dl2d{3},norm_time);
inst_norm_zg(1,:) = interp1(T_zldgoodn2/T_zldgoodn2(end),inst_zg{1},norm_time);
inst_norm_zg(2,:) = interp1(T_zldgoodn3/T_zldgoodn3(end),inst_zg{2},norm_time);
inst_norm_zg(3,:) = interp1(T_zldgoodn4/T_zldgoodn4(end),inst_zg{3},norm_time);

ave_inst_c = nanmean(inst_norm_c,1);
std_inst_c = nanstd(inst_norm_c);
ave_inst_b = nanmean(inst_norm_b,1);
std_inst_b = nanstd(inst_norm_b);
ave_inst_d = nanmean(inst_norm_d,1);
std_inst_d = nanstd(inst_norm_d);
ave_inst_dl2b = nanmean(inst_norm_dl2b,1);
std_inst_dl2b = nanstd(inst_norm_dl2b);
ave_inst_dl2d = nanmean(inst_norm_dl2d,1);
std_inst_dl2d = nanstd(inst_norm_dl2d);
ave_inst_zg = nanmean(inst_norm_zg,1);
std_inst_zg = nanstd(inst_norm_zg);


%het
for i = 1:length(T_hetcontroln1)
    inst_hetc{1}(i) = length(find(M_hetcontroln1(i,activeM_hetcontroln1) > .8*thresh_base_hetcontroln1));
end
for i = 1:length(T_hetcontroln2)
    inst_hetc{2}(i) = length(find(M_hetcontroln2(i,activeM_hetcontroln2) > .8*thresh_base_hetcontroln2));
end
for i = 1:length(T_hetcontroln3)
    inst_hetc{3}(i) = length(find(M_hetcontroln3(i,activeM_hetcontroln3) > .8*thresh_base_hetcontroln3));
end
for i = 1:length(T_hetbestn1)
    inst_hetb{1}(i) = length(find(M_hetbestn1(i,activeM_hetbestn1) > .8*thresh_base_hetbestn1));
end
for i = 1:length(T_hetbestn2)
    inst_hetb{2}(i) = length(find(M_hetbestn2(i,activeM_hetbestn2) > .8*thresh_base_hetbestn2));
end
for i = 1:length(T_hetbestn3)
    inst_hetb{3}(i) = length(find(M_hetbestn3(i,activeM_hetbestn3) > .8*thresh_base_hetbestn3));
end
for i = 1:length(T_hetdl2bestn1)
    inst_hetdl2b{1}(i) = length(find(M_hetdl2bestn1(i,activeM_hetdl2bestn1) > .8*thresh_base_hetdl2bestn1));
end
for i = 1:length(T_hetdl2bestn2)
    inst_hetdl2b{2}(i) = length(find(M_hetdl2bestn2(i,activeM_hetdl2bestn2) > .8*thresh_base_hetdl2bestn2));
end
for i = 1:length(T_hetdl2bestn3)
    inst_hetdl2b{3}(i) = length(find(M_hetdl2bestn3(i,activeM_hetdl2bestn3) > .8*thresh_base_hetdl2bestn3));
end
for i = 1:length(T_hetdl2deln1)
    inst_hetdl2d{1}(i) = length(find(M_hetdl2deln1(i,activeM_hetdl2deln1) > .8*thresh_base_hetdl2deln1));
end
for i = 1:length(T_hetdl2deln2)
    inst_hetdl2d{2}(i) = length(find(M_hetdl2deln2(i,activeM_hetdl2deln2) > .8*thresh_base_hetdl2deln2));
end
for i = 1:length(T_hetdl2deln3)
    inst_hetdl2d{3}(i) = length(find(M_hetdl2deln3(i,activeM_hetdl2deln3) > .8*thresh_base_hetdl2deln3));
end
for i = 1:length(T_hetzldbestn1)
    inst_hetzg{1}(i) = length(find(M_hetzldbestn1(i,activeM_hetzldbestn2) > .8*thresh_base_hetzldbestn2));
end
for i = 1:length(T_hetzldbestn2)
    inst_hetzg{2}(i) = length(find(M_hetzldbestn3(i,activeM_hetzldbestn3) > .8*thresh_base_hetzldbestn3));
end
for i = 1:length(T_hetzldbestn3)
    inst_hetzg{3}(i) = length(find(M_hetzldbestn3(i,activeM_hetzldbestn3) > .8*thresh_base_hetzldbestn3));
end

inst_norm_hetc(1,:) = interp1(T_hetcontroln1/T_hetcontroln1(end),inst_hetc{1},norm_time);
inst_norm_hetc(2,:) = interp1(T_hetcontroln2/T_hetcontroln2(end),inst_hetc{2},norm_time);
inst_norm_hetc(3,:) = interp1(T_hetcontroln3/T_hetcontroln3(end),inst_hetc{3},norm_time);
inst_norm_hetb(1,:) = interp1(T_hetbestn1/T_hetbestn1(end),inst_hetb{1},norm_time);
inst_norm_hetb(2,:) = interp1(T_hetbestn2/T_hetbestn2(end),inst_hetb{2},norm_time);
inst_norm_hetb(3,:) = interp1(T_hetbestn3/T_hetbestn3(end),inst_hetb{3},norm_time);
inst_norm_hetdl2b(1,:) = interp1(T_hetdl2bestn1/T_hetdl2bestn1(end),inst_hetdl2b{1},norm_time);
inst_norm_hetdl2b(2,:) = interp1(T_hetdl2bestn2/T_hetdl2bestn2(end),inst_hetdl2b{2},norm_time);
inst_norm_hetdl2b(3,:) = interp1(T_hetdl2bestn3/T_hetdl2bestn3(end),inst_hetdl2b{3},norm_time);
inst_norm_hetdl2d(1,:) = interp1(T_hetdl2deln1/T_hetdl2deln1(end),inst_hetdl2d{1},norm_time);
inst_norm_hetdl2d(2,:) = interp1(T_hetdl2deln2/T_hetdl2deln2(end),inst_hetdl2d{2},norm_time);
inst_norm_hetdl2d(3,:) = interp1(T_hetdl2deln3/T_hetdl2deln3(end),inst_hetdl2d{3},norm_time);
inst_norm_hetzg(1,:) = interp1(T_hetzldbestn1/T_hetzldbestn1(end),inst_hetzg{1},norm_time);
inst_norm_hetzg(2,:) = interp1(T_hetzldbestn2/T_hetzldbestn2(end),inst_hetzg{2},norm_time);
inst_norm_hetzg(3,:) = interp1(T_hetzldbestn3/T_hetzldbestn3(end),inst_hetzg{3},norm_time);

ave_inst_hetc = nanmean(inst_norm_hetc,1);
std_inst_hetc = nanstd(inst_norm_hetc);
ave_inst_hetb = nanmean(inst_norm_hetb,1);
std_inst_hetb = nanstd(inst_norm_hetb);
ave_inst_hetdl2b = nanmean(inst_norm_hetdl2b,1);
std_inst_hetdl2b = nanstd(inst_norm_hetdl2b);
ave_inst_hetdl2d = nanmean(inst_norm_hetdl2d,1);
std_inst_hetdl2d = nanstd(inst_norm_hetdl2d);
ave_inst_hetzg = nanmean(inst_norm_hetzg,1);
std_inst_hetzg = nanstd(inst_norm_hetzg);

%Revision
for i = 1:length(T_Rcn1)
    inst_Rc{1}(i) = length(find(M_Rcn1(i,activeM_Rcn1) > .8*thresh_base_Rcn1));
end
for i = 1:length(T_Rcn2)
    inst_Rc{2}(i) = length(find(M_Rcn2(i,activeM_Rcn2) > .8*thresh_base_Rcn2));
end
for i = 1:length(T_Rcn3)
    inst_Rc{3}(i) = length(find(M_Rcn3(i,activeM_Rcn3) > .8*thresh_base_Rcn3));
end

for i = 1:length(T_Rdzn1)
    inst_Rdz{1}(i) = length(find(M_Rdzn1(i,activeM_Rdzn1) > .8*thresh_base_Rdzn1));
end
for i = 1:length(T_Rdzn2)
    inst_Rdz{2}(i) = length(find(M_Rdzn2(i,activeM_Rdzn2) > .8*thresh_base_Rdzn2));
end
for i = 1:length(T_Rdzn3)
    inst_Rdz{3}(i) = length(find(M_Rdzn3(i,activeM_Rdzn3) > .8*thresh_base_Rdzn3));
end

for i = 1:length(T_Rdn1)
    inst_Rd{1}(i) = length(find(M_Rdn1(i,activeM_Rdn1) > .8*thresh_base_Rdn1));
end
for i = 1:length(T_Rdn2)
    inst_Rd{2}(i) = length(find(M_Rdn2(i,activeM_Rdn2) > .8*thresh_base_Rdn2));
end
for i = 1:length(T_Rdn3)
    inst_Rd{3}(i) = length(find(M_Rdn3(i,activeM_Rdn3) > .8*thresh_base_Rdn3));
end


inst_norm_Rc(1,:) = interp1(T_Rcn1/T_Rcn1(end),inst_Rc{1},norm_time);
inst_norm_Rc(2,:) = interp1(T_Rcn2/T_Rcn2(end),inst_Rc{2},norm_time);
inst_norm_Rc(3,:) = interp1(T_Rcn3/T_Rcn3(end),inst_Rc{3},norm_time);
inst_norm_Rdz(1,:) = interp1(T_Rdzn1/T_Rdzn1(end),inst_Rdz{1},norm_time);
inst_norm_Rdz(2,:) = interp1(T_Rdzn2/T_Rdzn2(end),inst_Rdz{2},norm_time);
inst_norm_Rdz(3,:) = interp1(T_Rdzn3/T_Rdzn3(end),inst_Rdz{3},norm_time);
inst_norm_Rd(1,:) = interp1(T_Rdn1/T_Rdn1(end),inst_Rd{1},norm_time);
inst_norm_Rd(2,:) = interp1(T_Rdn2/T_Rdn2(end),inst_Rd{2},norm_time);
inst_norm_Rd(3,:) = interp1(T_Rdn3/T_Rdn3(end),inst_Rd{3},norm_time);

ave_inst_Rc = nanmean(inst_norm_Rc,1);
std_inst_Rc = nanstd(inst_norm_Rc);
ave_inst_Rdz = nanmean(inst_norm_Rdz,1);
std_inst_Rdz = nanstd(inst_norm_Rdz);
ave_inst_Rd = nanmean(inst_norm_Rd,1);
std_inst_Rd = nanstd(inst_norm_Rd);


figure;
shadedErrorBar(norm_time,ave_inst_c,std_inst_c,{'color',rgb('black'),'LineWidth',2},1); hold on
shadedErrorBar(norm_time,ave_inst_b,std_inst_b,{'color',rgb('red'),'LineWidth',2},1); hold on
shadedErrorBar(norm_time,ave_inst_d,std_inst_d,{'--','color',rgb('lightcoral'),'LineWidth',2},1); hold on
title('Active Nuclei vs Time')
ylabel('Number Active Nuclei')
xlabel('t''')
ylim([0 450])
xlim([.1 1])
legend('control', 'Dl1 strong','Dl1 weak','Location', 'northwest')


figure;
shadedErrorBar(norm_time,ave_inst_c,std_inst_c,{'color',rgb('black'),'LineWidth',2},1); hold on
shadedErrorBar(norm_time,ave_inst_dl2b,std_inst_dl2b,{'color',rgb('blue'),'LineWidth',2},1); hold on
shadedErrorBar(norm_time,ave_inst_dl2d,std_inst_dl2d,{'--','color',rgb('cornflowerblue'),'LineWidth',2},1); hold on
title('Active Nuclei vs Time')
ylabel('Number Active Nuclei')
xlabel('t''')
ylim([0 450])
xlim([.1 1])
legend('control', 'Dl2 strong','Dl2 weak','Location', 'northwest')


figure;
shadedErrorBar(norm_time,ave_inst_c,std_inst_c,{'color',rgb('black'),'LineWidth',2},1); hold on
shadedErrorBar(norm_time,ave_inst_zg,std_inst_zg,{'color',rgb('green'),'LineWidth',2},1); hold on
title('Active Nuclei vs Time')
ylabel('Number Active Nuclei')
xlabel('t''')
ylim([0 450])
xlim([.1 1])
legend('control','Zld strong','Location', 'northwest')

%het
figure;
shadedErrorBar(norm_time,ave_inst_hetc,std_inst_hetc,{'color',rgb('black'),'LineWidth',2},1); hold on
shadedErrorBar(norm_time,ave_inst_hetb,std_inst_hetb,{'color',rgb('silver'),'LineWidth',2},1); hold on
title('Active Nuclei vs Time')
ylabel('Number Active Nuclei')
xlabel('t''')
ylim([0 450])
xlim([.1 1])
legend('het control', 'het Dl1 strong','Location', 'northwest')


figure;
shadedErrorBar(norm_time,ave_inst_hetc,std_inst_hetc,{'color',rgb('black'),'LineWidth',2},1); hold on
shadedErrorBar(norm_time,ave_inst_hetdl2b,std_inst_hetdl2b,{'color',rgb('purple'),'LineWidth',2},1); hold on
shadedErrorBar(norm_time,ave_inst_hetdl2d,std_inst_hetdl2d,{'--','color',rgb('plum'),'LineWidth',2},1); hold on
title('Active Nuclei vs Time')
ylabel('Number Active Nuclei')
xlabel('t''')
ylim([0 450])
xlim([.1 1])
legend('het control', 'het Dl2 strong','het Dl2 weak','Location', 'northwest')


figure;
shadedErrorBar(norm_time,ave_inst_hetc,std_inst_hetc,{'color',rgb('black'),'LineWidth',2},1); hold on
shadedErrorBar(norm_time,ave_inst_hetzg,std_inst_hetzg,{'color',rgb('teal'),'LineWidth',2},1); hold on
title('Active Nuclei vs Time')
ylabel('Number Active Nuclei')
xlabel('t''')
ylim([0 450])
xlim([.1 1])
legend('het control','het Zld strong','Location', 'northwest')

%% Nomralizing Signal and Output time

for i = 1:bins_bestn1
    ave_Signal_controln1(:,i) = interp1(time_controln1,aveSignal_bin_controln1(:,i),norm_time);
    ave_Signal_controln2(:,i) = interp1(time_controln2,aveSignal_bin_controln2(:,i),norm_time);
    ave_Signal_controln4(:,i) = interp1(time_controln4,aveSignal_bin_controln4(:,i),norm_time);
    ave_Signal_controln3(:,i) = interp1(time_controln3,aveSignal_bin_controln3(:,i),norm_time);
    ave_Signal_controln5(:,i) = interp1(time_controln5,aveSignal_bin_controln5(:,i),norm_time);
    ave_Signal_controln6(:,i) = interp1(time_controln6,aveSignal_bin_controln6(:,i),norm_time);
    
    ave_Signal_bestn1(:,i) = interp1(time_bestn1,aveSignal_bin_bestn1(:,i),norm_time);
    ave_Signal_bestn2(:,i) = interp1(time_bestn2,aveSignal_bin_bestn2(:,i),norm_time);
    ave_Signal_bestn4(:,i) = interp1(time_bestn4,aveSignal_bin_bestn4(:,i),norm_time);
    
    ave_Signal_deln1(:,i) = interp1(time_deln1,aveSignal_bin_deln1(:,i),norm_time);
    ave_Signal_deln2(:,i) = interp1(time_deln2,aveSignal_bin_deln2(:,i),norm_time);
    ave_Signal_deln4(:,i) = interp1(time_deln4,aveSignal_bin_deln4(:,i),norm_time);
    
    ave_Signal_dl2bestn1(:,i) = interp1(time_dl2bestn1,aveSignal_bin_dl2bestn1(:,i),norm_time);
    ave_Signal_dl2bestn2(:,i) = interp1(time_dl2bestn2,aveSignal_bin_dl2bestn2(:,i),norm_time);
    ave_Signal_dl2bestn5(:,i) = interp1(time_dl2bestn5,aveSignal_bin_dl2bestn5(:,i),norm_time);
    
    ave_Signal_dl2deln1(:,i) = interp1(time_dl2deln1,aveSignal_bin_dl2deln1(:,i),norm_time);
    ave_Signal_dl2deln2(:,i) = interp1(time_dl2deln2,aveSignal_bin_dl2deln2(:,i),norm_time);
    ave_Signal_dl2deln4(:,i) = interp1(time_dl2deln4,aveSignal_bin_dl2deln4(:,i),norm_time);
    
    ave_Signal_zldgoodn2(:,i) = interp1(time_zldgoodn2,aveSignal_bin_zldgoodn2(:,i),norm_time);
    ave_Signal_zldgoodn3(:,i) = interp1(time_zldgoodn3,aveSignal_bin_zldgoodn3(:,i),norm_time);
    ave_Signal_zldgoodn4(:,i) = interp1(time_zldgoodn4,aveSignal_bin_zldgoodn4(:,i),norm_time);
    
    ave_Signal_zldweakn1(:,i) = interp1(time_zldweakn1,aveSignal_bin_zldweakn1(:,i),norm_time);
    ave_Signal_zldweakn2(:,i) = interp1(time_zldweakn2,aveSignal_bin_zldweakn2(:,i),norm_time);
    ave_Signal_zldweakn5(:,i) = interp1(time_zldweakn5,aveSignal_bin_zldweakn5(:,i),norm_time);
    
    ave_Signal_mcphetn1(:,i) = interp1(time_mcphetn1,aveSignal_bin_mcphetn1(:,i),norm_time);
    ave_Signal_mcphetn2(:,i) = interp1(time_mcphetn2,aveSignal_bin_mcphetn2(:,i),norm_time);
    ave_Signal_mcphetn3(:,i) = interp1(time_mcphetn3,aveSignal_bin_mcphetn3(:,i),norm_time);
    
    ave_Signal_dlhetcontroln1(:,i) = interp1(time_dlhetcontroln1,aveSignal_bin_dlhetcontroln1(:,i),norm_time);
    ave_Signal_dlhetcontroln3(:,i) = interp1(time_dlhetcontroln3,aveSignal_bin_dlhetcontroln3(:,i),norm_time);
    ave_Signal_dlhetcontroln4(:,i) = interp1(time_dlhetcontroln4,aveSignal_bin_dlhetcontroln4(:,i),norm_time);

    ave_Signal_dlhetdl2bestn2(:,i) = interp1(time_dlhetdl2bestn2,aveSignal_bin_dlhetdl2bestn2(:,i),norm_time);
    ave_Signal_dlhetdl2bestn3(:,i) = interp1(time_dlhetdl2bestn3,aveSignal_bin_dlhetdl2bestn3(:,i),norm_time);
    ave_Signal_dlhetdl2bestn5(:,i) = interp1(time_dlhetdl2bestn5,aveSignal_bin_dlhetdl2bestn5(:,i),norm_time);

    ave_Signal_twideln1(:,i) = interp1(time_twideln1,aveSignal_bin_twideln1(:,i),norm_time);
    ave_Signal_twideln2(:,i) = interp1(time_twideln2,aveSignal_bin_twideln2(:,i),norm_time);

    %het signal
    ave_Signal_hetcontroln1(:,i) = interp1(time_hetcontroln1,aveSignal_bin_hetcontroln1(:,i),norm_time);
    ave_Signal_hetcontroln2(:,i) = interp1(time_hetcontroln2,aveSignal_bin_hetcontroln2(:,i),norm_time);
    ave_Signal_hetcontroln3(:,i) = interp1(time_hetcontroln3,aveSignal_bin_hetcontroln3(:,i),norm_time);
    
    ave_Signal_hetbestn1(:,i) = interp1(time_hetbestn1,aveSignal_bin_hetbestn1(:,i),norm_time);
    ave_Signal_hetbestn2(:,i) = interp1(time_hetbestn2,aveSignal_bin_hetbestn2(:,i),norm_time);
    ave_Signal_hetbestn3(:,i) = interp1(time_hetbestn3,aveSignal_bin_hetbestn3(:,i),norm_time);
        
    ave_Signal_hetdl2bestn1(:,i) = interp1(time_hetdl2bestn1,aveSignal_bin_hetdl2bestn1(:,i),norm_time);
    ave_Signal_hetdl2bestn2(:,i) = interp1(time_hetdl2bestn2,aveSignal_bin_hetdl2bestn2(:,i),norm_time);
    ave_Signal_hetdl2bestn3(:,i) = interp1(time_hetdl2bestn3,aveSignal_bin_hetdl2bestn3(:,i),norm_time);
    
    ave_Signal_hetdl2deln1(:,i) = interp1(time_hetdl2deln1,aveSignal_bin_hetdl2deln1(:,i),norm_time);
    ave_Signal_hetdl2deln2(:,i) = interp1(time_hetdl2deln2,aveSignal_bin_hetdl2deln2(:,i),norm_time);
    ave_Signal_hetdl2deln3(:,i) = interp1(time_hetdl2deln3,aveSignal_bin_hetdl2deln3(:,i),norm_time);
    
    ave_Signal_hetzldbestn1(:,i) = interp1(time_hetzldbestn1,aveSignal_bin_hetzldbestn1(:,i),norm_time);
    ave_Signal_hetzldbestn2(:,i) = interp1(time_hetzldbestn2,aveSignal_bin_hetzldbestn2(:,i),norm_time);
    ave_Signal_hetzldbestn3(:,i) = interp1(time_hetzldbestn3,aveSignal_bin_hetzldbestn3(:,i),norm_time);
    
    %Revision
    ave_Signal_Rcontroln1(:,i) = interp1(time_Rcontroln1,aveSignal_bin_Rcn1(:,i),norm_time);
    ave_Signal_Rcontroln2(:,i) = interp1(time_Rcontroln2,aveSignal_bin_Rcn2(:,i),norm_time);
    ave_Signal_Rcontroln3(:,i) = interp1(time_Rcontroln3,aveSignal_bin_Rcn3(:,i),norm_time);
    
    ave_Signal_Rdzn1(:,i) = interp1(time_Rdzn1,aveSignal_bin_Rdzn1(:,i),norm_time);
    ave_Signal_Rdzn2(:,i) = interp1(time_Rdzn2,aveSignal_bin_Rdzn2(:,i),norm_time);
    ave_Signal_Rdzn3(:,i) = interp1(time_Rdzn3,aveSignal_bin_Rdzn3(:,i),norm_time);
    
    ave_Signal_Rdn1(:,i) = interp1(time_Rdn1,aveSignal_bin_Rdn1(:,i),norm_time);
    ave_Signal_Rdn2(:,i) = interp1(time_Rdn2,aveSignal_bin_Rdn2(:,i),norm_time);
    ave_Signal_Rdn3(:,i) = interp1(time_Rdn3,aveSignal_bin_Rdn3(:,i),norm_time);
    
    %Output
    ave_Area_controln1(:,i) = interp1(time_controln1(2:end),aveArea_bin_controln1(:,i),norm_time);
    ave_Area_controln2(:,i) = interp1(time_controln2(2:end),aveArea_bin_controln2(:,i),norm_time);
    ave_Area_controln4(:,i) = interp1(time_controln4(2:end),aveArea_bin_controln4(:,i),norm_time);
    ave_Area_controln3(:,i) = interp1(time_controln3(2:end),aveArea_bin_controln3(:,i),norm_time);
    ave_Area_controln5(:,i) = interp1(time_controln5(2:end),aveArea_bin_controln5(:,i),norm_time);
    ave_Area_controln6(:,i) = interp1(time_controln6(2:end),aveArea_bin_controln6(:,i),norm_time);
    
    ave_Area_bestn1(:,i) = interp1(time_bestn1(2:end),aveArea_bin_bestn1(:,i),norm_time);
    ave_Area_bestn2(:,i) = interp1(time_bestn2(2:end),aveArea_bin_bestn2(:,i),norm_time);
    ave_Area_bestn4(:,i) = interp1(time_bestn4(2:end),aveArea_bin_bestn4(:,i),norm_time);
    
    ave_Area_deln1(:,i) = interp1(time_deln1(2:end),aveArea_bin_deln1(:,i),norm_time);
    ave_Area_deln2(:,i) = interp1(time_deln2(2:end),aveArea_bin_deln2(:,i),norm_time);
    ave_Area_deln4(:,i) = interp1(time_deln4(2:end),aveArea_bin_deln4(:,i),norm_time);
    
    ave_Area_dl2bestn1(:,i) = interp1(time_dl2bestn1(2:end),aveArea_bin_dl2bestn1(:,i),norm_time);
    ave_Area_dl2bestn2(:,i) = interp1(time_dl2bestn2(2:end),aveArea_bin_dl2bestn2(:,i),norm_time);
    ave_Area_dl2bestn5(:,i) = interp1(time_dl2bestn5(2:end),aveArea_bin_dl2bestn5(:,i),norm_time);
    
    ave_Area_dl2deln1(:,i) = interp1(time_dl2deln1(2:end),aveArea_bin_dl2deln1(:,i),norm_time);
    ave_Area_dl2deln2(:,i) = interp1(time_dl2deln2(2:end),aveArea_bin_dl2deln2(:,i),norm_time);
    ave_Area_dl2deln4(:,i) = interp1(time_dl2deln4(2:end),aveArea_bin_dl2deln4(:,i),norm_time);
    
    ave_Area_zldgoodn2(:,i) = interp1(time_zldgoodn2(2:end),aveArea_bin_zldgoodn2(:,i),norm_time);
    ave_Area_zldgoodn3(:,i) = interp1(time_zldgoodn3(2:end),aveArea_bin_zldgoodn3(:,i),norm_time);
    ave_Area_zldgoodn4(:,i) = interp1(time_zldgoodn4(2:end),aveArea_bin_zldgoodn4(:,i),norm_time);

    ave_Area_zldweakn1(:,i) = interp1(time_zldweakn1(2:end),aveArea_bin_zldweakn1(:,i),norm_time);
    ave_Area_zldweakn2(:,i) = interp1(time_zldweakn2(2:end),aveArea_bin_zldweakn2(:,i),norm_time);
    ave_Area_zldweakn5(:,i) = interp1(time_zldweakn5(2:end),aveArea_bin_zldweakn5(:,i),norm_time);
    
    ave_Area_mcphetn1(:,i) = interp1(time_mcphetn1(2:end),aveArea_bin_mcphetn1(:,i),norm_time);
    ave_Area_mcphetn2(:,i) = interp1(time_mcphetn2(2:end),aveArea_bin_mcphetn2(:,i),norm_time);
    ave_Area_mcphetn3(:,i) = interp1(time_mcphetn3(2:end),aveArea_bin_mcphetn3(:,i),norm_time);
    
    ave_Area_dlhetcontroln1(:,i) = interp1(time_dlhetcontroln1(2:end),aveArea_bin_dlhetcontroln1(:,i),norm_time);
    ave_Area_dlhetcontroln3(:,i) = interp1(time_dlhetcontroln3(2:end),aveArea_bin_dlhetcontroln3(:,i),norm_time);
    ave_Area_dlhetcontroln4(:,i) = interp1(time_dlhetcontroln4(2:end),aveArea_bin_dlhetcontroln4(:,i),norm_time);

    ave_Area_dlhetdl2bestn2(:,i) = interp1(time_dlhetdl2bestn2(2:end),aveArea_bin_dlhetdl2bestn2(:,i),norm_time);
    ave_Area_dlhetdl2bestn3(:,i) = interp1(time_dlhetdl2bestn3(2:end),aveArea_bin_dlhetdl2bestn3(:,i),norm_time);
    ave_Area_dlhetdl2bestn5(:,i) = interp1(time_dlhetdl2bestn5(2:end),aveArea_bin_dlhetdl2bestn5(:,i),norm_time);

    ave_Area_twideln1(:,i) = interp1(time_twideln1(2:end),aveArea_bin_twideln1(:,i),norm_time);
    ave_Area_twideln2(:,i) = interp1(time_twideln2(2:end),aveArea_bin_twideln2(:,i),norm_time);

    %het output
    ave_Area_hetcontroln1(:,i) = interp1(time_hetcontroln1(2:end),aveArea_bin_hetcontroln1(:,i),norm_time);
    ave_Area_hetcontroln2(:,i) = interp1(time_hetcontroln2(2:end),aveArea_bin_hetcontroln2(:,i),norm_time);
    ave_Area_hetcontroln3(:,i) = interp1(time_hetcontroln3(2:end),aveArea_bin_hetcontroln3(:,i),norm_time);
    
    ave_Area_hetbestn1(:,i) = interp1(time_hetbestn1(2:end),aveArea_bin_hetbestn1(:,i),norm_time);
    ave_Area_hetbestn2(:,i) = interp1(time_hetbestn2(2:end),aveArea_bin_hetbestn2(:,i),norm_time);
    ave_Area_hetbestn3(:,i) = interp1(time_hetbestn3(2:end),aveArea_bin_hetbestn3(:,i),norm_time);
        
    ave_Area_hetdl2bestn1(:,i) = interp1(time_hetdl2bestn1(2:end),aveArea_bin_hetdl2bestn1(:,i),norm_time);
    ave_Area_hetdl2bestn2(:,i) = interp1(time_hetdl2bestn2(2:end),aveArea_bin_hetdl2bestn2(:,i),norm_time);
    ave_Area_hetdl2bestn3(:,i) = interp1(time_hetdl2bestn3(2:end),aveArea_bin_hetdl2bestn3(:,i),norm_time);
    
    ave_Area_hetdl2deln1(:,i) = interp1(time_hetdl2deln1(2:end),aveArea_bin_hetdl2deln1(:,i),norm_time);
    ave_Area_hetdl2deln2(:,i) = interp1(time_hetdl2deln2(2:end),aveArea_bin_hetdl2deln2(:,i),norm_time);
    ave_Area_hetdl2deln3(:,i) = interp1(time_hetdl2deln3(2:end),aveArea_bin_hetdl2deln3(:,i),norm_time);
    
    ave_Area_hetzldbestn1(:,i) = interp1(time_hetzldbestn1(2:end),aveArea_bin_hetzldbestn1(:,i),norm_time);
    ave_Area_hetzldbestn2(:,i) = interp1(time_hetzldbestn2(2:end),aveArea_bin_hetzldbestn2(:,i),norm_time);
    ave_Area_hetzldbestn3(:,i) = interp1(time_hetzldbestn3(2:end),aveArea_bin_hetzldbestn3(:,i),norm_time);
    
    %Revision Output
    ave_Area_Rcontroln1(:,i) = interp1(time_Rcontroln1(2:end),aveArea_bin_Rcn1(:,i),norm_time);
    ave_Area_Rcontroln2(:,i) = interp1(time_Rcontroln2(2:end),aveArea_bin_Rcn2(:,i),norm_time);
    ave_Area_Rcontroln3(:,i) = interp1(time_Rcontroln3(2:end),aveArea_bin_Rcn3(:,i),norm_time);
    
    ave_Area_Rdzn1(:,i) = interp1(time_Rdzn1(2:end),aveArea_bin_Rdzn1(:,i),norm_time);
    ave_Area_Rdzn2(:,i) = interp1(time_Rdzn2(2:end),aveArea_bin_Rdzn2(:,i),norm_time);
    ave_Area_Rdzn3(:,i) = interp1(time_Rdzn3(2:end),aveArea_bin_Rdzn3(:,i),norm_time);
    
    ave_Area_Rdn1(:,i) = interp1(time_Rdn1(2:end),aveArea_bin_Rdn1(:,i),norm_time);
    ave_Area_Rdn2(:,i) = interp1(time_Rdn2(2:end),aveArea_bin_Rdn2(:,i),norm_time);
    ave_Area_Rdn3(:,i) = interp1(time_Rdn3(2:end),aveArea_bin_Rdn3(:,i),norm_time);
    
    %Nuclei Output
    nuc_output_controln1(i) = aveArea_bin_controln1(end,i)/length(bin_controln1{i});
    nuc_output_controln2(i) = aveArea_bin_controln2(end,i)/length(bin_controln2{i});
    nuc_output_controln4(i) = aveArea_bin_controln4(end,i)/length(bin_controln4{i});
    nuc_output_controln3(i) = aveArea_bin_controln3(end,i)/length(bin_controln3{i});
    nuc_output_controln5(i) = aveArea_bin_controln5(end,i)/length(bin_controln5{i});
    nuc_output_controln6(i) = aveArea_bin_controln6(end,i)/length(bin_controln6{i});
    
    nuc_output_bestn1(i) = aveArea_bin_bestn1(end,i)/length(bin_bestn1{i});
    nuc_output_bestn2(i) = aveArea_bin_bestn2(end,i)/length(bin_bestn2{i});
    nuc_output_bestn4(i) = aveArea_bin_bestn4(end,i)/length(bin_bestn4{i});
        
    nuc_output_deln1(i) = aveArea_bin_deln1(end,i)/length(bin_deln1{i});
    nuc_output_deln2(i) = aveArea_bin_deln2(end,i)/length(bin_deln2{i});
    nuc_output_deln4(i) = aveArea_bin_deln4(end,i)/length(bin_deln4{i});
    
    nuc_output_dl2bestn1(i) = aveArea_bin_dl2bestn1(end,i)/length(bin_dl2bestn1{i});
    nuc_output_dl2bestn2(i) = aveArea_bin_dl2bestn2(end,i)/length(bin_dl2bestn2{i});
    nuc_output_dl2bestn5(i) = aveArea_bin_dl2bestn5(end,i)/length(bin_dl2bestn5{i});
        
    nuc_output_dl2deln1(i) = aveArea_bin_dl2deln1(end,i)/length(bin_dl2deln1{i});
    nuc_output_dl2deln2(i) = aveArea_bin_dl2deln2(end,i)/length(bin_dl2deln2{i});
    nuc_output_dl2deln4(i) = aveArea_bin_dl2deln4(end,i)/length(bin_dl2deln4{i});
        
    nuc_output_zldgoodn2(i) = aveArea_bin_zldgoodn2(end,i)/length(bin_zldgoodn2{i});
    nuc_output_zldgoodn3(i) = aveArea_bin_zldgoodn3(end,i)/length(bin_zldgoodn3{i});
    nuc_output_zldgoodn4(i) = aveArea_bin_zldgoodn4(end,i)/length(bin_zldgoodn4{i});

    nuc_output_zldweakn1(i) = aveArea_bin_zldweakn1(end,i)/length(bin_zldweakn1{i});
    nuc_output_zldweakn2(i) = aveArea_bin_zldweakn2(end,i)/length(bin_zldweakn2{i});
    nuc_output_zldweakn5(i) = aveArea_bin_zldweakn5(end,i)/length(bin_zldweakn5{i});
    
    nuc_output_mcphetn1(i) = aveArea_bin_mcphetn1(end,i)/length(bin_mcphetn1{i});
    nuc_output_mcphetn2(i) = aveArea_bin_mcphetn2(end,i)/length(bin_mcphetn2{i});
    nuc_output_mcphetn3(i) = aveArea_bin_mcphetn3(end,i)/length(bin_mcphetn3{i});

    nuc_output_dlhetcontroln1(i) = aveArea_bin_dlhetcontroln1(end,i)/length(bin_dlhetcontroln1{i});
    nuc_output_dlhetcontroln3(i) = aveArea_bin_dlhetcontroln3(end,i)/length(bin_dlhetcontroln3{i});
    nuc_output_dlhetcontroln4(i) = aveArea_bin_dlhetcontroln4(end,i)/length(bin_dlhetcontroln4{i});

    nuc_output_dlhetdl2bestn2(i) = aveArea_bin_dlhetdl2bestn2(end,i)/length(bin_dlhetdl2bestn2{i});
    nuc_output_dlhetdl2bestn3(i) = aveArea_bin_dlhetdl2bestn3(end,i)/length(bin_dlhetdl2bestn3{i});
    nuc_output_dlhetdl2bestn5(i) = aveArea_bin_dlhetdl2bestn5(end,i)/length(bin_dlhetdl2bestn5{i});
    
end

%% Shift signal to middle
ave_Signal_controln1 = circshift(ave_Signal_controln1,0,2); 
ave_Signal_controln2 = circshift(ave_Signal_controln2,0,2);
ave_Signal_controln4 = circshift(ave_Signal_controln4,1,2);
ave_Signal_controln3 = circshift(ave_Signal_controln3,1,2);
ave_Signal_controln5 = circshift(ave_Signal_controln5,1,2);
ave_Signal_controln6 = circshift(ave_Signal_controln6,0,2);
ave_Signal_bestn1 = circshift(ave_Signal_bestn1,-1,2);
ave_Signal_bestn2 = circshift(ave_Signal_bestn2,-3,2);
ave_Signal_bestn4 = circshift(ave_Signal_bestn4,2,2);
ave_Signal_deln1 = circshift(ave_Signal_deln1,-3,2);
ave_Signal_deln2 = circshift(ave_Signal_deln2,-3,2);
ave_Signal_deln4 = circshift(ave_Signal_deln4,-1,2);
ave_Signal_dl2bestn1 = circshift(ave_Signal_dl2bestn1,0,2);
ave_Signal_dl2bestn2 = circshift(ave_Signal_dl2bestn2,1,2); %set back to 2
ave_Signal_dl2bestn5 = circshift(ave_Signal_dl2bestn5,2,2);
ave_Signal_dl2deln1 = circshift(ave_Signal_dl2deln1,0,2);
ave_Signal_dl2deln2 = circshift(ave_Signal_dl2deln2,0,2);
ave_Signal_dl2deln4 = circshift(ave_Signal_dl2deln4,1,2);
ave_Signal_zldgoodn2 = circshift(ave_Signal_zldgoodn2,2,2);
ave_Signal_zldgoodn3 = circshift(ave_Signal_zldgoodn3,2,2);
ave_Signal_zldgoodn4 = circshift(ave_Signal_zldgoodn4,2,2);
ave_Signal_zldweakn1 = circshift(ave_Signal_zldweakn1,1,2);
ave_Signal_zldweakn2 = circshift(ave_Signal_zldweakn2,1,2);
ave_Signal_zldweakn5 = circshift(ave_Signal_zldweakn5,1,2);
ave_Signal_mcphetn1 = circshift(ave_Signal_mcphetn1,1,2);
ave_Signal_mcphetn2 = circshift(ave_Signal_mcphetn2,0,2);
ave_Signal_mcphetn3 = circshift(ave_Signal_mcphetn3,-1,2);
ave_Signal_dlhetcontroln1 = circshift(ave_Signal_dlhetcontroln1,1,2);
ave_Signal_dlhetcontroln3 = circshift(ave_Signal_dlhetcontroln3,-1,2);
ave_Signal_dlhetcontroln4 = circshift(ave_Signal_dlhetcontroln4,0,2);
ave_Signal_dlhetdl2bestn2 = circshift(ave_Signal_dlhetdl2bestn2,1,2);
ave_Signal_dlhetdl2bestn3 = circshift(ave_Signal_dlhetdl2bestn3,0,2);
ave_Signal_dlhetdl2bestn5 = circshift(ave_Signal_dlhetdl2bestn5,1,2);
ave_Signal_twideln1 = circshift(ave_Signal_twideln1,-1,2);
ave_Signal_twideln2 = circshift(ave_Signal_twideln2,1,2);

ave_Area_controln1 = circshift(ave_Area_controln1,0,2); 
ave_Area_controln2 = circshift(ave_Area_controln2,0,2);
ave_Area_controln4 = circshift(ave_Area_controln4,1,2);
ave_Area_controln3 = circshift(ave_Area_controln3,1,2);
ave_Area_controln5 = circshift(ave_Area_controln5,1,2);
ave_Area_controln6 = circshift(ave_Area_controln6,0,2);
ave_Area_bestn1 = circshift(ave_Area_bestn1,-1,2);
ave_Area_bestn2 = circshift(ave_Area_bestn2,-3,2);
ave_Area_bestn4 = circshift(ave_Area_bestn4,2,2);
ave_Area_deln1 = circshift(ave_Area_deln1,-3,2);
ave_Area_deln2 = circshift(ave_Area_deln2,-3,2);
ave_Area_deln4 = circshift(ave_Area_deln4,-1,2);
ave_Area_dl2bestn1 = circshift(ave_Area_dl2bestn1,0,2);
ave_Area_dl2bestn2 = circshift(ave_Area_dl2bestn2,1,2); %set back to 2
ave_Area_dl2bestn5 = circshift(ave_Area_dl2bestn5,2,2);
ave_Area_dl2deln1 = circshift(ave_Area_dl2deln1,0,2);
ave_Area_dl2deln2 = circshift(ave_Area_dl2deln2,0,2);
ave_Area_dl2deln4 = circshift(ave_Area_dl2deln4,1,2);
ave_Area_zldgoodn2 = circshift(ave_Area_zldgoodn2,2,2);
ave_Area_zldgoodn3 = circshift(ave_Area_zldgoodn3,2,2);
ave_Area_zldgoodn4 = circshift(ave_Area_zldgoodn4,2,2);
ave_Area_zldweakn1 = circshift(ave_Area_zldweakn1,1,2);
ave_Area_zldweakn2 = circshift(ave_Area_zldweakn2,1,2);
ave_Area_zldweakn5 = circshift(ave_Area_zldweakn5,1,2);
ave_Area_mcphetn1 = circshift(ave_Area_mcphetn1,1,2);
ave_Area_mcphetn2 = circshift(ave_Area_mcphetn2,0,2);
ave_Area_mcphetn3 = circshift(ave_Area_mcphetn3,-1,2);
ave_Area_dlhetcontroln1 = circshift(ave_Area_dlhetcontroln1,1,2);
ave_Area_dlhetcontroln3 = circshift(ave_Area_dlhetcontroln3,-3,2);
ave_Area_dlhetcontroln4 = circshift(ave_Area_dlhetcontroln4,0,2);
ave_Area_dlhetdl2bestn2 = circshift(ave_Area_dlhetdl2bestn2,1,2);
ave_Area_dlhetdl2bestn3 = circshift(ave_Area_dlhetdl2bestn3,0,2);
ave_Area_dlhetdl2bestn5 = circshift(ave_Area_dlhetdl2bestn5,0,2);
ave_Area_twideln1 = circshift(ave_Area_twideln1,-1,2);
ave_Area_twideln2 = circshift(ave_Area_twideln2,1,2);

nuc_output_controln1 = circshift(nuc_output_controln1,-1,2); 
nuc_output_controln2 = circshift(nuc_output_controln2,0,2);
nuc_output_controln4 = circshift(nuc_output_controln4,1,2);
nuc_output_controln3 = circshift(nuc_output_controln3,1,2);
nuc_output_controln5 = circshift(nuc_output_controln5,1,2);
nuc_output_controln6 = circshift(nuc_output_controln6,0,2);
nuc_output_bestn1 = circshift(nuc_output_bestn1,-2,2);
nuc_output_bestn2 = circshift(nuc_output_bestn2,-1,2);
nuc_output_bestn4 = circshift(nuc_output_bestn4,3,2);
nuc_output_deln1 = circshift(nuc_output_deln1,-2,2);
nuc_output_deln2 = circshift(nuc_output_deln2,0,2);
nuc_output_deln4 = circshift(nuc_output_deln4,-1,2);
nuc_output_dl2bestn1 = circshift(nuc_output_dl2bestn1,0,2);
nuc_output_dl2bestn2 = circshift(nuc_output_dl2bestn2,2,2);
nuc_output_dl2bestn5 = circshift(nuc_output_dl2bestn5,2,2);
nuc_output_dl2deln1 = circshift(nuc_output_dl2deln1,1,2);
nuc_output_dl2deln2 = circshift(nuc_output_dl2deln2,0,2);
nuc_output_dl2deln4 = circshift(nuc_output_dl2deln4,0,2);
nuc_output_zldgoodn2 = circshift(nuc_output_zldgoodn2,2,2);
nuc_output_zldgoodn3 = circshift(nuc_output_zldgoodn3,1,2);
nuc_output_zldgoodn4 = circshift(nuc_output_zldgoodn4,2,2);
nuc_output_zldweakn1 = circshift(nuc_output_zldweakn1,0,2);
nuc_output_zldweakn2 = circshift(nuc_output_zldweakn2,1,2);
nuc_output_zldweakn5 = circshift(nuc_output_zldweakn5,1,2);
nuc_output_mcphetn1 = circshift(nuc_output_mcphetn1,1,2);
nuc_output_mcphetn2 = circshift(nuc_output_mcphetn2,0,2);
nuc_output_mcphetn3 = circshift(nuc_output_mcphetn3,-1,2);
nuc_output_dlhetcontroln1 = circshift(nuc_output_dlhetcontroln1,1,2);
nuc_output_dlhetcontroln3 = circshift(nuc_output_dlhetcontroln3,-3,2);
nuc_output_dlhetcontroln4 = circshift(nuc_output_dlhetcontroln4,0,2);
nuc_output_dlhetdl2bestn2 = circshift(nuc_output_dlhetdl2bestn2,1,2);
nuc_output_dlhetdl2bestn3 = circshift(nuc_output_dlhetdl2bestn3,0,2);
nuc_output_dlhetdl2bestn5 = circshift(nuc_output_dlhetdl2bestn5,0,2);

Area_bin_cn1 = circshift(Area_bin_cn1,0,2);
Area_bin_cn2 = circshift(Area_bin_cn2,0,2);
Area_bin_cn4 = circshift(Area_bin_cn4,1,2);
Area_bin_cn3 = circshift(Area_bin_cn3,1,2);
Area_bin_cn5 = circshift(Area_bin_cn5,1,2);
Area_bin_cn6 = circshift(Area_bin_cn6,0,2);
Area_bin_bn1 = circshift(Area_bin_bn1,-1,2);
Area_bin_bn2 = circshift(Area_bin_bn2,-3,2);
Area_bin_bn4 = circshift(Area_bin_bn4,2,2);
Area_bin_dn1 = circshift(Area_bin_dn1,-3,2);
Area_bin_dn2 = circshift(Area_bin_dn2,-3,2);
Area_bin_dn4 = circshift(Area_bin_dn4,-1,2);
Area_bin_dl2bn1 = circshift(Area_bin_dl2bn1,0,2);
Area_bin_dl2bn2 = circshift(Area_bin_dl2bn2,1,2); %set back to 2
Area_bin_dl2bn5 = circshift(Area_bin_dl2bn5,2,2);
Area_bin_dl2dn1 = circshift(Area_bin_dl2dn1,0,2);
Area_bin_dl2dn2 = circshift(Area_bin_dl2dn2,0,2);
Area_bin_dl2dn4 = circshift(Area_bin_dl2dn4,1,2);
Area_bin_zgn2 = circshift(Area_bin_zgn2,2,2);
Area_bin_zgn3 = circshift(Area_bin_zgn3,2,2);
Area_bin_zgn4 = circshift(Area_bin_zgn4,2,2);
Area_bin_tdn1 = circshift(Area_bin_tdn1,-1,2);
Area_bin_tdn2 = circshift(Area_bin_tdn2,1,2);

%het 
ave_Signal_hetcontroln1 = circshift(ave_Signal_hetcontroln1,-2,2);
ave_Signal_hetcontroln2 = circshift(ave_Signal_hetcontroln2,-1,2);
ave_Signal_hetcontroln3 = circshift(ave_Signal_hetcontroln3,0,2);
ave_Signal_hetbestn1 = circshift(ave_Signal_hetbestn1,2,2);
ave_Signal_hetbestn2 = circshift(ave_Signal_hetbestn2,2,2);
ave_Signal_hetbestn3 = circshift(ave_Signal_hetbestn3,0,2);
ave_Signal_hetdl2bestn1 = circshift(ave_Signal_hetdl2bestn1,1,2);
ave_Signal_hetdl2bestn2 = circshift(ave_Signal_hetdl2bestn2,3,2);
ave_Signal_hetdl2bestn3 = circshift(ave_Signal_hetdl2bestn3,0,2);
ave_Signal_hetdl2deln1 = circshift(ave_Signal_hetdl2deln1,0,2);
ave_Signal_hetdl2deln2 = circshift(ave_Signal_hetdl2deln2,1,2);
ave_Signal_hetdl2deln3 = circshift(ave_Signal_hetdl2deln3,2,2);
ave_Signal_hetzldbestn1 = circshift(ave_Signal_hetzldbestn1,-1,2);
ave_Signal_hetzldbestn2 = circshift(ave_Signal_hetzldbestn2,0,2);
ave_Signal_hetzldbestn3 = circshift(ave_Signal_hetzldbestn3,2,2);

ave_Area_hetcontroln1 = circshift(ave_Area_hetcontroln1,-2,2);
ave_Area_hetcontroln2 = circshift(ave_Area_hetcontroln2,-1,2);
ave_Area_hetcontroln3 = circshift(ave_Area_hetcontroln3,0,2);
ave_Area_hetbestn1 = circshift(ave_Area_hetbestn1,2,2);
ave_Area_hetbestn2 = circshift(ave_Area_hetbestn2,2,2);
ave_Area_hetbestn3 = circshift(ave_Area_hetbestn3,0,2);
ave_Area_hetdl2bestn1 = circshift(ave_Area_hetdl2bestn1,1,2);
ave_Area_hetdl2bestn2 = circshift(ave_Area_hetdl2bestn2,3,2);
ave_Area_hetdl2bestn3 = circshift(ave_Area_hetdl2bestn3,0,2);
ave_Area_hetdl2deln1 = circshift(ave_Area_hetdl2deln1,0,2);
ave_Area_hetdl2deln2 = circshift(ave_Area_hetdl2deln2,1,2);
ave_Area_hetdl2deln3 = circshift(ave_Area_hetdl2deln3,2,2);
ave_Area_hetzldbestn1 = circshift(ave_Area_hetzldbestn1,-1,2);
ave_Area_hetzldbestn2 = circshift(ave_Area_hetzldbestn2,0,2);
ave_Area_hetzldbestn3 = circshift(ave_Area_hetzldbestn3,2,2);

%Revision
ave_Signal_Rcontroln1 = circshift(ave_Signal_Rcontroln1,1,2); 
ave_Signal_Rcontroln2 = circshift(ave_Signal_Rcontroln2,-1,2);
ave_Signal_Rcontroln3 = circshift(ave_Signal_Rcontroln3,2,2);
ave_Signal_Rdzn1 = circshift(ave_Signal_Rdzn1,2,2); 
ave_Signal_Rdzn2 = circshift(ave_Signal_Rdzn2,2,2);
ave_Signal_Rdzn3 = circshift(ave_Signal_Rdzn3,-2,2);
ave_Signal_Rdn1 = circshift(ave_Signal_Rdn1,2,2); 
ave_Signal_Rdn2 = circshift(ave_Signal_Rdn2,2,2);
ave_Signal_Rdn3 = circshift(ave_Signal_Rdn3,0,2);

ave_Area_Rcontroln1 = circshift(ave_Area_Rcontroln1,1,2); 
ave_Area_Rcontroln2 = circshift(ave_Area_Rcontroln2,-1,2);
ave_Area_Rcontroln3 = circshift(ave_Area_Rcontroln3,2,2);
ave_Area_Rdzn1 = circshift(ave_Area_Rdzn1,2,2); 
ave_Area_Rdzn2 = circshift(ave_Area_Rdzn2,2,2);
ave_Area_Rdzn3 = circshift(ave_Area_Rdzn3,-2,2);
ave_Area_Rdn1 = circshift(ave_Area_Rdn1,2,2); 
ave_Area_Rdn2 = circshift(ave_Area_Rdn2,2,2);
ave_Area_Rdn3 = circshift(ave_Area_Rdn3,0,2);

Area_bin_Rcn1 = circshift(Area_bin_Rcn1,1,2);
Area_bin_Rcn2 = circshift(Area_bin_Rcn2,-1,2);
Area_bin_Rcn3 = circshift(Area_bin_Rcn3,2,2);
Area_bin_Rdzn1 = circshift(Area_bin_Rdzn1,2,2);
Area_bin_Rdzn2 = circshift(Area_bin_Rdzn2,2,2);
Area_bin_Rdzn3 = circshift(Area_bin_Rdzn3,-2,2);
Area_bin_Rdn1 = circshift(Area_bin_Rdn1,2,2);
Area_bin_Rdn2 = circshift(Area_bin_Rdn2,2,2);
Area_bin_Rdn3 = circshift(Area_bin_Rdn3,0,2);

%Standard Error
stdSignal_bin_controln1 = circshift(stdSignal_bin_controln1,0,2); 
stdSignal_bin_controln2 = circshift(stdSignal_bin_controln2,0,2);
stdSignal_bin_controln4 = circshift(stdSignal_bin_controln4,0,2);
stdSignal_bin_controln3 = circshift(stdSignal_bin_controln3,1,2);
stdSignal_bin_controln5 = circshift(stdSignal_bin_controln5,1,2);
stdSignal_bin_controln6 = circshift(stdSignal_bin_controln6,0,2);
stdSignal_bin_bestn1 = circshift(stdSignal_bin_bestn1,-1,2);
stdSignal_bin_bestn2 = circshift(stdSignal_bin_bestn2,-3,2);
stdSignal_bin_bestn4 = circshift(stdSignal_bin_bestn4,2,2);
stdSignal_bin_deln1 = circshift(stdSignal_bin_deln1,-3,2);
stdSignal_bin_deln2 = circshift(stdSignal_bin_deln2,-3,2);
stdSignal_bin_deln4 = circshift(stdSignal_bin_deln4,-1,2);
stdSignal_bin_dl2bestn1 = circshift(stdSignal_bin_dl2bestn1,0,2);
stdSignal_bin_dl2bestn2 = circshift(stdSignal_bin_dl2bestn2,1,2); %set back to 2
stdSignal_bin_dl2bestn5 = circshift(stdSignal_bin_dl2bestn5,2,2);
stdSignal_bin_dl2deln1 = circshift(stdSignal_bin_dl2deln1,0,2);
stdSignal_bin_dl2deln2 = circshift(stdSignal_bin_dl2deln2,0,2);
stdSignal_bin_dl2deln4 = circshift(stdSignal_bin_dl2deln4,1,2);
stdSignal_bin_zldgoodn2 = circshift(stdSignal_bin_zldgoodn2,2,2);
stdSignal_bin_zldgoodn3 = circshift(stdSignal_bin_zldgoodn3,2,2);
stdSignal_bin_zldgoodn4 = circshift(stdSignal_bin_zldgoodn4,2,2);
stdSignal_bin_zldweakn1 = circshift(stdSignal_bin_zldweakn1,1,2);
stdSignal_bin_zldweakn2 = circshift(stdSignal_bin_zldweakn2,1,2);
stdSignal_bin_zldweakn5 = circshift(stdSignal_bin_zldweakn5,1,2);
stdSignal_mcphetn1 = circshift(stdSignal_bin_mcphetn1,1,2);
stdSignal_mcphetn2 = circshift(stdSignal_bin_mcphetn2,0,2);
stdSignal_mcphetn3 = circshift(stdSignal_bin_mcphetn3,-1,2);
stdSignal_dlhetcontroln1 = circshift(stdSignal_bin_dlhetcontroln1,1,2);
stdSignal_dlhetcontroln3 = circshift(stdSignal_bin_dlhetcontroln3,-1,2);
stdSignal_dlhetcontroln4 = circshift(stdSignal_bin_dlhetcontroln4,0,2);
stdSignal_dlhetdl2bestn2 = circshift(stdSignal_bin_dlhetdl2bestn2,1,2);
stdSignal_dlhetdl2bestn3 = circshift(stdSignal_bin_dlhetdl2bestn3,0,2);
stdSignal_dlhetdl2bestn5 = circshift(stdSignal_bin_dlhetdl2bestn5,1,2);
stdSignal_bin_twideln1 = circshift(stdSignal_bin_twideln1,-1,2);
stdSignal_bin_twideln2 = circshift(stdSignal_bin_twideln2,1,2);

stdArea_bin_controln1 = circshift(stdArea_bin_controln1,0,2); 
stdArea_bin_controln2 = circshift(stdArea_bin_controln2,0,2);
stdArea_bin_controln4 = circshift(stdArea_bin_controln4,1,2);
stdArea_bin_controln3 = circshift(stdArea_bin_controln3,1,2);
stdArea_bin_controln5 = circshift(stdArea_bin_controln5,1,2);
stdArea_bin_controln6 = circshift(stdArea_bin_controln6,0,2);
stdArea_bin_bestn1 = circshift(stdArea_bin_bestn1,-1,2);
stdArea_bin_bestn2 = circshift(stdArea_bin_bestn2,-3,2);
stdArea_bin_bestn4 = circshift(stdArea_bin_bestn4,2,2);
stdArea_bin_deln1 = circshift(stdArea_bin_deln1,-3,2);
stdArea_bin_deln2 = circshift(stdArea_bin_deln2,-3,2);
stdArea_bin_deln4 = circshift(stdArea_bin_deln4,-1,2);
stdArea_bin_dl2bestn1 = circshift(stdArea_bin_dl2bestn1,0,2);
stdArea_bin_dl2bestn2 = circshift(stdArea_bin_dl2bestn2,1,2); %set back to 2
stdArea_bin_dl2bestn5 = circshift(stdArea_bin_dl2bestn5,2,2);
stdArea_bin_dl2deln1 = circshift(stdArea_bin_dl2deln1,0,2);
stdArea_bin_dl2deln2 = circshift(stdArea_bin_dl2deln2,0,2);
stdArea_bin_dl2deln4 = circshift(stdArea_bin_dl2deln4,1,2);
stdArea_bin_zldgoodn2 = circshift(stdArea_bin_zldgoodn2,2,2);
stdArea_bin_zldgoodn3 = circshift(stdArea_bin_zldgoodn3,2,2);
stdArea_bin_zldgoodn4 = circshift(stdArea_bin_zldgoodn4,2,2);
stdArea_bin_zldweakn1 = circshift(stdArea_bin_zldweakn1,1,2);
stdArea_bin_zldweakn2 = circshift(stdArea_bin_zldweakn2,1,2);
stdArea_bin_zldweakn5 = circshift(stdArea_bin_zldweakn5,1,2);
stdArea_mcphetn1 = circshift(stdArea_bin_mcphetn1,1,2);
stdArea_mcphetn2 = circshift(stdArea_bin_mcphetn2,0,2);
stdArea_mcphetn3 = circshift(stdArea_bin_mcphetn3,-1,2);
stdArea_dlhetcontroln1 = circshift(stdArea_bin_dlhetcontroln1,1,2);
stdArea_dlhetcontroln3 = circshift(stdArea_bin_dlhetcontroln3,-3,2);
stdArea_dlhetcontroln4 = circshift(stdArea_bin_dlhetcontroln4,0,2);
stdArea_dlhetdl2bestn2 = circshift(stdArea_bin_dlhetdl2bestn2,1,2);
stdArea_dlhetdl2bestn3 = circshift(stdArea_bin_dlhetdl2bestn3,0,2);
stdArea_dlhetdl2bestn5 = circshift(stdArea_bin_dlhetdl2bestn5,0,2);
stdArea_bin_twideln1 = circshift(stdArea_bin_twideln1,-1,2);
stdArea_bin_twideln2 = circshift(stdArea_bin_twideln2,1,2);

%het
stdSignal_bin_hetcontroln1 = circshift(stdSignal_bin_hetcontroln1,-2,2);
stdSignal_bin_hetcontroln2 = circshift(stdSignal_bin_hetcontroln2,-1,2);
stdSignal_bin_hetcontroln3 = circshift(stdSignal_bin_hetcontroln3,0,2);
stdSignal_bin_hetbestn1 = circshift(stdSignal_bin_hetbestn1,2,2);
stdSignal_bin_hetbestn2 = circshift(stdSignal_bin_hetbestn2,2,2);
stdSignal_bin_hetbestn3 = circshift(stdSignal_bin_hetbestn3,0,2);
stdSignal_bin_hetdl2bestn1 = circshift(stdSignal_bin_hetdl2bestn1,1,2);
stdSignal_bin_hetdl2bestn2 = circshift(stdSignal_bin_hetdl2bestn2,3,2);
stdSignal_bin_hetdl2bestn3 = circshift(stdSignal_bin_hetdl2bestn3,0,2);
stdSignal_bin_hetdl2deln1 = circshift(stdSignal_bin_hetdl2deln1,0,2);
stdSignal_bin_hetdl2deln2 = circshift(stdSignal_bin_hetdl2deln2,1,2);
stdSignal_bin_hetdl2deln3 = circshift(stdSignal_bin_hetdl2deln3,2,2);
stdSignal_bin_hetzldbestn1 = circshift(stdSignal_bin_hetzldbestn1,-1,2);
stdSignal_bin_hetzldbestn2 = circshift(stdSignal_bin_hetzldbestn2,0,2);
stdSignal_bin_hetzldbestn3 = circshift(stdSignal_bin_hetzldbestn3,2,2);

stdArea_bin_hetcontroln1 = circshift(stdArea_bin_hetcontroln1,-2,2);
stdArea_bin_hetcontroln2 = circshift(stdArea_bin_hetcontroln2,-1,2);
stdArea_bin_hetcontroln3 = circshift(stdArea_bin_hetcontroln3,0,2);
stdArea_bin_hetbestn1 = circshift(stdArea_bin_hetbestn1,2,2);
stdArea_bin_hetbestn2 = circshift(stdArea_bin_hetbestn2,2,2);
stdArea_bin_hetbestn3 = circshift(stdArea_bin_hetbestn3,0,2);
stdArea_bin_hetdl2bestn1 = circshift(stdArea_bin_hetdl2bestn1,1,2);
stdArea_bin_hetdl2bestn2 = circshift(stdArea_bin_hetdl2bestn2,3,2);
stdArea_bin_hetdl2bestn3 = circshift(stdArea_bin_hetdl2bestn3,0,2);
stdArea_bin_hetdl2deln1 = circshift(stdArea_bin_hetdl2deln1,0,2);
stdArea_bin_hetdl2deln2 = circshift(stdArea_bin_hetdl2deln2,1,2);
stdArea_bin_hetdl2deln3 = circshift(stdArea_bin_hetdl2deln3,2,2);
stdArea_bin_hetzldbestn1 = circshift(stdArea_bin_hetzldbestn1,-1,2);
stdArea_bin_hetzldbestn2 = circshift(stdArea_bin_hetzldbestn2,0,2);
stdArea_bin_hetzldbestn3 = circshift(stdArea_bin_hetzldbestn3,2,2);

%Revision Std Dv
stdSignal_bin_Rcn1 = circshift(stdSignal_bin_Rcn1,1,2); 
stdSignal_bin_Rcn2 = circshift(stdSignal_bin_Rcn2,-1,2);
stdSignal_bin_Rcn3 = circshift(stdSignal_bin_Rcn3,2,2);
stdSignal_bin_Rdzn1 = circshift(stdSignal_bin_Rdzn1,2,2); 
stdSignal_bin_Rdzn2 = circshift(stdSignal_bin_Rdzn2,2,2);
stdSignal_bin_Rdzn3 = circshift(stdSignal_bin_Rdzn3,-2,2);
stdSignal_bin_Rdn1 = circshift(stdSignal_bin_Rdn1,2,2); 
stdSignal_bin_Rdn2 = circshift(stdSignal_bin_Rdn2,2,2);
stdSignal_bin_Rdn3 = circshift(stdSignal_bin_Rdn3,0,2);

stdArea_bin_Rcn1 = circshift(stdArea_bin_Rcn1,1,2); 
stdArea_bin_Rcn2 = circshift(stdArea_bin_Rcn2,-1,2);
stdArea_bin_Rcn3 = circshift(stdArea_bin_Rcn3,2,2);
stdArea_bin_Rdzn1 = circshift(stdArea_bin_Rdzn1,2,2); 
stdArea_bin_Rdzn2 = circshift(stdArea_bin_Rdzn2,2,2);
stdArea_bin_Rdzn3 = circshift(stdArea_bin_Rdzn3,-2,2);
stdArea_bin_Rdn1 = circshift(stdArea_bin_Rdn1,2,2); 
stdArea_bin_Rdn2 = circshift(stdArea_bin_Rdn2,2,2);
stdArea_bin_Rdn3 = circshift(stdArea_bin_Rdn3,0,2);

%Bins
shift_bin_controln1 = circshift(bin_controln1,0,2); 
shift_bin_controln2 = circshift(bin_controln2,0,2);
shift_bin_controln4 = circshift(bin_controln4,1,2);
shift_bin_controln3 = circshift(bin_controln3,1,2);
shift_bin_controln5 = circshift(bin_controln5,1,2);
shift_bin_controln6 = circshift(bin_controln6,0,2);
shift_bin_bestn1 = circshift(bin_bestn1,-1,2);
shift_bin_bestn2 = circshift(bin_bestn2,-3,2);
shift_bin_bestn4 = circshift(bin_bestn4,2,2);
shift_bin_deln1 = circshift(bin_deln1,-3,2);
shift_bin_deln2 = circshift(bin_deln2,-3,2);
shift_bin_deln4 = circshift(bin_deln4,-1,2);
shift_bin_dl2bestn1 = circshift(bin_dl2bestn1,0,2);
shift_bin_dl2bestn2 = circshift(bin_dl2bestn2,1,2); %set back to 2
shift_bin_dl2bestn5 = circshift(bin_dl2bestn5,2,2);
shift_bin_dl2deln1 = circshift(bin_dl2deln1,0,2);
shift_bin_dl2deln2 = circshift(bin_dl2deln2,0,2);
shift_bin_dl2deln4 = circshift(bin_dl2deln4,1,2);
shift_bin_zldgoodn2 = circshift(bin_zldgoodn2,2,2);
shift_bin_zldgoodn3 = circshift(bin_zldgoodn3,2,2);
shift_bin_zldgoodn4 = circshift(bin_zldgoodn4,2,2);
shift_bin_zldweakn1 = circshift(bin_zldweakn1,1,2);
shift_bin_zldweakn2 = circshift(bin_zldweakn2,1,2);
shift_bin_zldweakn5 = circshift(bin_zldweakn5,1,2);
shift_bin_mcphetn1 = circshift(bin_mcphetn1,1,2);
shift_bin_mcphetn2 = circshift(bin_mcphetn2,0,2);
shift_bin_mcphetn3 = circshift(bin_mcphetn3,-1,2);
shift_bin_dlhetcontroln1 = circshift(bin_dlhetcontroln1,1,2);
shift_bin_dlhetcontroln3 = circshift(bin_dlhetcontroln3,-1,2);
shift_bin_dlhetcontroln4 = circshift(bin_dlhetcontroln4,0,2);
shift_bin_dlhetdl2bestn2 = circshift(bin_dlhetdl2bestn2,1,2);
shift_bin_dlhetdl2bestn3 = circshift(bin_dlhetdl2bestn3,0,2);
shift_bin_dlhetdl2bestn5 = circshift(bin_dlhetdl2bestn5,1,2);
shift_bin_twideln1 = circshift(bin_twideln1,-1,2);
shift_bin_twideln2 = circshift(bin_twideln2,1,2);

%het
shift_bin_hetcontroln1 = circshift(bin_hetcontroln1,-2,2);
shift_bin_hetcontroln2 = circshift(bin_hetcontroln2,-1,2);
shift_bin_hetcontroln3 = circshift(bin_hetcontroln3,0,2);
shift_bin_hetbestn1 = circshift(bin_hetbestn1,2,2);
shift_bin_hetbestn2 = circshift(bin_hetbestn2,2,2);
shift_bin_hetbestn3 = circshift(bin_hetbestn3,0,2);
shift_bin_hetdl2bestn1 = circshift(bin_hetdl2bestn1,1,2);
shift_bin_hetdl2bestn2 = circshift(bin_hetdl2bestn2,3,2);
shift_bin_hetdl2bestn3 = circshift(bin_hetdl2bestn3,0,2);
shift_bin_hetdl2deln1 = circshift(bin_hetdl2deln1,0,2);
shift_bin_hetdl2deln2 = circshift(bin_hetdl2deln2,1,2);
shift_bin_hetdl2deln3 = circshift(bin_hetdl2deln3,2,2);
shift_bin_hetzldbestn1 = circshift(bin_hetzldbestn1,-1,2);
shift_bin_hetzldbestn2 = circshift(bin_hetzldbestn2,0,2);
shift_bin_hetzldbestn3 = circshift(bin_hetzldbestn3,2,2);

%Revision
shift_bin_Rcn1 = circshift(bin_Rcn1,1,2); 
shift_bin_Rcn2 = circshift(bin_Rcn2,-1,2);
shift_bin_Rcn3 = circshift(bin_Rcn3,2,2);
shift_bin_Rdzn1 = circshift(bin_Rdzn1,2,2); 
shift_bin_Rdzn2 = circshift(bin_Rdzn2,2,2);
shift_bin_Rdzn3 = circshift(bin_Rdzn3,-2,2);
shift_bin_Rdn1 = circshift(bin_Rdn1,2,2); 
shift_bin_Rdn2 = circshift(bin_Rdn2,2,2);
shift_bin_Rdn3 = circshift(bin_Rdn3,0,2);

%% Signal Plot
point = 150;
% figure, plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_controln1(point,:)); hold on
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_controln2(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_controln4(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_controln3(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_controln5(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_controln6(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_bestn1(point,:));hold on
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_bestn2(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_bestn4(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_deln1(point,:)); hold on
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_deln2(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_deln4(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2bestn1(point,:));hold on
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2bestn2(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2bestn5(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2deln1(point,:));hold on
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2deln2(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2deln4(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldgoodn2(point,:));hold on
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldgoodn3(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldgoodn4(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldweakn1(point,:));hold on
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldweakn2(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldweakn5(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_mcphetn1(point,:));hold on
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_mcphetn2(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_mcphetn3(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dlhetcontroln1(point,:));hold on
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dlhetcontroln3(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dlhetcontroln4(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dlhetdl2bestn2(point,:));hold on
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dlhetdl2bestn3(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dlhetdl2bestn5(point,:))
% figure, plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_Rcontroln1(point,:)); hold on
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_Rcontroln2(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_Rcontroln3(point,:))
% figure, plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_Rdzn1(point,:)); hold on
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_Rdzn2(point,:))
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_Rdzn3(point,:))

% title('Signal Output Across Embryo')
% xlabel('Location')
% ylabel('Nuclei Signal')


ave_Signal_control = (ave_Signal_controln1 + ave_Signal_controln2 + ave_Signal_controln4 + ave_Signal_controln3 + ave_Signal_controln5 + ave_Signal_controln6)./6;
ave_Signal_best = (ave_Signal_bestn1 + ave_Signal_bestn2 + ave_Signal_bestn4)./rep;
ave_Signal_del = (ave_Signal_deln1 + ave_Signal_deln2 + ave_Signal_deln4)./rep;

ave_Signal_dl2best = (ave_Signal_dl2bestn1 + ave_Signal_dl2bestn2 + ave_Signal_dl2bestn5)./rep;
ave_Signal_dl2del = (ave_Signal_dl2deln1 + ave_Signal_dl2deln2 + ave_Signal_dl2deln4)./rep;

ave_Signal_zldgood = (ave_Signal_zldgoodn2 + ave_Signal_zldgoodn3 + ave_Signal_zldgoodn4)./rep;
ave_Signal_zldweak = (ave_Signal_zldweakn1 + ave_Signal_zldweakn2 + ave_Signal_zldweakn5)./rep;

% ave_Signal_mcphet = (ave_Signal_mcphetn1 + ave_Signal_mcphetn2 + ave_Signal_mcphetn3)./rep;
ave_Signal_mcphet = (ave_Signal_mcphetn2 + ave_Signal_mcphetn3)./2;

ave_Signal_dlhetcontrol = (ave_Signal_dlhetcontroln1 + ave_Signal_dlhetcontroln3 + ave_Signal_dlhetcontroln4)./rep;
ave_Signal_dlhetdl2best = (ave_Signal_dlhetdl2bestn2 + ave_Signal_dlhetdl2bestn3 + ave_Signal_dlhetdl2bestn5)./rep;


%het
ave_Signal_hetcontrol = (ave_Signal_hetcontroln1 + ave_Signal_hetcontroln2 + ave_Signal_hetcontroln3)./rep;
ave_Signal_hetbest = (ave_Signal_hetbestn1 + ave_Signal_hetbestn2 + ave_Signal_hetbestn3)./rep;

ave_Signal_hetdl2best = (ave_Signal_hetdl2bestn1 + ave_Signal_hetdl2bestn2 + ave_Signal_hetdl2bestn3)./rep;
ave_Signal_hetdl2del = (ave_Signal_hetdl2deln1 + ave_Signal_hetdl2deln2 + ave_Signal_hetdl2deln3)./rep;

ave_Signal_hetzldbest = (ave_Signal_hetzldbestn1 + ave_Signal_hetzldbestn3 + ave_Signal_hetzldbestn3)./rep;

ave_Signal_Rcontrol = (ave_Signal_Rcontroln1 + ave_Signal_Rcontroln2 + ave_Signal_Rcontroln3)./3;

%calculate standard deviation
for i = 1:length(norm_time)
    for j = 1:bins_bestn1
        std_Signal_control(i,j) = std([ave_Signal_controln1(i,j) ave_Signal_controln2(i,j) ave_Signal_controln4(i,j) ave_Signal_controln3(i,j) ave_Signal_controln5(i,j) ave_Signal_controln6(i,j)]);
        std_Signal_best(i,j) = std([ave_Signal_bestn1(i,j) ave_Signal_bestn2(i,j) ave_Signal_bestn4(i,j)]);
        std_Signal_del(i,j) = std([ave_Signal_deln1(i,j) ave_Signal_deln2(i,j) ave_Signal_deln4(i,j)]);
        std_Signal_dl2best(i,j) = std([ave_Signal_dl2bestn1(i,j) ave_Signal_dl2bestn2(i,j) ave_Signal_dl2bestn5(i,j)]);
        std_Signal_dl2del(i,j) = std([ave_Signal_dl2deln1(i,j) ave_Signal_dl2deln2(i,j) ave_Signal_dl2deln4(i,j)]);
        std_Signal_zldgood(i,j) = std([ave_Signal_zldgoodn2(i,j) ave_Signal_zldgoodn3(i,j) ave_Signal_zldgoodn4(i,j)]);
        std_Signal_zldweak(i,j) = std([ave_Signal_zldweakn1(i,j) ave_Signal_zldweakn2(i,j) ave_Signal_zldweakn5(i,j)]);
%         std_Signal_mcphet(i,j) = std([ave_Signal_mcphetn1(i,j) ave_Signal_mcphetn2(i,j) ave_Signal_mcphetn3(i,j)]);
        std_Signal_mcphet(i,j) = std([ave_Signal_mcphetn2(i,j) ave_Signal_mcphetn3(i,j)]);
        std_Signal_dlhetcontrol(i,j) = std([ave_Signal_dlhetcontroln1(i,j) ave_Signal_dlhetcontroln3(i,j) ave_Signal_dlhetcontroln4(i,j)]);
        std_Signal_dlhetdl2best(i,j) = std([ave_Signal_dlhetdl2bestn2(i,j) ave_Signal_dlhetdl2bestn3(i,j) ave_Signal_dlhetdl2bestn5(i,j)]);
        
        std_Signal_Rcontrol(i,j) = std([ave_Signal_controln1(i,j) ave_Signal_controln2(i,j) ave_Signal_controln3(i,j)]);
    end
end

%het
for i = 1:length(norm_time)
    for j = 1:bins_hetbestn1
        std_Signal_hetcontrol(i,j) = std([ave_Signal_hetcontroln1(i,j) ave_Signal_hetcontroln2(i,j) ave_Signal_hetcontroln3(i,j)]);
        std_Signal_hetbest(i,j) = std([ave_Signal_hetbestn1(i,j) ave_Signal_hetbestn2(i,j) ave_Signal_hetbestn3(i,j)]);
        std_Signal_hetdl2best(i,j) = std([ave_Signal_hetdl2bestn1(i,j) ave_Signal_hetdl2bestn2(i,j) ave_Signal_hetdl2bestn3(i,j)]);
        std_Signal_hetdl2del(i,j) = std([ave_Signal_hetdl2deln1(i,j) ave_Signal_hetdl2deln2(i,j) ave_Signal_hetdl2deln3(i,j)]);
        std_Signal_hetzldbest(i,j) = std([ave_Signal_hetzldbestn1(i,j) ave_Signal_hetzldbestn2(i,j) ave_Signal_hetzldbestn3(i,j)]);
    end
end

% ave_Signal_control(:,[1,2,16]) = NaN;
ave_Signal_del(:,15) = NaN;
ave_Signal_dl2del(:,[1,2,16]) = NaN;

for i = 1:181
ave_combo_del(i,:) = nansum([ave_Signal_del(i,:)-170;ave_Signal_dl2del(i,:)-170]);
std_combo_del(i,:) = nansum([std_Signal_del(i,:);std_Signal_dl2del(i,:)]);
end
%{        
for i = 1:length(norm_time)
    h = figure;
    errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_control(i,:)-170,std_Signal_control(i,:),'color',rgb('black'),'LineWidth',2); hold on
    %errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_best(i,:),std_Signal_best(i,:),'--','color',rgb('red')); 
    errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_del(i,:)-170,std_Signal_del(i,:),'color',rgb('lightcoral'),'LineWidth',2);
    %errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2best(i,:),std_Signal_dl2best(i,:),'color',rgb('blue'));
    errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2del(i,:)-170,std_Signal_dl2del(i,:),'--','color',rgb('cornflowerblue'),'LineWidth',2);
    errorbar([0:bins_bestn1-1]/bins_bestn1,ave_combo_del(i,:),std_combo_del(i,:),'color',rgb('purple'),'LineWidth',2);
    %plot([0:bins_bestn1-1]/bins_bestn1,ave_combo_del(i,:),'color',rgb('purple'),'LineWidth',2);
    %errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldgood(i,:),std_Signal_zldgood(i,:),'color',rgb('green'));
    %errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldweak(i,:),std_Signal_zldweak(i,:),'color',rgb('brown'));
    %plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_control(i,:),'color',rgb('black'),'LineWidth',2); hold on
    %plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_del(i,:),'--','color',rgb('lightcoral'),'LineWidth',2)
    %plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_best(i,:),'color',rgb('red'),'LineWidth',2);
    %plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2best(i,:),'color',rgb('blue'),'LineWidth',2);
    %plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2del(i,:),'--','color',rgb('cornflowerblue'),'LineWidth',2);
    %plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldgood(i,:),'color',rgb('green'),'LineWidth',2);
    %plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldweak(i,:),'color',rgb('brown'),'LineWidth',2); hold off
    title('Signal Across Embryo')
    xlabel('Location')
    ylabel('Nuclei Signal')
    legend('control','Dl1 weak','Dl2 weak','Sum','Location', 'northwest')% 'dl2 strong', 'dl2 weak','Location', 'northwest')
    ylim([0 3000])
    xlim([0.1875 .6875])
    a = annotation('textbox',[.14 .10 .4 .5],'String',['Time = ',num2str(norm_time(i)*45)],'FitBoxToText','on');
    a.EdgeColor = 'none';
    F(i) = getframe(gcf);
    close(h)
end
video = VideoWriter('G:\My Drive\Research\Affinity Project\Short_True_Synergy_Signal.mp4');
video.FrameRate = 10;
video.Quality = 100;
open(video)
writeVideo(video,F);
close(video)
%}

%half dorsal

% for i = 1:length(norm_time)
%     h = figure;
%     errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_mcphet(i,:),std_Signal_mcphet(i,:),'color',rgb('blue')); hold on
%     errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dlhetcontrol(i,:),std_Signal_dlhetcontrol(i,:),'color',rgb('black')); 
%     errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dlhetdl2best(i,:),std_Signal_dlhetdl2best(i,:),'color',rgb('red'));
%     plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_mcphet(i,:),'color',rgb('blue'),'LineWidth',2); hold on
%     plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dlhetcontrol(i,:),'color',rgb('black'),'LineWidth',2)
%     plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dlhetdl2best(i,:),'color',rgb('red'),'LineWidth',2); hold off
%     title('Signal Across Embryo')
%     xlabel('Location')
%     ylabel('Nuclei Signal')
%     legend('mcp control', 'dlhet control','dlhet dl2 best','Location', 'northwest')
%     ylim([0 2000])
%     xlim([0 .9375])
%     a = annotation('textbox',[.14 .15 .4 .5],'String',['% NC14 = ',num2str(norm_time(i)*100)],'FitBoxToText','on');
%     a.EdgeColor = 'none';
%     F(i) = getframe(gcf);
%     close(h)
% end
% video = VideoWriter('G:\My Drive\Research\Affinity Project\test.mp4');
% video.FrameRate = 10;
% video.Quality = 100;
% open(video)
% writeVideo(video,F);
% close(video)

i = 181;
figure;
plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_control(i,:),'color',rgb('blue')); hold on
plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_best(i,:),'color',rgb('black'));
plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_del(i,:),'color',rgb('red'))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2best(i,:),'color',rgb('darkorange'));
plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2del(i,:),'color',rgb('purple'));
plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldgood(i,:),'color',rgb('green')); 
plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldweak(i,:),'color',rgb('brown')); 
legend('control', 'dl1 strong','dl1 weak', 'dl2 strong', 'dl2 weak','zld good','zld weak','Location', 'northwest')

i = 181;
figure;
plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_mcphet(i,:),'color',rgb('blue')); hold on
plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dlhetcontrol(i,:),'color',rgb('black'));
plot([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dlhetdl2best(i,:),'color',rgb('red'))
legend('mcp control', 'dl het control','Location', 'northwest')

%% Signal Replicate Plot

percent = .5;
sig = find(norm_time == percent);
c_point(1,:) = find(round(time_controln1,2) == percent);    %find closest point where std. dev. matches normalized time
c_point(2,:) = find(round(time_controln2,2) == percent);
c_point(3,:) = find(round(time_controln4,2) == percent);
c_point(4,:) = find(round(time_controln3,2) == percent);
c_point(5,:) = find(round(time_controln5,2) == percent);
c_point(6,:) = find(round(time_controln6,2) == percent);

b_point(1,:) = find(round(time_bestn1,2) == percent);
b_point(2,:) = find(round(time_bestn2,2) == percent);
b_point(3,:) = find(round(time_bestn4,2) == percent);

d_point(1,:) = find(round(time_deln1,2) == percent);
d_point(2,:) = find(round(time_deln2,2) == percent);
d_point(3,:) = find(round(time_deln4,2) == percent);

dl2b_point(1,:) = find(round(time_dl2bestn1,2) == percent);
dl2b_point(2,:) = find(round(time_dl2bestn2,2) == percent);
dl2b_point(3,:) = find(round(time_dl2bestn5,2) == percent);

dl2d_point(1,:) = find(round(time_dl2deln1,2) == percent);
dl2d_point(2,:) = find(round(time_dl2deln2,2) == percent);
dl2d_point(3,:) = find(round(time_dl2deln4,2) == percent);

zg_point(1,:) = find(round(time_zldgoodn2,2) == percent);
zg_point(2,:) = find(round(time_zldgoodn3,2) == percent);
zg_point(3,:) = find(round(time_zldgoodn4,2) == percent);

zw_point(1,:) = find(round(time_zldweakn1,2) == percent);
zw_point(2,:) = find(round(time_zldweakn2,2) == percent);
zw_point(3,:) = find(round(time_zldweakn5,2) == percent);


cR_point(1,1:2) = find(round(time_Rcontroln1,2) == percent);    %find closest point where std. dev. matches normalized time
cR_point(2,:) = find(round(time_Rcontroln2,2) == percent);
cR_point(3,:) = find(round(time_Rcontroln3,2) == percent);

%Error bars are standard error that is closest in time to interpolated
%signal

figure, errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_Rcontroln1(sig,:),stdSignal_bin_Rcn1(cR_point(1,end),:),'color',rgb('blue')); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_Rcontroln2(sig,:),stdSignal_bin_Rcn2(cR_point(2,end),:),'color',rgb('red'));
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_Rcontroln3(sig,:),stdSignal_bin_Rcn3(cR_point(3,end),:),'color',rgb('black')); 
title('Signal Across Embryo Control')
xlabel('Location')
ylabel('Nuclei Signal')
legend('n1','n2','n4','Location', 'northwest')
xlim([0 .9375])
ylim([0 2500])

figure, errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_bestn1(sig,:),stdSignal_bin_bestn1(b_point(1,end),:),'color',rgb('blue')); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_bestn2(sig,:),stdSignal_bin_bestn2(b_point(2,end),:),'color',rgb('red'));
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_bestn4(sig,:),stdSignal_bin_bestn4(b_point(3,end),:),'color',rgb('black')); 
title('Signal Across Embryo Strong')
xlabel('Location')
ylabel('Nuclei Signal')
legend('n1','n2','n4','Location', 'northwest')
xlim([0 .9375])
ylim([0 2500])

figure, errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_deln1(sig,:),stdSignal_bin_deln1(d_point(1,end),:),'color',rgb('blue')); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_deln2(sig,:),stdSignal_bin_deln2(d_point(2,end),:),'color',rgb('red'));
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_deln4(sig,:),stdSignal_bin_deln4(d_point(3,end),:),'color',rgb('black')); 
title('Signal Across Embryo Weak')
xlabel('Location')
ylabel('Nuclei Signal')
legend('n1','n2','n4','Location', 'northwest')
xlim([0 .9375])
ylim([0 2500])

figure, errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2bestn1(sig,:),stdSignal_bin_dl2bestn1(dl2b_point(1,end),:),'color',rgb('blue')); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2bestn2(sig,:),stdSignal_bin_dl2bestn2(dl2b_point(2,end),:),'color',rgb('red'));
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2bestn5(sig,:),stdSignal_bin_dl2bestn5(dl2b_point(3,end),:),'color',rgb('black')); 
title('Signal Across Embryo Dl2 Strong')
xlabel('Location')
ylabel('Nuclei Signal')
legend('n1','n2','n5','Location', 'northwest')
xlim([0 .9375])
ylim([0 2500])

figure, errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2deln1(sig,:),stdSignal_bin_dl2deln1(dl2d_point(1,end),:),'color',rgb('blue')); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2deln2(sig,:),stdSignal_bin_dl2deln2(dl2d_point(2,end),:),'color',rgb('red'));
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2deln4(sig,:),stdSignal_bin_dl2deln4(dl2d_point(3,end),:),'color',rgb('black')); 
title('Signal Across Embryo Dl2 Weak')
xlabel('Location')
ylabel('Nuclei Signal')
legend('n1','n2','n4','Location', 'northwest')
xlim([0 .9375])
ylim([0 2500])

figure, errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldgoodn2(sig,:),stdSignal_bin_zldgoodn2(zg_point(1,end),:),'color',rgb('blue')); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldgoodn3(sig,:),stdSignal_bin_zldgoodn3(zg_point(2,end),:),'color',rgb('red'));
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldgoodn4(sig,:),stdSignal_bin_zldgoodn4(zg_point(3,end),:),'color',rgb('black')); 
title('Signal Across Embryo Zld Good')
xlabel('Location')
ylabel('Nuclei Signal')
legend('n2','n3','n4','Location', 'northwest')
xlim([0 .9375])
ylim([0 2500])

figure, errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldweakn1(sig,:),stdSignal_bin_zldweakn1(zw_point(1,end),:),'color',rgb('blue')); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldweakn2(sig,:),stdSignal_bin_zldweakn2(zw_point(2,end),:),'color',rgb('red'));
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldweakn5(sig,:),stdSignal_bin_zldweakn5(zw_point(3,end),:),'color',rgb('black')); 
title('Signal Across Embryo Zld Weak')
xlabel('Location')
ylabel('Nuclei Signal')
legend('n1','n2','n5','Location', 'northwest')
xlim([0 .9375])
ylim([0 2500])


%% Output Nuclei
%{
figure, plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_controln1(end,:)); hold on
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_controln2(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_controln4(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_controln3(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_controln5(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_controln6(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetbestn1(end,:)); hold on
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetbestn2(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetbestn3(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_deln1(end,:)); hold on
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetdeln2(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetdeln4(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetdl2bestn1(end,:)); hold on
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetdl2bestn2(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetdl2bestn3(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetdl2deln1(end,:)); hold on
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetdl2deln2(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetdl2deln3(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetzldbestn1(end,:)); hold on
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetzldbestn2(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetzldbestn3(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_zldweakn1(end,:)); hold on
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_zldweakn2(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_zldweakn5(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_mcphetn1(point,:));hold on
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_mcphetn2(point,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_mcphetn3(point,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_dlhetcontroln1(point,:));hold on
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_dlhetcontroln3(point,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_dlhetcontroln4(point,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_dlhetdl2bestn2(point,:));hold on
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_dlhetdl2bestn3(point,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_dlhetdl2bestn5(point,:))
title('Signal Output Across Embryo')
xlabel('Location')
ylabel('Nuclei Signal')

%Revision
figure, plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_Rcontroln1(end,:)); hold on
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_Rcontroln2(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_Rcontroln3(end,:))
figure, plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_Rdzn1(end,:)); hold on
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_Rdzn2(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_Rdzn3(end,:))
figure, plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_Rdn1(end,:)); hold on
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_Rdn2(end,:))
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_Rdn3(end,:))
%}

%{
%SEM output test


ave_Area_control = nansum(cat(3,ave_Area_controln1.*lengthc(1,:), ave_Area_controln2.*lengthc(2,:), ave_Area_controln4.*lengthc(3,:)),3)./sum(lengthc);
ave_Area_best = nansum(cat(3,ave_Area_bestn1.*lengthb(1,:), ave_Area_bestn2.*lengthb(2,:), ave_Area_bestn4.*lengthb(3,:)),3)./sum(lengthb);
ave_Area_del = nansum(cat(3,ave_Area_deln1.*lengthd(1,:), ave_Area_deln2.*lengthd(2,:), ave_Area_deln4.*lengthd(3,:)),3)./sum(lengthd);
ave_Area_dl2best = nansum(cat(3,ave_Area_dl2bestn1.*lengthdl2b(1,:), ave_Area_dl2bestn2.*lengthdl2b(2,:), ave_Area_dl2bestn5.*lengthdl2b(3,:)),3)./sum(lengthdl2b);
ave_Area_dl2del = nansum(cat(3,ave_Area_dl2deln1.*lengthdl2d(1,:), ave_Area_dl2deln2.*lengthdl2d(2,:), ave_Area_dl2deln4.*lengthdl2d(3,:)),3)./sum(lengthdl2d);

ave_Area_zldgood = nansum(cat(3,ave_Area_zldgoodn2.*lengthzg(1,:), ave_Area_zldgoodn3.*lengthzg(2,:), ave_Area_zldgoodn4.*lengthzg(3,:)),3)./sum(lengthzg);
ave_Area_zldweak = (ave_Area_zldweakn1 + ave_Area_zldweakn2 + ave_Area_zldweakn5)./rep;


ave_Area_control = (ave_Area_controln1 + ave_Area_controln2 + ave_Area_controln4)./rep;
ave_Area_best = (ave_Area_bestn1 + ave_Area_bestn2 + ave_Area_bestn4)./rep;
ave_Area_del = (ave_Area_deln1 + ave_Area_deln2 + ave_Area_deln4)./rep;

ave_Area_dl2best = (ave_Area_dl2bestn1 + ave_Area_dl2bestn2 + ave_Area_dl2bestn5)./rep;
ave_Area_dl2del = (ave_Area_dl2deln1 + ave_Area_dl2deln2 + ave_Area_dl2deln4)./rep;

ave_Area_zldgood = (ave_Area_zldgoodn2 + ave_Area_zldgoodn3 + ave_Area_zldgoodn4)./rep;
ave_Area_zldweak = (ave_Area_zldweakn1 + ave_Area_zldweakn2 + ave_Area_zldweakn5)./rep;

% ave_Area_mcphet = (ave_Area_mcphetn1 + ave_Area_mcphetn2 + ave_Area_mcphetn3)./rep;
ave_Area_mcphet = (ave_Area_mcphetn2 + ave_Area_mcphetn3)./2;

ave_Area_dlhetcontrol = (ave_Area_dlhetcontroln1 + ave_Area_dlhetcontroln3 + ave_Area_dlhetcontroln4)./rep;
ave_Area_dlhetdl2best = (ave_Area_dlhetdl2bestn2 + ave_Area_dlhetdl2bestn3 + ave_Area_dlhetdl2bestn5)./rep;


%het
ave_Area_hetcontrol = (ave_Area_hetcontroln1 + ave_Area_hetcontroln2 + ave_Area_hetcontroln3)./rep;
ave_Area_hetbest = (ave_Area_hetbestn1 + ave_Area_hetbestn2 + ave_Area_hetbestn3)./rep;

ave_Area_hetdl2best = (ave_Area_hetdl2bestn1 + ave_Area_hetdl2bestn2 + ave_Area_hetdl2bestn3)./rep;
ave_Area_hetdl2del = (ave_Area_hetdl2deln1 + ave_Area_hetdl2deln2 + ave_Area_hetdl2deln3)./rep;

ave_Area_hetzldbest = (ave_Area_hetzldbestn1 + ave_Area_hetzldbestn2 + ave_Area_hetzldbestn3)./rep;

%calculate standard deviation
for i = 1:length(norm_time)
    for j = 1:bins_bestn1
        std_Area_control(i,j) = std([ave_Area_controln1(i,j) ave_Area_controln2(i,j) ave_Area_controln4(i,j)]);
        std_Area_best(i,j) = std([ave_Area_bestn1(i,j) ave_Area_bestn2(i,j) ave_Area_bestn4(i,j)]);
        std_Area_del(i,j) = std([ave_Area_deln1(i,j) ave_Area_deln2(i,j) ave_Area_deln4(i,j)]);
        std_Area_dl2best(i,j) = std([ave_Area_dl2bestn1(i,j) ave_Area_dl2bestn2(i,j) ave_Area_dl2bestn5(i,j)]);
        std_Area_dl2del(i,j) = std([ave_Area_dl2deln1(i,j) ave_Area_dl2deln2(i,j) ave_Area_dl2deln4(i,j)]);
        std_Area_zldgood(i,j) = std([ave_Area_zldgoodn2(i,j) ave_Area_zldgoodn3(i,j) ave_Area_zldgoodn4(i,j)]);
        std_Area_zldweak(i,j) = std([ave_Area_zldweakn1(i,j) ave_Area_zldweakn2(i,j) ave_Area_zldweakn5(i,j)]);
%         std_Area_mcphet(i,j) = std([ave_Area_mcphetn1(i,j) ave_Area_mcphetn2(i,j) ave_Area_mcphetn3(i,j)]);
        std_Area_mcphet(i,j) = std([ave_Area_mcphetn2(i,j) ave_Area_mcphetn3(i,j)]);
        std_Area_dlhetcontrol(i,j) = std([ave_Area_dlhetcontroln1(i,j) ave_Area_dlhetcontroln3(i,j) ave_Area_dlhetcontroln4(i,j)]);
        std_Area_dlhetdl2best(i,j) = std([ave_Area_dlhetdl2bestn2(i,j) ave_Area_dlhetdl2bestn3(i,j) ave_Area_dlhetdl2bestn5(i,j)]);
    end
end

%het
for i = 1:length(norm_time)
    for j = 1:bins_hetbestn1
        std_Area_hetcontrol(i,j) = std([ave_Area_hetcontroln1(i,j) ave_Area_hetcontroln2(i,j) ave_Area_hetcontroln3(i,j)]);
        std_Area_hetbest(i,j) = std([ave_Area_hetbestn1(i,j) ave_Area_hetbestn2(i,j) ave_Area_hetbestn3(i,j)]);
        std_Area_hetdl2best(i,j) = std([ave_Area_hetdl2bestn1(i,j) ave_Area_hetdl2bestn2(i,j) ave_Area_hetdl2bestn3(i,j)]);
        std_Area_hetdl2del(i,j) = std([ave_Area_hetdl2deln1(i,j) ave_Area_hetdl2deln2(i,j) ave_Area_hetdl2deln3(i,j)]);
        std_Area_hetzldbest(i,j) = std([ave_Area_hetzldbestn1(i,j) ave_Area_hetzldbestn2(i,j) ave_Area_hetzldbestn3(i,j)]);

    end
end
%}

% Have 0 become NaN
ave_Area_controln1(ave_Area_controln1 == 0) = nan;
ave_Area_controln2(ave_Area_controln2 == 0) = nan;
ave_Area_controln4(ave_Area_controln4 == 0) = nan;
ave_Area_controln3(ave_Area_controln3 == 0) = nan;
ave_Area_controln5(ave_Area_controln5 == 0) = nan;
ave_Area_controln6(ave_Area_controln6 == 0) = nan;
ave_Area_bestn1(ave_Area_bestn1 == 0) = nan;
ave_Area_bestn2(ave_Area_bestn2 == 0) = nan;
ave_Area_bestn4(ave_Area_bestn4 == 0) = nan;
ave_Area_deln1(ave_Area_bestn1 == 0) = nan;
ave_Area_deln2(ave_Area_deln2 == 0) = nan;
ave_Area_deln4(ave_Area_deln4 == 0) = nan;
ave_Area_dl2bestn1(ave_Area_dl2bestn1 == 0) = nan;
ave_Area_dl2bestn2(ave_Area_dl2bestn2 == 0) = nan;
ave_Area_dl2bestn5(ave_Area_dl2bestn5 == 0) = nan;
ave_Area_dl2deln1(ave_Area_dl2deln1 == 0) = nan;
ave_Area_dl2deln2(ave_Area_dl2deln2 == 0) = nan;
ave_Area_dl2deln4(ave_Area_dl2deln4 == 0) = nan;
ave_Area_zldgoodn2(ave_Area_zldgoodn2 == 0) = nan;
ave_Area_zldgoodn3(ave_Area_zldgoodn3 == 0) = nan;
ave_Area_zldgoodn4(ave_Area_zldgoodn4 == 0) = nan;
ave_Area_twideln1(ave_Area_twideln1 == 0) = nan;

%het
ave_Area_hetcontroln1(ave_Area_hetcontroln1 == 0) = nan;
ave_Area_hetcontroln2(ave_Area_hetcontroln2 == 0) = nan;
ave_Area_hetcontroln3(ave_Area_hetcontroln3 == 0) = nan;
ave_Area_hetbestn1(ave_Area_hetbestn1 == 0) = nan;
ave_Area_hetbestn2(ave_Area_hetbestn2 == 0) = nan;
ave_Area_hetbestn3(ave_Area_hetbestn3 == 0) = nan;
ave_Area_hetdl2bestn1(ave_Area_hetdl2bestn1 == 0) = nan;
ave_Area_hetdl2bestn2(ave_Area_hetdl2bestn2 == 0) = nan;
ave_Area_hetdl2bestn3(ave_Area_hetdl2bestn3 == 0) = nan;
ave_Area_hetdl2deln1(ave_Area_hetdl2deln1 == 0) = nan;
ave_Area_hetdl2deln2(ave_Area_hetdl2deln2 == 0) = nan;
ave_Area_hetdl2deln3(ave_Area_hetdl2deln3 == 0) = nan;
ave_Area_hetzldbestn1(ave_Area_hetzldbestn1 == 0) = nan;
ave_Area_hetzldbestn2(ave_Area_hetzldbestn2 == 0) = nan;
ave_Area_hetzldbestn3(ave_Area_hetzldbestn3 == 0) = nan;

%Revision
ave_Area_Rcontroln1(ave_Area_Rcontroln1 == 0) = nan;
ave_Area_Rcontroln2(ave_Area_Rcontroln2 == 0) = nan;
ave_Area_Rcontroln3(ave_Area_Rcontroln3 == 0) = nan;
ave_Area_Rdzn1(ave_Area_Rdzn1 == 0) = nan;
ave_Area_Rdzn2(ave_Area_Rdzn2 == 0) = nan;
ave_Area_Rdzn3(ave_Area_Rdzn3 == 0) = nan;
ave_Area_Rdn1(ave_Area_Rdn1 == 0) = nan;
ave_Area_Rdn2(ave_Area_Rdn2 == 0) = nan;
ave_Area_Rdn3(ave_Area_Rdn3 == 0) = nan;

% Non SEM
% ave_Area_control = nanmean(cat(3,ave_Area_controln1(activeM_controln1), ave_Area_controln2(activeM_controln2), ave_Area_controln4(activeM_controln4)),3);
% ave_Area_best = nanmean(cat(3,ave_Area_bestn1, ave_Area_bestn2, ave_Area_bestn4),3);
% ave_Area_del = nanmean(cat(3,ave_Area_deln1, ave_Area_deln2, ave_Area_deln4),3);
% ave_Area_dl2best = nanmean(cat(3,ave_Area_dl2bestn1, ave_Area_dl2bestn2, ave_Area_dl2bestn5),3);
% ave_Area_dl2del = nanmean(cat(3,ave_Area_dl2deln1, ave_Area_dl2deln2, ave_Area_dl2deln4),3);
% ave_Area_zldgood = nanmean(cat(3,ave_Area_zldgoodn2, ave_Area_zldgoodn3, ave_Area_zldgoodn4),3);
% ave_Area_zldweak = nanmean(cat(3,ave_Area_zldweakn1, ave_Area_zldweakn2, ave_Area_zldweakn5),3);
% 
% std_Area_control = nanstd(cat(3,ave_Area_controln1, ave_Area_controln2, ave_Area_controln4),[],3);
% std_Area_best = nanstd(cat(3,ave_Area_bestn1, ave_Area_bestn2, ave_Area_bestn4),[],3);
% std_Area_del = nanstd(cat(3,ave_Area_deln1, ave_Area_deln2, ave_Area_deln4),[],3);
% std_Area_dl2best = nanstd(cat(3,ave_Area_dl2bestn1, ave_Area_dl2bestn2, ave_Area_dl2bestn5),[],3);
% std_Area_dl2del = nanstd(cat(3,ave_Area_dl2deln1, ave_Area_dl2deln2, ave_Area_dl2deln4),[],3);
% std_Area_zldgood = nanstd(cat(3,ave_Area_zldgoodn2, ave_Area_zldgoodn3, ave_Area_zldgoodn4),[],3);
% std_Area_zldweak = nanstd(cat(3,ave_Area_zldweakn1, ave_Area_zldweakn2, ave_Area_zldweakn5),[],3);

%SEM



%Active Nuclei Only
% for i = 1:16
%     %For getting active nuclei
%     [~,ic1{i}] = intersect(shift_bin_controln1{i},activeM_controln1);
%     [~,ic2{i}] = intersect(shift_bin_controln2{i},activeM_controln2);
%     [~,ic4{i}] = intersect(shift_bin_controln4{i},activeM_controln4);
%     [~,ib1{i}] = intersect(shift_bin_bestn1{i},activeM_bestn1);
%     [~,ib2{i}] = intersect(shift_bin_bestn2{i},activeM_bestn2);
%     [~,ib4{i}] = intersect(shift_bin_bestn4{i},activeM_bestn4);
%     [~,id1{i}] = intersect(shift_bin_deln1{i},activeM_deln1);
%     [~,id2{i}] = intersect(shift_bin_deln2{i},activeM_deln2);
%     [~,id4{i}] = intersect(shift_bin_deln4{i},activeM_deln4);
%     [~,idl2b1{i}] = intersect(shift_bin_dl2bestn1{i},activeM_dl2bestn1);
%     [~,idl2b2{i}] = intersect(shift_bin_dl2bestn2{i},activeM_dl2bestn2);
%     [~,idl2b5{i}] = intersect(shift_bin_dl2bestn5{i},activeM_dl2bestn5);
%     [~,idl2d1{i}] = intersect(shift_bin_dl2deln1{i},activeM_dl2deln1);
%     [~,idl2d2{i}] = intersect(shift_bin_dl2deln2{i},activeM_dl2deln2);
%     [~,idl2d4{i}] = intersect(shift_bin_dl2deln4{i},activeM_dl2deln4);
%     [~,izg2{i}] = intersect(shift_bin_zldgoodn2{i},activeM_zldgoodn2);
%     [~,izg3{i}] = intersect(shift_bin_zldgoodn3{i},activeM_zldgoodn3);
%     [~,izg4{i}] = intersect(shift_bin_zldgoodn4{i},activeM_zldgoodn4);
%     
%     Area_c{i} = [Area_bin_cn1{end,i}(ic1{i}) Area_bin_cn2{end,i}(ic2{i}) Area_bin_cn4{end,i}(ic4{i})];
%     Area_b{i} = [Area_bin_bn1{end,i}(ib1{i}) Area_bin_bn2{end,i}(ib2{i}) Area_bin_bn4{end,i}(ib4{i})];
%     Area_d{i} = [Area_bin_dn1{end,i}(id1{i}) Area_bin_dn2{end,i}(id2{i}) Area_bin_dn4{end,i}(id4{i})];
%     Area_dl2b{i} = [Area_bin_dl2bn1{end,i}(idl2b1{i}) Area_bin_dl2bn2{end,i}(idl2b2{i}) Area_bin_dl2bn5{end,i}(idl2b5{i})];
%     Area_dl2d{i} = [Area_bin_dl2dn1{end,i}(idl2d1{i}) Area_bin_dl2dn2{end,i}(idl2d2{i}) Area_bin_dl2dn4{end,i}(idl2d4{i})];
%     Area_zg{i} = [Area_bin_zgn2{end,i}(izg2{i}) Area_bin_zgn3{end,i}(izg3{i}) Area_bin_zgn4{end,i}(izg4{i})];
%     
%     if isempty(Area_c{i})
%         Area_c{i} = [];
%     end
%     if isempty(Area_b{i})
%         Area_b{i} = [];
%     end
%     if isempty(Area_d{i})
%         Area_d{i} = [];
%     end
%     if isempty(Area_dl2b{i})
%         Area_dl2b{i} = [];
%     end
%     if isempty(Area_dl2d{i})
%         Area_dl2d{i} = [];
%     end
%     if isempty(Area_zg{i})
%         Area_zg{i} = [];
%     end
% end



%Combine Output into one matrix

%All nuclei 
for i = 1:16
    Area_c{i} = [Area_bin_cn1{end,i} Area_bin_cn2{end,i} Area_bin_cn4{end,i} Area_bin_cn3{end,i} Area_bin_cn5{end,i} Area_bin_cn6{end,i}];
    Area_b{i} = [Area_bin_bn1{end,i} Area_bin_bn2{end,i} Area_bin_bn4{end,i}];
    Area_d{i} = [Area_bin_dn1{end,i} Area_bin_dn2{end,i} Area_bin_dn4{end,i}];
    Area_dl2b{i} = [Area_bin_dl2bn1{end,i} Area_bin_dl2bn2{end,i} Area_bin_dl2bn5{end,i}];
    Area_dl2d{i} = [Area_bin_dl2dn1{end,i} Area_bin_dl2dn2{end,i} Area_bin_dl2dn4{end,i}];
    Area_zg{i} = [Area_bin_zgn2{end,i} Area_bin_zgn3{end,i} Area_bin_zgn4{end,i}];  
    Area_td{i} = [Area_bin_tdn1{end,i} Area_bin_tdn2{end,i}];

    %Revision
    Area_Rc{i} = [Area_bin_Rcn1{end,i} Area_bin_Rcn2{end,i} Area_bin_Rcn3{end,i}];
    Area_Rdz{i} = [Area_bin_Rdzn1{end,i} Area_bin_Rdzn2{end,i} Area_bin_Rdzn3{end,i}];
    Area_Rd{i} = [Area_bin_Rdn1{end,i} Area_bin_Rdn2{end,i} Area_bin_Rdn3{end,i}];
end

for i = 1:16
    %Mean
    ave_Area_control(i) = mean(Area_c{i});
    ave_Area_best(i) = mean(Area_b{i});
    ave_Area_del(i) = mean(Area_d{i});
    ave_Area_dl2best(i) = mean(Area_dl2b{i});
    ave_Area_dl2del(i) = mean(Area_dl2d{i});
    ave_Area_zldgood(i) = mean(Area_zg{i});
    ave_Area_twidel(i) = mean(Area_td{i});

    %Revision
    ave_Area_Rcontrol(i) = mean(Area_Rc{i});
    ave_Area_Rdz(i) = mean(Area_Rdz{i});
    ave_Area_Rd(i) = mean(Area_Rd{i});

    
    %Standard Error of the Mean for Output
    std_Area_control(i) = std(Area_c{i})./sqrt(length(Area_c{i}));
    std_Area_best(i) = std(Area_b{i})./sqrt(length(Area_b{i}));
    std_Area_del(i) = std(Area_d{i})./sqrt(length(Area_d{i}));
    std_Area_dl2best(i) = std(Area_dl2b{i})./sqrt(length(Area_dl2b{i}));
    std_Area_dl2del(i) = std(Area_dl2d{i})./sqrt(length(Area_dl2d{i}));
    std_Area_zldgood(i) = std(Area_zg{i})./sqrt(length(Area_zg{i}));
    std_Area_twidel(i) = std(Area_td{i})./sqrt(length(Area_td{i}));

    %Revision
    std_Area_Rcontrol(i) = std(Area_Rc{i})./sqrt(length(Area_Rc{i}));
    std_Area_Rdz(i) = std(Area_Rdz{i})./sqrt(length(Area_Rdz{i}));
    std_Area_Rd(i) = std(Area_Rd{i})./sqrt(length(Area_Rd{i}));


end





%het
ave_Area_hetcontrol = nanmean(cat(3,ave_Area_hetcontroln1, ave_Area_hetcontroln2, ave_Area_hetcontroln3),3);
ave_Area_hetbest = nanmean(cat(3,ave_Area_hetbestn1, ave_Area_hetbestn2, ave_Area_hetbestn3),3);
ave_Area_hetdl2best = nanmean(cat(3,ave_Area_hetdl2bestn1, ave_Area_hetdl2bestn2, ave_Area_hetdl2bestn3),3);
ave_Area_hetdl2del = nanmean(cat(3,ave_Area_hetdl2deln1, ave_Area_hetdl2deln2, ave_Area_hetdl2deln3),3);
ave_Area_hetzldbest = nanmean(cat(3,ave_Area_hetzldbestn1, ave_Area_hetzldbestn2, ave_Area_hetzldbestn3),3);

std_Area_hetcontrol = nanstd(cat(3,ave_Area_hetcontroln1, ave_Area_hetcontroln2, ave_Area_hetcontroln3),[],3);
std_Area_hetbest = nanstd(cat(3,ave_Area_hetbestn1, ave_Area_hetbestn2, ave_Area_hetbestn3),[],3);
std_Area_hetdl2best = nanstd(cat(3,ave_Area_hetdl2bestn1, ave_Area_hetdl2bestn2, ave_Area_hetdl2bestn3),[],3);
std_Area_hetdl2del = nanstd(cat(3,ave_Area_hetdl2deln1, ave_Area_hetdl2deln2, ave_Area_hetdl2deln3),[],3);
std_Area_hetzldbest = nanstd(cat(3,ave_Area_hetzldbestn1, ave_Area_hetzldbestn2, ave_Area_hetzldbestn3),[],3);


ave_Area_combo = nansum([ave_Area_del(end,:)-5.97e3; ave_Area_dl2del(end,:)-5.97e3]);
std_Area_combo = nansum([std_Area_del; std_Area_dl2del]);
%{
%Poster Plots
%go output
figure;
% subplot(1,2,1)
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_control(point,:),std_Area_control(:),'color',rgb('black'),'LineWidth',2); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_best(point,:),std_Area_best(:),'color',rgb('red'),'LineWidth',2);  hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_del(point,:),std_Area_del(:),'--','color',rgb('lightcoral'),'LineWidth',2); hold on
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_control(point,:),'color',rgb('black')); hold on
% plot([0:bins_bestn1-1]/bins_bestn1, ave_Area_best(point,:),'color',rgb('red'))
% plot([0:bins_bestn1-1]/bins_bestn1, ave_Area_del(point,:),'--','color',rgb('red'))
title('Average Nuclei Output Across Embryo')
xlabel('Location')
ylabel('Nuclei Output (A.U.)')
legend('control', 'Dl1 strong','Dl1 weak','Location', 'northwest')
xlim([0 .9375])
ylim([0 6e4])

figure,
% subplot(1,2,2)
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_control(point,:),std_Area_control(point,:),'color',rgb('black'),'LineWidth',2); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_dl2best(point,:),std_Area_dl2best(point,:),'color',rgb('blue'),'LineWidth',2); 
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_dl2del(point,:),std_Area_dl2del(point,:),'--','color',rgb('cornflowerblue'),'LineWidth',2); %-5.97e3
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_control(point,:),'color',rgb('black')); hold on
% plot([0:bins_bestn1-1]/bins_bestn1, ave_Area_dl2best(point,:),'color',rgb('blue'))
% plot([0:bins_bestn1-1]/bins_bestn1, ave_Area_dl2del(point,:),'--','color',rgb('blue'))
% shadedErrorBar([.3 .4 .5],[-500 -500 -500],[1 1 1],{'LineWidth',2},1)
title('Average Nuclei Output Across Embryo')
xlabel('Location')
ylabel('Nuclei Output (A.U.)')
legend('control', 'Dl2 strong', 'Dl2 weak','Location', 'northwest')
xlim([0 .9375])
ylim([0 6e4])

subplot(1,2,2)
figure,
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_control(point,:),std_Area_control(point,:),'color',rgb('black')); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_zldgood(point,:),std_Area_zldgood(point,:),'color',rgb('green'),'LineWidth',2);
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_combo,std_Area_combo,'color', rgb('purple'),'LineWidth',2);
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_control(point,:),'color',rgb('black')); hold on
% plot([0:bins_bestn1-1]/bins_bestn1, ave_Area_zldgood(point,:),'color',rgb('green'))
title('Average Nuclei Output Across Embryo')
xlabel('Location')
ylabel('Nuclei Output')
legend('control','Dl1 Strong','Dl2 Strong','Zld strong','Location', 'northwest')
xlim([0 .9375])
ylim([0 6e4])

legend('control','Dl1 weak','Dl2 weak','Combo','Location', 'northwest')

%Revision
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_control(point,:),std_Area_control(:),'color',rgb('gray'),'LineWidth',2); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_Rcontrol(point,:),std_Area_Rcontrol(:),'color',rgb('black'),'LineWidth',2); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_Rdz(point,:),std_Area_Rdz(:),'color',rgb('red'),'LineWidth',2); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_Rd(point,:),std_Area_Rd(:),'color',rgb('blue'),'LineWidth',2); hold on

title('Average Nuclei Output Across Embryo')
xlabel('Location')
ylabel('Nuclei Output (A.U.)')
legend('control', 'Revision Control','Location', 'northwest')
xlim([0 .9375])
ylim([0 6e4])
%}


point = 1;
figure; 
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_control(point,:),std_Area_control(point,:),'color',rgb('blue')); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_del(point,:),std_Area_del(point,:),'color',rgb('red'));
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_best(point,:),std_Area_best(point,:),'color',rgb('black')); 
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_dl2best(point,:),std_Area_dl2best(point,:),'color',rgb('darkorange')); 
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_dl2del(point,:),std_Area_dl2del(point,:),'color',rgb('purple')); 
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_zldgood(point,:),std_Area_zldgood(point,:),'color',rgb('green'));
% errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_zldweak(point,:),std_Area_zldweak(point,:),'color',rgb('brown'));
plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_control(point,:),'color',rgb('blue')); hold on
plot([0:bins_bestn1-1]/bins_bestn1, ave_Area_del(point,:),'color',rgb('red'))
plot([0:bins_bestn1-1]/bins_bestn1, ave_Area_best(point,:),'color',rgb('black'))
plot([0:bins_bestn1-1]/bins_bestn1, ave_Area_dl2best(point,:),'color',rgb('darkorange'))
plot([0:bins_bestn1-1]/bins_bestn1, ave_Area_dl2del(point,:),'color',rgb('purple'))
plot([0:bins_bestn1-1]/bins_bestn1, ave_Area_zldgood(point,:),'color',rgb('green'))
% plot([0:bins_bestn1-1]/bins_bestn1, ave_Area_zldweak(point,:),'color',rgb('brown'))
title(['Average Nuclei Output Across Embryo at ' num2str(norm_time(point)*100) '% through NC14'])
xlabel('Location')
ylabel('Nuclei Output')
legend('control', 'dl1 strong','dl1 weak', 'dl2 strong', 'dl2 weak','zld good','zld weak','Location', 'northwest')
xlim([0 .9375])
ylim([0 6.8e4])

%half dorsal
% point = 181;
% figure; 
% errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_mcphet(point,:),std_Area_mcphet(point,:),'color',rgb('blue')); hold on
% errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_dlhetcontrol(point,:),std_Area_dlhetcontrol(point,:),'color',rgb('red'));
% errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_dlhetdl2best(point,:),std_Area_dlhetdl2best(point,:),'color',rgb('black')); 
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_mcphet(point,:),'color',rgb('blue')); hold on
% plot([0:bins_bestn1-1]/bins_bestn1, ave_Area_dlhetcontrol(point,:),'color',rgb('red'))
% plot([0:bins_bestn1-1]/bins_bestn1, ave_Area_dlhetdl2best(point,:),'color',rgb('black'))
% title(['Average Nuclei Output Across Embryo at ' num2str(norm_time(point)*100) '% through NC14'])
% xlabel('Location')
% ylabel('Nuclei Output')
% legend('mcp control', 'dl het control','dl het dl2 best','Location', 'northwest')
% xlim([0 .9375])
% ylim([0 3.5e4])

%het
figure;
errorbar([0:bins_hetbestn1-1]/bins_hetbestn1,ave_Area_hetcontrol(point,:),std_Area_hetcontrol(point,:),'color',rgb('black')); hold on
errorbar([0:bins_hetbestn1-1]/bins_hetbestn1,ave_Area_hetbest(point,:),std_Area_hetbest(point,:),'color',rgb('red'));  hold on
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_control(point,:),'color',rgb('black')); hold on
% plot([0:bins_bestn1-1]/bins_bestn1, ave_Area_best(point,:),'color',rgb('red'))
% plot([0:bins_bestn1-1]/bins_bestn1, ave_Area_del(point,:),'--','color',rgb('red'))
title('Average Nuclei Output Across Embryo')
xlabel('Location')
ylabel('Nuclei Output (A.U.)')
legend('het control', 'het Dl1 strong','Location', 'northwest')
xlim([0 .9375])
ylim([0 6.8e4])

figure,
errorbar([0:bins_hetbestn1-1]/bins_hetbestn1,ave_Area_hetcontrol(point,:),std_Area_hetcontrol(point,:),'color',rgb('black')); hold on
errorbar([0:bins_hetbestn1-1]/bins_hetbestn1,ave_Area_hetdl2best(point,:),std_Area_hetdl2best(point,:),'color',rgb('blue')); 
errorbar([0:bins_hetbestn1-1]/bins_hetbestn1,ave_Area_hetdl2del(point,:),std_Area_hetdl2del(point,:),'--','color',rgb('cornflowerblue')); 
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_control(point,:),'color',rgb('black')); hold on
% plot([0:bins_bestn1-1]/bins_bestn1, ave_Area_dl2best(point,:),'color',rgb('blue'))
% plot([0:bins_bestn1-1]/bins_bestn1, ave_Area_dl2del(point,:),'--','color',rgb('blue'))
shadedErrorBar([.3 .4 .5],[-500 -500 -500],[1 1 1],{'LineWidth',2},1)
title('Average Nuclei Output Across Embryo')
xlabel('Location')
ylabel('Nuclei Output (A.U.)')
legend('het control', 'het Dl2 strong', 'het Dl2 weak','Location', 'northwest')
xlim([0 .9375])
ylim([0 6.8e4])

figure,
errorbar([0:bins_hetbestn1-1]/bins_hetbestn1,ave_Area_hetcontrol(point,:),std_Area_hetcontrol(point,:),'color',rgb('black')); hold on
errorbar([0:bins_hetbestn1-1]/bins_hetbestn1,ave_Area_hetzldbest(point,:),std_Area_hetzldbest(point,:),'color',rgb('green'));
% plot([0:bins_bestn1-1]/bins_bestn1,ave_Area_control(point,:),'color',rgb('black')); hold on
% plot([0:bins_bestn1-1]/bins_bestn1, ave_Area_hetzldbest(point,:),'color',rgb('green'))
title('Average Nuclei Output Across Embryo')
xlabel('Location')
ylabel('Nuclei Output')
legend('het control','het Dl1 Strong','het Dl2 Strong','het Zld strong','Location', 'northwest')
xlim([0 .9375])
ylim([0 6.8e4])


%% Output Replicate Plot

percent = 1;
sig = find(norm_time == percent);
c_point(1,:) = find(round(time_controln1,2) == percent);
c_point(2,:) = find(round(time_controln2,2) == percent);
c_point(3,:) = find(round(time_controln4,2) == percent);
c_point(4,:) = find(round(time_controln3,2) == percent);
c_point(5,:) = find(round(time_controln5,2) == percent);
c_point(6,:) = find(round(time_controln6,2) == percent);

b_point(1,:) = find(round(time_bestn1,2) == percent);
b_point(2,:) = find(round(time_bestn2,2) == percent);
b_point(3,:) = find(round(time_bestn4,2) == percent);

d_point(1,:) = find(round(time_deln1,2) == percent);
d_point(2,:) = find(round(time_deln2,2) == percent);
d_point(3,:) = find(round(time_deln4,2) == percent);

dl2b_point(1,:) = find(round(time_dl2bestn1,2) == percent);
dl2b_point(2,:) = find(round(time_dl2bestn2,2) == percent);
dl2b_point(3,:) = find(round(time_dl2bestn5,2) == percent);

dl2d_point(1,:) = find(round(time_dl2deln1,2) == percent);
dl2d_point(2,:) = find(round(time_dl2deln2,2) == percent);
dl2d_point(3,:) = find(round(time_dl2deln4,2) == percent);

zg_point(1,:) = find(round(time_zldgoodn2,2) == percent);
zg_point(2,:) = find(round(time_zldgoodn3,2) == percent);
zg_point(3,:) = find(round(time_zldgoodn4,2) == percent);

zw_point(1,:) = find(round(time_zldweakn1,2) == percent);
zw_point(2,:) = find(round(time_zldweakn2,2) == percent);
zw_point(3,:) = find(round(time_zldweakn5,2) == percent);

figure, errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_controln1(sig,:),stdArea_bin_controln1(c_point(1,end)-1,:),'color',rgb('blue')); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_controln2(sig,:),stdArea_bin_controln2(c_point(2,end)-1,:),'color',rgb('red'));
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_controln3(sig,:),stdArea_bin_controln3(c_point(3,end)-1,:),'color',rgb('black')); 
title('Output Across Embryo Het Control')
xlabel('Location')
ylabel('Nuclei Signal')
legend('n1','n2','n3','Location', 'northwest')
xlim([0 .9375])
ylim([0 8e4])

figure, errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetbestn1(sig,:),stdArea_bin_hetbestn1(b_point(1,end)-1,:),'color',rgb('blue')); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetbestn2(sig,:),stdArea_bin_hetbestn2(b_point(2,end)-1,:),'color',rgb('red'));
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetbestn3(sig,:),stdArea_bin_hetbestn3(b_point(3,end)-1,:),'color',rgb('black')); 
title('Output Across Embryo Het Dl1 Strong')
xlabel('Location')
ylabel('Nuclei Signal')
legend('n1','n2','n3','Location', 'northwest')
xlim([0 .9375])
ylim([0 8e4])

figure, errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_deln1(sig,:),stdArea_bin_deln1(d_point(1,end)-1,:),'color',rgb('blue')); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_deln2(sig,:),stdArea_bin_deln2(d_point(2,end)-1,:),'color',rgb('red'));
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_deln4(sig,:),stdArea_bin_deln4(d_point(3,end)-1,:),'color',rgb('black')); 
title('Output Across Embryo Weak')
xlabel('Location')
ylabel('Nuclei Signal')
legend('n1','n2','n4','Location', 'northwest')
xlim([0 .9375])
ylim([0 8e4])

figure, errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetdl2bestn1(sig,:),stdArea_bin_hetdl2bestn1(dl2b_point(1,end)-1,:),'color',rgb('blue')); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetdl2bestn2(sig,:),stdArea_bin_hetdl2bestn2(dl2b_point(2,end)-1,:),'color',rgb('red'));
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetdl2bestn3(sig,:),stdArea_bin_hetdl2bestn3(dl2b_point(3,end)-1,:),'color',rgb('black')); 
title('Output Across Embryo het Dl2 Strong')
xlabel('Location')
ylabel('Nuclei Signal')
legend('n1','n2','n3','Location', 'northwest')
xlim([0 .9375])
ylim([0 8e4])

figure, errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetdl2deln1(sig,:),stdArea_bin_hetdl2deln1(dl2d_point(1,end)-1,:),'color',rgb('blue')); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetdl2deln2(sig,:),stdArea_bin_hetdl2deln2(dl2d_point(2,end)-1,:),'color',rgb('red'));
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetdl2deln3(sig,:),stdArea_bin_hetdl2deln3(dl2d_point(3,end)-1,:),'color',rgb('black')); 
title('Output Across Embryo het Dl2 Weak')
xlabel('Location')
ylabel('Nuclei Signal')
legend('n1','n2','n3','Location', 'northwest')
xlim([0 .9375])
ylim([0 8e4])

figure, errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetzldbestn1(sig,:),stdArea_bin_hetzldbestn1(zg_point(1,end)-1,:),'color',rgb('blue')); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetzldbestn2(sig,:),stdArea_bin_hetzldbestn2(zg_point(2,end)-1,:),'color',rgb('red'));
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_hetzldbestn3(sig,:),stdArea_bin_hetzldbestn3(zg_point(3,end)-1,:),'color',rgb('black')); 
title('Output Across Embryo het Zld Good')
xlabel('Location')
ylabel('Nuclei Signal')
legend('n1','n2','n3','Location', 'northwest')
xlim([0 .9375])
ylim([0 8e4])

figure, errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_zldweakn1(sig,:),stdArea_bin_zldweakn1(zw_point(1,end)-1,:),'color',rgb('blue')); hold on
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_zldweakn2(sig,:),stdArea_bin_zldweakn2(zw_point(2,end)-1,:),'color',rgb('red'));
errorbar([0:bins_bestn1-1]/bins_bestn1,ave_Area_zldweakn5(sig,:),stdArea_bin_zldweakn5(zw_point(3,end)-1,:),'color',rgb('black')); 
title('Output Across Embryo Zld Good')
xlabel('Location')
ylabel('Nuclei Signal')
legend('n1','n2','n5','Location', 'northwest')
xlim([0 .9375])
ylim([0 8e4])


%% Active Domain
% Non SEM
% ave_Area_control = nanmean(cat(3,ave_Area_controln1, ave_Area_controln2, ave_Area_controln4, ave_Area_controln3, ave_Area_controln5, ave_Area_controln6),3);
% ave_Area_best = nanmean(cat(3,ave_Area_bestn1, ave_Area_bestn2, ave_Area_bestn4),3);
% ave_Area_del = nanmean(cat(3,ave_Area_deln1, ave_Area_deln2, ave_Area_deln4),3);
% ave_Area_dl2best = nanmean(cat(3,ave_Area_dl2bestn1, ave_Area_dl2bestn2, ave_Area_dl2bestn5),3);
% ave_Area_dl2del = nanmean(cat(3,ave_Area_dl2deln1, ave_Area_dl2deln2, ave_Area_dl2deln4),3);
% ave_Area_zldgood = nanmean(cat(3,ave_Area_zldgoodn2, ave_Area_zldgoodn3, ave_Area_zldgoodn4),3);
% ave_Area_zldweak = nanmean(cat(3,ave_Area_zldweakn1, ave_Area_zldweakn2, ave_Area_zldweakn5),3);

% std_Area_control = nanstd(cat(3,ave_Area_controln1, ave_Area_controln2, ave_Area_controln4),[],3);
% std_Area_best = nanstd(cat(3,ave_Area_bestn1, ave_Area_bestn2, ave_Area_bestn4),[],3);
% std_Area_del = nanstd(cat(3,ave_Area_deln1, ave_Area_deln2, ave_Area_deln4),[],3);
% std_Area_dl2best = nanstd(cat(3,ave_Area_dl2bestn1, ave_Area_dl2bestn2, ave_Area_dl2bestn5),[],3);
% std_Area_dl2del = nanstd(cat(3,ave_Area_dl2deln1, ave_Area_dl2deln2, ave_Area_dl2deln4),[],3);
% std_Area_zldgood = nanstd(cat(3,ave_Area_zldgoodn2, ave_Area_zldgoodn3, ave_Area_zldgoodn4),[],3);
% std_Area_zldweak = nanstd(cat(3,ave_Area_zldweakn1, ave_Area_zldweakn2, ave_Area_zldweakn5),[],3);

domain = 0:.005:.9375;
act_thresh = 260; %threshold adjusts for base signal (inactive nuclei)

act_threshR = 200;

for i = 1:length(norm_time)
    %Ave Signal Domain - interpolating average signal
    signal_domain_controln1(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_controln1(i,:),domain);
    signal_domain_controln2(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_controln2(i,:),domain);
    signal_domain_controln4(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_controln4(i,:),domain);
    signal_domain_controln3(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_controln3(i,:),domain);
    signal_domain_controln5(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_controln5(i,:),domain);
    signal_domain_controln6(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_controln6(i,:),domain);
    signal_domain_bestn1(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_bestn1(i,:),domain);
    signal_domain_bestn2(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_bestn2(i,:),domain);
    signal_domain_bestn4(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_bestn4(i,:),domain);
    signal_domain_deln1(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_deln1(i,:),domain);
    signal_domain_deln2(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_deln2(i,:),domain);
    signal_domain_deln4(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_deln4(i,:),domain);
    signal_domain_dl2bestn1(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2bestn1(i,:),domain);
    signal_domain_dl2bestn2(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2bestn2(i,:),domain);
    signal_domain_dl2bestn5(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2bestn5(i,:),domain);
    signal_domain_dl2deln1(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2deln1(i,:),domain);
    signal_domain_dl2deln2(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2deln2(i,:),domain);
    signal_domain_dl2deln4(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dl2deln4(i,:),domain);
    signal_domain_zldgoodn2(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldgoodn2(i,:),domain);
    signal_domain_zldgoodn3(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldgoodn3(i,:),domain);
    signal_domain_zldgoodn4(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldgoodn4(i,:),domain);
%     signal_domain_zldweakn1(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldweakn1(i,:),domain);
%     signal_domain_zldweakn2(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldweakn2(i,:),domain);
%     signal_domain_zldweakn5(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_zldweakn5(i,:),domain);
    signal_domain_twideln1(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_twideln1(i,:),domain);
    signal_domain_twideln2(i,:) = interp1([0:bins_bestn1-1]/bins_bestn2,ave_Signal_twideln2(i,:),domain);

    signal_domain_mcphet(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_mcphet(i,:),domain);
    signal_domain_dlhetcontrol(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dlhetcontrol(i,:),domain);
    signal_domain_dlhetdl2best(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_dlhetdl2best(i,:),domain);
    
    signal_domain_Rcontroln1(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_Rcontroln1(i,:),domain);
    signal_domain_Rcontroln2(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_Rcontroln2(i,:),domain);
    signal_domain_Rcontroln3(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_Rcontroln3(i,:),domain);
    
    signal_domain_Rdzn1(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_Rdzn1(i,:),domain);
    signal_domain_Rdzn2(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_Rdzn2(i,:),domain);
    signal_domain_Rdzn3(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_Rdzn3(i,:),domain);

    signal_domain_Rdn1(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_Rdn1(i,:),domain);
    signal_domain_Rdn2(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_Rdn2(i,:),domain);
    signal_domain_Rdn3(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Signal_Rdn3(i,:),domain);


    %Calculating domain size
    
    [height_cn1(i), tmp] = max(signal_domain_controln1(i,:)-act_thresh);    %get signal peak and subtract threshold
    height_cn1(i) = height_cn1(i)/2;    %find half height
    [~, idx1(i)] = min(abs(signal_domain_controln1(i,1:tmp)-act_thresh - height_cn1(i)));   %get the two points where signal is equal to half height
    [~, idx2(i)] = min(abs(signal_domain_controln1(i,tmp:end)-act_thresh - height_cn1(i)));
    expression_span_cn1(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));  %calculate distance between those two points
    if height_cn1(i) < 0
        expression_span_cn1(i) = 0;
    end
    if isnan(signal_domain_controln1(i,1)) == 1 
        expression_span_cn1(i) = 0;
    end
    
    [height_cn2(i), tmp] = max(signal_domain_controln2(i,:)-act_thresh);    %get signal peak and subtract threshold
    height_cn2(i) = height_cn2(i)/2;    %find half height
    [~, idx1(i)] = min(abs(signal_domain_controln2(i,1:tmp)-act_thresh - height_cn2(i)));   %get the two points where signal is equal to half height
    [~, idx2(i)] = min(abs(signal_domain_controln2(i,tmp:end)-act_thresh - height_cn2(i)));
    expression_span_cn2(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));  %calculate distance between those two points
    if height_cn2(i) < 0
        expression_span_cn2(i) = 0;
    end
    if isnan(signal_domain_controln2(i,1)) == 1 
        expression_span_cn2(i) = 0;
    end
    
    [height_cn4(i), tmp] = max(signal_domain_controln4(i,:)-act_thresh);    %get signal peak and subtract threshold
    height_cn4(i) = height_cn4(i)/2;    %find half height
    [~, idx1(i)] = min(abs(signal_domain_controln4(i,1:tmp)-act_thresh - height_cn4(i)));   %get the two points where signal is equal to half height
    [~, idx2(i)] = min(abs(signal_domain_controln4(i,tmp:end)-act_thresh - height_cn4(i)));
    expression_span_cn4(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));  %calculate distance between those two points
    if height_cn4(i) < 0
        expression_span_cn4(i) = 0;
    end
    if isnan(signal_domain_controln4(i,1)) == 1 
        expression_span_cn4(i) = 0;
    end
    
    [height_cn3(i), tmp] = max(signal_domain_controln3(i,:)-act_thresh);    %get signal peak and subtract threshold
    height_cn3(i) = height_cn3(i)/2;    %find half height
    [~, idx1(i)] = min(abs(signal_domain_controln3(i,1:tmp)-act_thresh - height_cn3(i)));   %get the two points where signal is equal to half height
    [~, idx2(i)] = min(abs(signal_domain_controln3(i,tmp:end)-act_thresh - height_cn3(i)));
    expression_span_cn3(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));  %calculate distance between those two points
    if height_cn3(i) < 0
        expression_span_cn3(i) = 0;
    end
    if isnan(signal_domain_controln3(i,1)) == 1 
        expression_span_cn3(i) = 0;
    end
    
    [height_cn5(i), tmp] = max(signal_domain_controln5(i,:)-act_thresh);    %get signal peak and subtract threshold
    height_cn5(i) = height_cn5(i)/2;    %find half height
    [~, idx1(i)] = min(abs(signal_domain_controln5(i,1:tmp)-act_thresh - height_cn5(i)));   %get the two points where signal is equal to half height
    [~, idx2(i)] = min(abs(signal_domain_controln5(i,tmp:end)-act_thresh - height_cn5(i)));
    expression_span_cn5(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));  %calculate distance between those two points
    if height_cn5(i) < 0
        expression_span_cn5(i) = 0;
    end
    if isnan(signal_domain_controln5(i,1)) == 1 
        expression_span_cn5(i) = 0;
    end
    
    [height_cn6(i), tmp] = max(signal_domain_controln6(i,:)-act_thresh);    %get signal peak and subtract threshold
    height_cn6(i) = height_cn6(i)/2;    %find half height
    [~, idx1(i)] = min(abs(signal_domain_controln6(i,1:tmp)-act_thresh - height_cn6(i)));   %get the two points where signal is equal to half height
    [~, idx2(i)] = min(abs(signal_domain_controln6(i,tmp:end)-act_thresh - height_cn6(i)));
    expression_span_cn6(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));  %calculate distance between those two points
    if height_cn6(i) < 0
        expression_span_cn6(i) = 0;
    end
    if isnan(signal_domain_controln6(i,1)) == 1 
        expression_span_cn6(i) = 0;
    end
    
    [height_bn1(i), tmp] = max(signal_domain_bestn1(i,:)-act_thresh);
    height_bn1(i) = height_bn1(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_bestn1(i,1:tmp)-act_thresh - height_bn1(i)));
    [~, idx2(i)] = min(abs(signal_domain_bestn1(i,tmp:end)-act_thresh - height_bn1(i)));
    expression_span_bn1(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_bn1(i) < 0
        expression_span_bn1(i) = 0;
    end
    if isnan(signal_domain_bestn1(i,1)) == 1 
        expression_span_bn1(i) = 0;
    end
    
    [height_bn2(i), tmp] = max(signal_domain_bestn2(i,:)-act_thresh);
    height_bn2(i) = height_bn2(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_bestn2(i,1:tmp)-act_thresh - height_bn2(i)));
    [~, idx2(i)] = min(abs(signal_domain_bestn2(i,tmp:end)-act_thresh - height_bn2(i)));
    expression_span_bn2(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_bn2(i) < 0
        expression_span_bn2(i) = 0;
    end
    if isnan(signal_domain_bestn2(i,1)) == 1 
        expression_span_bn2(i) = 0;
    end
    
    [height_bn4(i), tmp] = max(signal_domain_bestn4(i,:)-act_thresh);
    height_bn4(i) = height_bn4(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_bestn4(i,1:tmp)-act_thresh - height_bn4(i)));
    [~, idx2(i)] = min(abs(signal_domain_bestn4(i,tmp:end)-act_thresh - height_bn4(i)));
    expression_span_bn4(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_bn4(i) < 0
        expression_span_bn4(i) = 0;
    end
    if isnan(signal_domain_bestn4(i,1)) == 1 
        expression_span_bn4(i) = 0;
    end
    
    [height_dn1(i), tmp] = max(signal_domain_deln1(i,:)-act_thresh);
    height_dn1(i) = height_dn1(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_deln1(i,1:tmp)-act_thresh - height_dn1(i)));
    [~, idx2(i)] = min(abs(signal_domain_deln1(i,tmp:end)-act_thresh - height_dn1(i)));
    expression_span_dn1(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_dn1(i) < 0
        expression_span_dn1(i) = 0;
    end
    if isnan(signal_domain_deln1(i,1)) == 1 
        expression_span_dn1(i) = 0;
    end
    
    [height_dn2(i), tmp] = max(signal_domain_deln2(i,:)-act_thresh);
    height_dn2(i) = height_dn2(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_deln2(i,1:tmp)-act_thresh - height_dn2(i)));
    [~, idx2(i)] = min(abs(signal_domain_deln2(i,tmp:end)-act_thresh - height_dn2(i)));
    expression_span_dn2(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_dn2(i) < 0
        expression_span_dn2(i) = 0;
    end
    if isnan(signal_domain_deln2(i,1)) == 1 
        expression_span_dn2(i) = 0;
    end
    
    [height_dn4(i), tmp] = max(signal_domain_deln4(i,:)-act_thresh);
    height_dn4(i) = height_dn4(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_deln4(i,1:tmp)-act_thresh - height_dn4(i)));
    [~, idx2(i)] = min(abs(signal_domain_deln4(i,tmp:end)-act_thresh - height_dn4(i)));
    expression_span_dn4(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_dn4(i) < 0
        expression_span_dn4(i) = 0;
    end
    if isnan(signal_domain_deln4(i,1)) == 1 
        expression_span_dn4(i) = 0;
    end
    
    [height_dl2bn1(i), tmp] = max(signal_domain_dl2bestn1(i,:)-act_thresh);
    height_dl2bn1(i) = height_dl2bn1(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_dl2bestn1(i,1:tmp)-act_thresh - height_dl2bn1(i)));
    [~, idx2(i)] = min(abs(signal_domain_dl2bestn1(i,tmp:end)-act_thresh - height_dl2bn1(i)));
    expression_span_dl2bn1(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_dl2bn1(i) < 0
        expression_span_dl2bn1(i) = 0;
    end
    if isnan(signal_domain_dl2bestn1(i,1)) == 1 
        expression_span_dl2bn1(i) = 0;
    end
    
    [height_dl2bn2(i), tmp] = max(signal_domain_dl2bestn2(i,:)-act_thresh);
    height_dl2bn2(i) = height_dl2bn2(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_dl2bestn2(i,1:tmp)-act_thresh - height_dl2bn2(i)));
    [~, idx2(i)] = min(abs(signal_domain_dl2bestn2(i,tmp:end)-act_thresh - height_dl2bn2(i)));
    expression_span_dl2bn2(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_dl2bn2(i) < 0
        expression_span_dl2bn2(i) = 0;
    end
    if isnan(signal_domain_dl2bestn2(i,1)) == 1 
        expression_span_dl2bn2(i) = 0;
    end
    
    [height_dl2bn5(i), tmp] = max(signal_domain_dl2bestn5(i,:)-act_thresh);
    height_dl2bn5(i) = height_dl2bn5(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_dl2bestn5(i,1:tmp)-act_thresh - height_dl2bn5(i)));
    [~, idx2(i)] = min(abs(signal_domain_dl2bestn5(i,tmp:end)-act_thresh - height_dl2bn5(i)));
    expression_span_dl2bn5(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_dl2bn5(i) < 0
        expression_span_dl2bn5(i) = 0;
    end
    if isnan(signal_domain_dl2bestn5(i,1)) == 1 
        expression_span_dl2bn5(i) = 0;
    end
    
    
    [height_dl2dn1(i), tmp] = max(signal_domain_dl2deln1(i,:)-act_thresh);
    height_dl2dn1(i) = height_dl2dn1(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_dl2deln1(i,1:tmp)-act_thresh - height_dl2dn1(i)));
    [~, idx2(i)] = min(abs(signal_domain_dl2deln1(i,tmp:end)-act_thresh - height_dl2dn1(i)));
    expression_span_dl2dn1(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_dl2dn1(i) < 0
        expression_span_dl2dn1(i) = 0;
    end
    if isnan(signal_domain_dl2deln1(i,1)) == 1 
        expression_span_dl2dn1(i) = 0;
    end
    
    [height_dl2dn2(i), tmp] = max(signal_domain_dl2deln2(i,:)-act_thresh);
    height_dl2dn2(i) = height_dl2dn1(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_dl2deln2(i,1:tmp)-act_thresh - height_dl2dn2(i)));
    [~, idx2(i)] = min(abs(signal_domain_dl2deln2(i,tmp:end)-act_thresh - height_dl2dn2(i)));
    expression_span_dl2dn2(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_dl2dn2(i) < 0
        expression_span_dl2dn2(i) = 0;
    end
    if isnan(signal_domain_dl2deln2(i,1)) == 1 
        expression_span_dl2dn2(i) = 0;
    end
    
    [height_dl2dn4(i), tmp] = max(signal_domain_dl2deln4(i,:)-act_thresh);
    height_dl2dn4(i) = height_dl2dn4(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_dl2deln4(i,1:tmp)-act_thresh - height_dl2dn4(i)));
    [~, idx2(i)] = min(abs(signal_domain_dl2deln4(i,tmp:end)-act_thresh - height_dl2dn4(i)));
    expression_span_dl2dn4(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_dl2dn4(i) < 0
        expression_span_dl2dn4(i) = 0;
    end
    if isnan(signal_domain_dl2deln4(i,1)) == 1 
        expression_span_dl2dn4(i) = 0;
    end
    
    [height_zgn2(i), tmp] = max(signal_domain_zldgoodn2(i,:)-act_thresh);
    height_zgn2(i) = height_zgn2(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_zldgoodn2(i,1:tmp)-act_thresh - height_zgn2(i)));
    [~, idx2(i)] = min(abs(signal_domain_zldgoodn2(i,tmp:end)-act_thresh - height_zgn2(i)));
    expression_span_zldgoodn2(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_zgn2(i) < 0
        expression_span_zldgoodn2(i) = 0;
    end
    if isnan(signal_domain_zldgoodn2(i,1)) == 1 
        expression_span_zldgoodn2(i) = 0;
    end
    
    [height_zgn3(i), tmp] = max(signal_domain_zldgoodn3(i,:)-act_thresh);
    height_zgn3(i) = height_zgn3(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_zldgoodn3(i,1:tmp)-act_thresh - height_zgn3(i)));
    [~, idx2(i)] = min(abs(signal_domain_zldgoodn3(i,tmp:end)-act_thresh - height_zgn3(i)));
    expression_span_zldgoodn3(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_zgn3(i) < 0
        expression_span_zldgoodn3(i) = 0;
    end
    if isnan(signal_domain_zldgoodn3(i,1)) == 1 
        expression_span_zldgoodn3(i) = 0;
    end
    
    [height_zgn4(i), tmp] = max(signal_domain_zldgoodn4(i,:)-act_thresh);
    height_zgn4(i) = height_zgn4(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_zldgoodn4(i,1:tmp)-act_thresh - height_zgn4(i)));
    [~, idx2(i)] = min(abs(signal_domain_zldgoodn4(i,tmp:end)-act_thresh - height_zgn4(i)));
    expression_span_zldgoodn4(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_zgn4(i) < 0
        expression_span_zldgoodn4(i) = 0;
    end
    if isnan(signal_domain_zldgoodn4(i,1)) == 1 
        expression_span_zldgoodn4(i) = 0;
    end
    
    [height_tdn1(i), tmp] = max(signal_domain_twideln1(i,:)-act_thresh);    %get signal peak and subtract threshold
    height_tdn1(i) = height_tdn1(i)/2;    %find half height
    [~, idx1(i)] = min(abs(signal_domain_twideln1(i,1:tmp)-act_thresh - height_tdn1(i)));   %get the two points where signal is equal to half height
    [~, idx2(i)] = min(abs(signal_domain_twideln1(i,tmp:end)-act_thresh - height_tdn1(i)));
    expression_span_tdn1(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));  %calculate distance between those two points
    if height_tdn1(i) < 0
        expression_span_tdn1(i) = 0;
    end
    if isnan(signal_domain_twideln1(i,1)) == 1 
        expression_span_tdn1(i) = 0;
    end
    
    [height_tdn2(i), tmp] = max(signal_domain_twideln2(i,:)-act_thresh);    %get signal peak and subtract threshold
    height_tdn2(i) = height_tdn2(i)/2;    %find half height
    [~, idx1(i)] = min(abs(signal_domain_twideln2(i,1:tmp)-act_thresh - height_tdn2(i)));   %get the two points where signal is equal to half height
    [~, idx2(i)] = min(abs(signal_domain_twideln2(i,tmp:end)-act_thresh - height_tdn2(i)));
    expression_span_tdn2(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));  %calculate distance between those two points
    if height_tdn2(i) < 0
        expression_span_tdn2(i) = 0;
    end
    if isnan(signal_domain_twideln2(i,1)) == 1 
        expression_span_tdn2(i) = 0;
    end
    
    
    %Revision
    [height_Rcn1(i), tmp] = max(signal_domain_Rcontroln1(i,:)-act_threshR);    %get signal peak and subtract threshold
    height_Rcn1(i) = height_Rcn1(i)/2;    %find half height
    [~, idx1(i)] = min(abs(signal_domain_Rcontroln1(i,1:tmp)-act_threshR - height_Rcn1(i)));   %get the two points where signal is equal to half height
    [~, idx2(i)] = min(abs(signal_domain_Rcontroln1(i,tmp:end)-act_threshR - height_Rcn1(i)));
    expression_span_Rcn1(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));  %calculate distance between those two points
    if height_Rcn1(i) < 0
        expression_span_Rcn1(i) = 0;
    end
    if isnan(signal_domain_Rcontroln1(i,1)) == 1 
        expression_span_Rcn1(i) = 0;
    end
    
    [height_Rcn2(i), tmp] = max(signal_domain_Rcontroln2(i,:)-act_threshR);    %get signal peak and subtract threshold
    height_Rcn2(i) = height_Rcn2(i)/2;    %find half height
    [~, idx1(i)] = min(abs(signal_domain_Rcontroln2(i,1:tmp)-act_threshR - height_Rcn2(i)));   %get the two points where signal is equal to half height
    [~, idx2(i)] = min(abs(signal_domain_Rcontroln2(i,tmp:end)-act_threshR - height_Rcn2(i)));
    expression_span_Rcn2(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));  %calculate distance between those two points
    if height_Rcn2(i) < 0
        expression_span_Rcn2(i) = 0;
    end
    if isnan(signal_domain_Rcontroln2(i,1)) == 1 
        expression_span_Rcn2(i) = 0;
    end
    
    [height_Rcn3(i), tmp] = max(signal_domain_Rcontroln3(i,:)-act_threshR);    %get signal peak and subtract threshold
    height_Rcn3(i) = height_Rcn3(i)/2;    %find half height
    [~, idx1(i)] = min(abs(signal_domain_Rcontroln3(i,1:tmp)-act_threshR - height_Rcn3(i)));   %get the two points where signal is equal to half height
    [~, idx2(i)] = min(abs(signal_domain_Rcontroln3(i,tmp:end)-act_threshR - height_Rcn3(i)));
    expression_span_Rcn3(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));  %calculate distance between those two points
    if height_Rcn3(i) < 0
        expression_span_Rcn3(i) = 0;
    end
    if isnan(signal_domain_Rcontroln3(i,1)) == 1 
        expression_span_Rcn3(i) = 0;
    end
    
    [height_Rdzn1(i), tmp] = max(signal_domain_Rdzn1(i,:)-act_threshR);    %get signal peak and subtract threshold
    height_Rdzn1(i) = height_Rdzn1(i)/2;    %find half height
    [~, idx1(i)] = min(abs(signal_domain_Rdzn1(i,1:tmp)-act_threshR - height_Rdzn1(i)));   %get the two points where signal is equal to half height
    [~, idx2(i)] = min(abs(signal_domain_Rdzn1(i,tmp:end)-act_threshR - height_Rdzn1(i)));
    expression_span_Rdzn1(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));  %calculate distance between those two points
    if height_Rdzn1(i) < 0
        expression_span_Rdzn1(i) = 0;
    end
    if isnan(signal_domain_Rdzn1(i,1)) == 1 
        expression_span_Rdzn1(i) = 0;
    end
    
    [height_Rdzn2(i), tmp] = max(signal_domain_Rdzn2(i,:)-act_threshR);    %get signal peak and subtract threshold
    height_Rdzn2(i) = height_Rdzn2(i)/2;    %find half height
    [~, idx1(i)] = min(abs(signal_domain_Rdzn2(i,1:tmp)-act_threshR - height_Rdzn2(i)));   %get the two points where signal is equal to half height
    [~, idx2(i)] = min(abs(signal_domain_Rdzn2(i,tmp:end)-act_threshR - height_Rdzn2(i)));
    expression_span_Rdzn2(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));  %calculate distance between those two points
    if height_Rdzn2(i) < 0
        expression_span_Rdzn2(i) = 0;
    end
    if isnan(signal_domain_Rdzn2(i,1)) == 1 
        expression_span_Rdzn2(i) = 0;
    end
    
    [height_Rdzn3(i), tmp] = max(signal_domain_Rdzn3(i,:)-act_threshR);    %get signal peak and subtract threshold
    height_Rdzn3(i) = height_Rdzn3(i)/2;    %find half height
    [~, idx1(i)] = min(abs(signal_domain_Rdzn3(i,1:tmp)-act_threshR - height_Rdzn3(i)));   %get the two points where signal is equal to half height
    [~, idx2(i)] = min(abs(signal_domain_Rdzn3(i,tmp:end)-act_threshR - height_Rdzn3(i)));
    expression_span_Rdzn3(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));  %calculate distance between those two points
    if height_Rdzn3(i) < 0
        expression_span_Rdzn3(i) = 0;
    end
    if isnan(signal_domain_Rdzn3(i,1)) == 1 
        expression_span_Rdzn3(i) = 0;
    end
    
    [height_Rdn1(i), tmp] = max(signal_domain_Rdn1(i,:)-act_threshR);    %get signal peak and subtract threshold
    height_Rdn1(i) = height_Rdn1(i)/2;    %find half height
    [~, idx1(i)] = min(abs(signal_domain_Rdn1(i,1:tmp)-act_threshR - height_Rdn1(i)));   %get the two points where signal is equal to half height
    [~, idx2(i)] = min(abs(signal_domain_Rdn1(i,tmp:end)-act_threshR - height_Rdn1(i)));
    expression_span_Rdn1(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));  %calculate distance between those two points
    if height_Rdn1(i) < 0
        expression_span_Rdn1(i) = 0;
    end
    if isnan(signal_domain_Rdn1(i,1)) == 1 
        expression_span_Rdn1(i) = 0;
    end
    
    [height_Rdn2(i), tmp] = max(signal_domain_Rdn2(i,:)-act_threshR);    %get signal peak and subtract threshold
    height_Rdn2(i) = height_Rdn2(i)/2;    %find half height
    [~, idx1(i)] = min(abs(signal_domain_Rdn2(i,1:tmp)-act_threshR - height_Rdn2(i)));   %get the two points where signal is equal to half height
    [~, idx2(i)] = min(abs(signal_domain_Rdn2(i,tmp:end)-act_threshR - height_Rdn2(i)));
    expression_span_Rdn2(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));  %calculate distance between those two points
    if height_Rdn2(i) < 0
        expression_span_Rdn2(i) = 0;
    end
    if isnan(signal_domain_Rdn2(i,1)) == 1 
        expression_span_Rdn2(i) = 0;
    end
    
    [height_Rdn3(i), tmp] = max(signal_domain_Rdn3(i,:)-act_threshR);    %get signal peak and subtract threshold
    height_Rdn3(i) = height_Rdn3(i)/2;    %find half height
    [~, idx1(i)] = min(abs(signal_domain_Rdn3(i,1:tmp)-act_threshR - height_Rdn3(i)));   %get the two points where signal is equal to half height
    [~, idx2(i)] = min(abs(signal_domain_Rdn3(i,tmp:end)-act_threshR - height_Rdn3(i)));
    expression_span_Rdn3(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));  %calculate distance between those two points
    if height_Rdzn3(i) < 0
        expression_span_Rdn3(i) = 0;
    end
    if isnan(signal_domain_Rdn3(i,1)) == 1 
        expression_span_Rdn3(i) = 0;
    end
    
%{
    [height_zwn1(i), tmp] = max(signal_domain_zldweakn1(i,:)-act_thresh);
    height_zwn1(i) = height_zwn1(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_zldweakn1(i,1:tmp)-act_thresh - height_zwn1(i)));
    [~, idx2(i)] = min(abs(signal_domain_zldweakn1(i,tmp:end)-act_thresh - height_zwn1(i)));
    expression_span_zldweakn1(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_zwn1(i) < 0
        expression_span_zldweakn1(i) = 0;
    end
    if isnan(signal_domain_zldweakn1(i,1)) == 1 
        expression_span_zldweakn1(i) = 0;
    end
    
    [height_zwn2(i), tmp] = max(signal_domain_zldweakn2(i,:)-act_thresh);
    height_zwn2(i) = height_zwn2(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_zldweakn2(i,1:tmp)-act_thresh - height_zwn2(i)));
    [~, idx2(i)] = min(abs(signal_domain_zldweakn2(i,tmp:end)-act_thresh - height_zwn2(i)));
    expression_span_zldweakn2(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_zwn2(i) < 0
        expression_span_zldweakn2(i) = 0;
    end
    if isnan(signal_domain_zldweakn2(i,1)) == 1 
        expression_span_zldweakn2(i) = 0;
    end
    
    [height_zwn5(i), tmp] = max(signal_domain_zldweakn5(i,:)-act_thresh);
    height_zwn5(i) = height_zwn5(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_zldweakn5(i,1:tmp)-act_thresh - height_zwn5(i)));
    [~, idx2(i)] = min(abs(signal_domain_zldweakn5(i,tmp:end)-act_thresh - height_zwn5(i)));
    expression_span_zldweakn5(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_zwn5(i) < 0
        expression_span_zldweakn5(i) = 0;
    end
    if isnan(signal_domain_zldweakn5(i,1)) == 1 
        expression_span_zldweakn5(i) = 0;
    end
    
    [height_mcphet(i), tmp] = max(signal_domain_mcphet(i,:)-act_thresh);
    height_mcphet(i) = height_mcphet(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_mcphet(i,1:tmp)-act_thresh - height_mcphet(i)));
    [~, idx2(i)] = min(abs(signal_domain_mcphet(i,tmp:end)-act_thresh - height_mcphet(i)));
    expression_span_mcphet(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_zgn2(i) < 0
        expression_span_mcphet(i) = 0;
    end
    if isnan(signal_domain_mcphet(i,1)) == 1 
        expression_span_mcphet(i) = 0;
    end
    
    [height_dlhetcontrol(i), tmp] = max(signal_domain_dlhetcontrol(i,:)-act_thresh);
    height_dlhetcontrol(i) = height_dlhetcontrol(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_dlhetcontrol(i,1:tmp)-act_thresh - height_dlhetcontrol(i)));
    [~, idx2(i)] = min(abs(signal_domain_dlhetcontrol(i,tmp:end)-act_thresh - height_dlhetcontrol(i)));
    expression_span_dlhetcontrol(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_zgn2(i) < 0
        expression_span_dlhetcontrol(i) = 0;
    end
    if isnan(signal_domain_dlhetcontrol(i,1)) == 1 
        expression_span_dlhetcontrol(i) = 0;
    end
    
    [height_dlhetdl2best(i), tmp] = max(signal_domain_dlhetdl2best(i,:)-act_thresh);
    height_dlhetdl2best(i) = height_dlhetdl2best(i)/2;
    [~, idx1(i)] = min(abs(signal_domain_dlhetdl2best(i,1:tmp)-act_thresh - height_dlhetdl2best(i)));
    [~, idx2(i)] = min(abs(signal_domain_dlhetdl2best(i,tmp:end)-act_thresh - height_dlhetdl2best(i)));
    expression_span_dlhetdl2best(i) = abs(domain(tmp-1+idx2(i))-domain(idx1(i)));
    if height_zgn2(i) < 0
        expression_span_dlhetdl2best(i) = 0;
    end
    if isnan(signal_domain_dlhetdl2best(i,1)) == 1 
        expression_span_dlhetdl2best(i) = 0;
    end
%}
    %Ave Output Domain
%     domain_control(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Area_control(i,:),domain);
%     domain_best(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Area_best(i,:),domain);
%     domain_del(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Area_del(i,:),domain);
%     domain_dl2best(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Area_dl2best(i,:),domain);
%     domain_dl2del(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Area_dl2del(i,:),domain);
%     domain_zldgood(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Area_zldgood(i,:),domain);
%     domain_zldweak(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Area_zldweak(i,:),domain);
%     domain_mcphet(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Area_mcphet(i,:),domain);
%     domain_dlhetcontrol(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Area_dlhetcontrol(i,:),domain);
%     domain_dlhetdl2best(i,:) = interp1([0:bins_bestn1-1]/bins_bestn1,ave_Area_dlhetdl2best(i,:),domain);
end

expression_span_c = (expression_span_cn1 + expression_span_cn2 + expression_span_cn4 + expression_span_cn3 + expression_span_cn5 + expression_span_cn6)./6;
expression_span_b = (expression_span_bn1 + expression_span_bn2 + expression_span_bn4)./rep;
expression_span_d = (expression_span_dn1 + expression_span_dn2 + expression_span_dn4)./rep;
expression_span_dl2b = (expression_span_dl2bn1 + expression_span_dl2bn2 + expression_span_dl2bn5)./rep;
expression_span_dl2d = (expression_span_dl2dn1 + expression_span_dl2dn2 + expression_span_dl2dn4)./rep;
expression_span_zldgood = (expression_span_zldgoodn2 + expression_span_zldgoodn3 + expression_span_zldgoodn4)./rep;
% expression_span_zldweak = (expression_span_zldweakn1 + expression_span_zldweakn2 + expression_span_zldweakn5)./rep;

expression_span_Rc = (expression_span_Rcn1 + expression_span_Rcn2 + expression_span_Rcn3)./3;
expression_span_Rdz = (expression_span_Rdzn1 + expression_span_Rdzn2 + expression_span_Rdzn3)./3;


std_expression_c = std([expression_span_cn1; expression_span_cn2; expression_span_cn4; expression_span_cn3; expression_span_cn5; expression_span_cn6]);
std_expression_b = std([expression_span_bn1; expression_span_bn2; expression_span_bn4]);
std_expression_d = std([expression_span_dn1; expression_span_dn2; expression_span_dn4]);
std_expression_dl2b = std([expression_span_dl2bn1; expression_span_dl2bn2; expression_span_dl2bn5]);
std_expression_dl2d = std([expression_span_dl2dn1; expression_span_dl2dn2; expression_span_dl2dn4]);
std_expression_zldgood = std([expression_span_zldgoodn2; expression_span_zldgoodn3; expression_span_zldgoodn4]);

% std_expression_zldweak = std([expression_span_zldweakn1; expression_span_zldweakn2; expression_span_zldweakn5]);
std_expression_Rc = std([expression_span_Rcn1; expression_span_Rcn2; expression_span_Rcn3]);
std_expression_Rdz = std([expression_span_Rdzn1; expression_span_Rdzn2; expression_span_Rdzn3]);


figure;
errorbar(norm_time,expression_span_c, std_expression_c,'Color', rgb('blue')); hold on
errorbar(norm_time,expression_span_b, std_expression_b,'Color', rgb('black')); 
errorbar(norm_time,expression_span_d, std_expression_d,'Color', rgb('red')); 
errorbar(norm_time,expression_span_dl2b, std_expression_dl2b,'Color', rgb('darkorange')); 
errorbar(norm_time,expression_span_dl2d, std_expression_dl2d,'Color', rgb('purple')); 
errorbar(norm_time,expression_span_zldgood, std_expression_zldgood,'Color', rgb('green')); 
% errorbar(norm_time,expression_span_zldweak, std_expression_zldweak,'Color', rgb('brown')); 
plot(norm_time,expression_span_c,'Color', rgb('blue'),'LineWidth', 2); hold on
plot(norm_time,expression_span_b,'Color', rgb('black'),'LineWidth', 2); 
plot(norm_time,expression_span_d,'Color', rgb('red'),'LineWidth', 2); 
plot(norm_time,expression_span_dl2b,'Color', rgb('darkorange'),'LineWidth', 2); 
plot(norm_time,expression_span_dl2d,'Color', rgb('purple'),'LineWidth', 2); 
plot(norm_time,expression_span_zldgood,'Color', rgb('green'),'LineWidth', 2);
% plot(norm_time,expression_span_zldweak,'Color', rgb('brown'),'LineWidth', 2); 
title(['Expression Domain Over Time'])
xlabel('t''')
ylabel('Expression Domain (AU)')
legend('control', 'dl1 strong','dl1 weak', 'dl2 strong', 'dl2 weak','zld good','zld weak','Location', 'northwest')
ylim([0 .65])

%half dorsal
% figure;
% plot(norm_time,expression_span_mcphet,'Color', rgb('blue')); hold on
% plot(norm_time,expression_span_dlhetcontrol,'Color', rgb('black')); 
% plot(norm_time,expression_span_dlhetdl2best,'Color', rgb('red')); 
% title(['Expression Domain Over Time'])
% xlabel('t''')
% ylabel('Expression Domain (AU)')
% legend('mcp control', 'dl het control','dl het dl2 best','Location', 'northwest')
% 

%{ 
%Poster Plots
figure
shadedErrorBar(norm_time,expression_span_c, std_expression_c,{'Color', rgb('black'),'LineWidth',2}); hold on
shadedErrorBar(norm_time,expression_span_b, std_expression_b,{'Color', rgb('red'),'LineWidth',2}); 
shadedErrorBar(norm_time,expression_span_d, std_expression_d,{'--','Color', rgb('lightcoral'),'LineWidth',2}); 
title(['Expression Domain Over Time'])
xlabel('t''')
ylabel('Expression Domain (AU)')
legend('control', 'dl1 strong','dl1 weak','Location','northwest')
ylim([0 .7])

figure
shadedErrorBar(norm_time,expression_span_c, std_expression_c,{'Color', rgb('black'),'LineWidth',2}); hold on
shadedErrorBar(norm_time,expression_span_dl2b, std_expression_dl2b,{'Color', rgb('blue'),'LineWidth',2}); 
shadedErrorBar(norm_time,expression_span_dl2d, std_expression_dl2d,{'--','Color', rgb('cornflowerblue'),'LineWidth',2}); 
title(['Expression Domain Over Time'])
xlabel('t''')
ylabel('Expression Domain (AU)')
legend('control', 'dl2 strong', 'dl2 weak','Location', 'northwest')
ylim([0 .7])

figure
shadedErrorBar(norm_time,expression_span_c, std_expression_c,{'Color', rgb('black'),'LineWidth',2}); hold on
shadedErrorBar(norm_time,expression_span_zldgood, std_expression_zldgood,{'Color', rgb('green'),'LineWidth',2}); 
title(['Expression Domain Over Time'])
xlabel('t''')
ylabel('Expression Domain (AU)')
legend('control','zld strong','Location', 'northwest')
ylim([0 .7])


%Revision
figure
shadedErrorBar(norm_time,expression_span_c, std_expression_c,{'Color', rgb('gray'),'LineWidth',2}); hold on
shadedErrorBar(norm_time,expression_span_Rc, std_expression_Rc,{'Color', rgb('black'),'LineWidth',2}); hold on
shadedErrorBar(norm_time,expression_span_Rdz, std_expression_Rdz,{'Color', rgb('red'),'LineWidth',2}); hold on
shadedErrorBar(norm_time,expression_span_Rd, std_expression_Rd,{'Color', rgb('blue'),'LineWidth',2}); hold on

title(['Expression Domain Over Time'])
xlabel('t''')
ylabel('Expression Domain (AU)')
legend('control', 'Revision Control','Location','northwest')
ylim([0 .7])
%}
%% Active Nuclei Position

%Rep 1
for i = 1:length(activeM_controln1)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_controln1(i) == [shift_bin_controln1{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_controln1{j}) + position;
        if tmp <= position
            loc_controln1(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 2
for i = 1:length(activeM_controln2)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_controln2(i) == [shift_bin_controln2{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_controln2{j}) + position;
        if tmp <= position
            loc_controln2(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 3
for i = 1:length(activeM_controln4)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_controln4(i) == [shift_bin_controln4{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_controln4{j}) + position;
        if tmp <= position
            loc_controln4(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 4
for i = 1:length(activeM_controln3)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_controln3(i) == [shift_bin_controln3{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_controln3{j}) + position;
        if tmp <= position
            loc_controln3(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 5
for i = 1:length(activeM_controln5)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_controln5(i) == [shift_bin_controln5{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_controln5{j}) + position;
        if tmp <= position
            loc_controln5(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 6
for i = 1:length(activeM_controln6)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_controln6(i) == [shift_bin_controln6{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_controln6{j}) + position;
        if tmp <= position
            loc_controln6(i) = (j-1)/16;
            check = 1;
        end
    end
end
lowSignal_controln1 = find(active_Signal_c{1} <= 500);
medSignal_controln1 = find(500 < active_Signal_c{1} <= 1000);
highSignal_controln1 = find(active_Signal_c{1} > 1000);
lowSignal_controln2 = find(active_Signal_c{2} <= 500);
medSignal_controln2 = find(500 < active_Signal_c{2} <= 1000);
highSignal_controln2 = find(active_Signal_c{2} > 1000);
lowSignal_controln4 = find(active_Signal_c{3} <= 500);
medSignal_controln4 = find(500 < active_Signal_c{3} <= 1000);
highSignal_controln4 = find(active_Signal_c{3} > 1000);


%Rep 1
for i = 1:length(activeM_bestn1)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_bestn1(i) == [shift_bin_bestn1{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_bestn1{j}) + position;
        if tmp <= position
            loc_bestn1(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 2
for i = 1:length(activeM_bestn2)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_bestn2(i) == [shift_bin_bestn2{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_bestn2{j}) + position;
        if tmp <= position
            loc_bestn2(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 3
for i = 1:length(activeM_bestn4)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_bestn4(i) == [shift_bin_bestn4{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_bestn4{j}) + position;
        if tmp <= position
            loc_bestn4(i) = (j-1)/16;
            check = 1;
        end
    end
end

lowSignal_bestn1 = find(active_Signal_b{1} <= 500);
medSignal_bestn1 = find(500 < active_Signal_b{1} <= 1000);
highSignal_bestn1 = find(active_Signal_b{1} > 1000);
lowSignal_bestn2 = find(active_Signal_b{2} <= 500);
medSignal_bestn2 = find(500 < active_Signal_b{2} <= 1000);
highSignal_bestn2 = find(active_Signal_b{2} > 1000);
lowSignal_bestn4 = find(active_Signal_b{3} <= 500);
medSignal_bestn4 = find(500 < active_Signal_b{3} <= 1000);
highSignal_bestn4 = find(active_Signal_b{3} > 1000);


%Rep 1
for i = 1:length(activeM_deln1)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_deln1(i) == [shift_bin_deln1{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_deln1{j}) + position;
        if tmp <= position
            loc_deln1(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 2
for i = 1:length(activeM_deln2)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_deln2(i) == [shift_bin_deln2{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_deln2{j}) + position;
        if tmp <= position
            loc_deln2(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 3
for i = 1:length(activeM_deln4)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_deln4(i) == [shift_bin_deln4{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_deln4{j}) + position;
        if tmp <= position
            loc_deln4(i) = (j-1)/16;
            check = 1;
        end
    end
end

lowSignal_deln1 = find(active_Signal_d{1} <= 500);
medSignal_deln1 = find(500 < active_Signal_d{1} <= 1000);
highSignal_deln1 = find(active_Signal_d{1} > 1000);
lowSignal_deln2 = find(active_Signal_d{2} <= 500);
medSignal_deln2 = find(500 < active_Signal_d{2} <= 1000);
highSignal_deln2 = find(active_Signal_d{2} > 1000);
lowSignal_deln4 = find(active_Signal_d{3} <= 500);
medSignal_deln4 = find(500 < active_Signal_d{3} <= 1000);
highSignal_deln4 = find(active_Signal_d{3} > 1000);

%Rep 1
for i = 1:length(activeM_dl2bestn1)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_dl2bestn1(i) == [shift_bin_dl2bestn1{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_dl2bestn1{j}) + position;
        if tmp <= position
            loc_dl2bestn1(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 2
for i = 1:length(activeM_dl2bestn2)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_dl2bestn2(i) == [shift_bin_dl2bestn2{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_dl2bestn2{j}) + position;
        if tmp <= position
            loc_dl2bestn2(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 3
for i = 1:length(activeM_dl2bestn5)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_dl2bestn5(i) == [shift_bin_dl2bestn5{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_dl2bestn5{j}) + position;
        if tmp <= position
            loc_dl2bestn5(i) = (j-1)/16;
            check = 1;
        end
    end
end

lowSignal_dl2bestn1 = find(active_Signal_dl2b{1} <= 500);
medSignal_dl2bestn1 = find(500 < active_Signal_dl2b{1} <= 1000);
highSignal_dl2bestn1 = find(active_Signal_dl2b{1} > 1000);
lowSignal_dl2bestn2 = find(active_Signal_dl2b{2} <= 500);
medSignal_dl2bestn2 = find(500 < active_Signal_dl2b{2} <= 1000);
highSignal_dl2bestn2 = find(active_Signal_dl2b{2} > 1000);
lowSignal_dl2bestn5 = find(active_Signal_dl2b{3} <= 500);
medSignal_dl2bestn5 = find(500 < active_Signal_dl2b{3} <= 1000);
highSignal_dl2bestn5 = find(active_Signal_dl2b{3} > 1000);


%Rep 1
for i = 1:length(activeM_dl2deln1)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_dl2deln1(i) == [shift_bin_dl2deln1{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_dl2deln1{j}) + position;
        if tmp <= position
            loc_dl2deln1(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 2
for i = 1:length(activeM_dl2deln2)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_dl2deln2(i) == [shift_bin_dl2deln2{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_dl2deln2{j}) + position;
        if tmp <= position
            loc_dl2deln2(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 3
for i = 1:length(activeM_dl2deln4)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_dl2deln4(i) == [shift_bin_dl2deln4{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_dl2deln4{j}) + position;
        if tmp <= position
            loc_dl2deln4(i) = (j-1)/16;
            check = 1;
        end
    end
end

lowSignal_dl2deln1 = find(active_Signal_dl2d{1} <= 500);
medSignal_dl2deln1 = find(500 < active_Signal_dl2d{1} <= 1000);
highSignal_dl2deln1 = find(active_Signal_dl2d{1} > 1000);
lowSignal_dl2deln2 = find(active_Signal_dl2d{2} <= 500);
medSignal_dl2deln2 = find(500 < active_Signal_dl2d{2} <= 1000);
highSignal_dl2deln2 = find(active_Signal_dl2d{2} > 1000);
lowSignal_dl2deln4 = find(active_Signal_dl2d{3} <= 500);
medSignal_dl2deln4 = find(500 < active_Signal_dl2d{3} <= 1000);
highSignal_dl2deln4 = find(active_Signal_dl2d{3} > 1000);


%Rep 1
for i = 1:length(activeM_zldgoodn2)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_zldgoodn2(i) == [shift_bin_zldgoodn2{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_zldgoodn2{j}) + position;
        if tmp <= position
            loc_zldgoodn2(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 2
for i = 1:length(activeM_zldgoodn3)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_zldgoodn3(i) == [shift_bin_zldgoodn3{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_zldgoodn3{j}) + position;
        if tmp <= position
            loc_zldgoodn3(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 3
for i = 1:length(activeM_zldgoodn4)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_zldgoodn4(i) == [shift_bin_zldgoodn4{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_zldgoodn4{j}) + position;
        if tmp <= position
            loc_zldgoodn4(i) = (j-1)/16;
            check = 1;
        end
    end
end

lowSignal_zldgoodn2 = find(active_Signal_zldgood{1} <= 500);
medSignal_zldgoodn2 = find(500 < active_Signal_zldgood{1} <= 1000);
highSignal_zldgoodn2 = find(active_Signal_zldgood{1} > 1000);
lowSignal_zldgoodn3 = find(active_Signal_zldgood{2} <= 500);
medSignal_zldgoodn3 = find(500 < active_Signal_zldgood{2} <= 1000);
highSignal_zldgoodn3 = find(active_Signal_zldgood{2} > 1000);
lowSignal_zldgoodn4 = find(active_Signal_zldgood{3} <= 500);
medSignal_zldgoodn4 = find(500 < active_Signal_zldgood{3} <= 1000);
highSignal_zldgoodn4 = find(active_Signal_zldgood{3} > 1000);

%Rep 1
for i = 1:length(activeM_zldweakn1)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_zldweakn1(i) == [shift_bin_zldweakn1{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_zldweakn1{j}) + position;
        if tmp <= position
            loc_zldweakn1(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 2
for i = 1:length(activeM_zldweakn2)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_zldweakn2(i) == [shift_bin_zldweakn2{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_zldweakn2{j}) + position;
        if tmp <= position
            loc_zldweakn2(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 3
for i = 1:length(activeM_zldweakn5)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_zldweakn5(i) == [shift_bin_zldweakn5{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_zldweakn5{j}) + position;
        if tmp <= position
            loc_zldweakn5(i) = (j-1)/16;
            check = 1;
        end
    end
end

lowSignal_zldweakn1 = find(active_Signal_zldweak{1} <= 500);
medSignal_zldweakn1 = find(500 < active_Signal_zldweak{1} <= 1000);
highSignal_zldweakn1 = find(active_Signal_zldweak{1} > 1000);
lowSignal_zldweakn2 = find(active_Signal_zldweak{2} <= 500);
medSignal_zldweakn2 = find(500 < active_Signal_zldweak{2} <= 1000);
highSignal_zldweakn2 = find(active_Signal_zldweak{2} > 1000);
lowSignal_zldweakn5 = find(active_Signal_zldweak{3} <= 500);
medSignal_zldweakn5 = find(500 < active_Signal_zldweak{3} <= 1000);
highSignal_zldweakn5 = find(active_Signal_zldweak{3} > 1000);

for i = 1:length(activeM_twideln1)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_twideln1(i) == [shift_bin_twideln1{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_twideln1{j}) + position;
        if tmp <= position
            loc_twideln1(i) = (j-1)/16;
            check = 1;
        end
    end
end

for i = 1:length(activeM_twideln2)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_twideln2(i) == [shift_bin_twideln2{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_twideln2{j}) + position;
        if tmp <= position
            loc_twideln2(i) = (j-1)/16;
            check = 1;
        end
    end
end

%Revision
%Rep 1
for i = 1:length(activeM_Rcn1)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_Rcn1(i) == [shift_bin_Rcn1{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_Rcn1{j}) + position;
        if tmp <= position
            loc_Rcontroln1(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 2
for i = 1:length(activeM_Rcn2)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_Rcn2(i) == [shift_bin_Rcn2{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_Rcn2{j}) + position;
        if tmp <= position
            loc_Rcontroln2(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 3
for i = 1:length(activeM_Rcn3)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_Rcn3(i) == [shift_bin_Rcn3{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_Rcn3{j}) + position;
        if tmp <= position
            loc_Rcontroln3(i) = (j-1)/16;
            check = 1;
        end
    end
end

%Rep 1
for i = 1:length(activeM_Rdzn1)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_Rdzn1(i) == [shift_bin_Rdzn1{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_Rdzn1{j}) + position;
        if tmp <= position
            loc_Rdzn1(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 2
for i = 1:length(activeM_Rdzn2)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_Rdzn2(i) == [shift_bin_Rdzn2{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_Rdzn2{j}) + position;
        if tmp <= position
            loc_Rdzn2(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 3
for i = 1:length(activeM_Rdzn3)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_Rdzn3(i) == [shift_bin_Rdzn3{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_Rdzn3{j}) + position;
        if tmp <= position
            loc_Rdzn3(i) = (j-1)/16;
            check = 1;
        end
    end
end

%Rep 1
for i = 1:length(activeM_Rdn1)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_Rdn1(i) == [shift_bin_Rdn1{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_Rdn1{j}) + position;
        if tmp <= position
            loc_Rdn1(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 2
for i = 1:length(activeM_Rdn2)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_Rdn2(i) == [shift_bin_Rdn2{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_Rdn2{j}) + position;
        if tmp <= position
            loc_Rdn2(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 3
for i = 1:length(activeM_Rdn3)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_Rdn3(i) == [shift_bin_Rdn3{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_Rdn3{j}) + position;
        if tmp <= position
            loc_Rdn3(i) = (j-1)/16;
            check = 1;
        end
    end
end

figure;
subplot(3,3,1)
histogram(loc_dl2bestn1(lowSignal_dl2bestn1),16);
title('Low Signal Zld Good - Rep 1')
xlim([0 1])
ylim([0 100])
subplot(3,3,2)
histogram(loc_dl2bestn1(medSignal_dl2bestn1),16);
title('Mid Signal Zld Good - Rep 1')
xlim([0 1])
ylim([0 100])
subplot(3,3,3)
histogram(loc_dl2bestn1(highSignal_dl2bestn1),16);
title('High Signal Zld Good - Rep 1')
xlim([0 1])
ylim([0 100])

subplot(3,3,4)
histogram(loc_zldweakn2(lowSignal_zldweakn2),16);
title('Low Signal Zld Good - Rep 2')
xlim([0 1])
ylim([0 100])
subplot(3,3,5)
histogram(loc_zldweakn2(medSignal_zldweakn2),16);
title('Mid Signal Zld Good - Rep 2')
xlim([0 1])
ylim([0 100])
subplot(3,3,6)
histogram(loc_zldweakn2(highSignal_zldweakn2),16);
title('High Signal Zld Good - Rep 2')
xlim([0 1])
ylim([0 100])

subplot(3,3,7)
histogram(loc_zldweakn5(lowSignal_zldweakn5),16);
title('Low Signal Zld Good - Rep 3')
xlim([0 1])
ylim([0 100])
subplot(3,3,8)
histogram(loc_zldweakn5(medSignal_zldweakn5),16);
title('Mid Signal Zld Good - Rep 3')
xlim([0 1])
ylim([0 100])
subplot(3,3,9)
histogram(loc_zldweakn5(highSignal_zldweakn5),16);
title('High Signal Zld Good - Rep 3')
xlim([0 1])
ylim([0 100])

%% Active Nuclei Position HET

%Rep 1
for i = 1:length(activeM_hetcontroln1)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_hetcontroln1(i) == [shift_bin_hetcontroln1{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_hetcontroln1{j}) + position;
        if tmp <= position
            loc_hetcontroln1(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 2
for i = 1:length(activeM_hetcontroln2)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_hetcontroln2(i) == [shift_bin_hetcontroln2{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_hetcontroln2{j}) + position;
        if tmp <= position
            loc_hetcontroln2(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 3
for i = 1:length(activeM_hetcontroln3)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_hetcontroln3(i) == [shift_bin_hetcontroln3{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_hetcontroln3{j}) + position;
        if tmp <= position
            loc_hetcontroln3(i) = (j-1)/16;
            check = 1;
        end
    end
end


%Rep 1
for i = 1:length(activeM_hetbestn1)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_hetbestn1(i) == [shift_bin_hetbestn1{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_hetbestn1{j}) + position;
        if tmp <= position
            loc_hetbestn1(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 2
for i = 1:length(activeM_hetbestn2)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_hetbestn2(i) == [shift_bin_hetbestn2{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_hetbestn2{j}) + position;
        if tmp <= position
            loc_hetbestn2(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 3
for i = 1:length(activeM_hetbestn3)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_hetbestn3(i) == [shift_bin_hetbestn3{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_hetbestn3{j}) + position;
        if tmp <= position
            loc_hetbestn3(i) = (j-1)/16;
            check = 1;
        end
    end
end


%Rep 1
for i = 1:length(activeM_hetdl2bestn1)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_hetdl2bestn1(i) == [shift_bin_hetdl2bestn1{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_hetdl2bestn1{j}) + position;
        if tmp <= position
            loc_hetdl2bestn1(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 2
for i = 1:length(activeM_hetdl2bestn2)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_hetdl2bestn2(i) == [shift_bin_hetdl2bestn2{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_hetdl2bestn2{j}) + position;
        if tmp <= position
            loc_hetdl2bestn2(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 3
for i = 1:length(activeM_hetdl2bestn3)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_hetdl2bestn3(i) == [shift_bin_hetdl2bestn3{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_hetdl2bestn3{j}) + position;
        if tmp <= position
            loc_hetdl2bestn3(i) = (j-1)/16;
            check = 1;
        end
    end
end


%Rep 1
for i = 1:length(activeM_hetdl2deln1)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_hetdl2deln1(i) == [shift_bin_hetdl2deln1{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_hetdl2deln1{j}) + position;
        if tmp <= position
            loc_hetdl2deln1(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 2
for i = 1:length(activeM_hetdl2deln2)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_hetdl2deln2(i) == [shift_bin_hetdl2deln2{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_hetdl2deln2{j}) + position;
        if tmp <= position
            loc_hetdl2deln2(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 3
for i = 1:length(activeM_hetdl2deln3)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_hetdl2deln3(i) == [shift_bin_hetdl2deln3{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_hetdl2deln3{j}) + position;
        if tmp <= position
            loc_hetdl2deln3(i) = (j-1)/16;
            check = 1;
        end
    end
end


%Rep 1
for i = 1:length(activeM_hetzldbestn1)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_hetzldbestn1(i) == [shift_bin_hetzldbestn1{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_hetzldbestn1{j}) + position;
        if tmp <= position
            loc_hetzldbestn1(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 2
for i = 1:length(activeM_hetzldbestn2)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_hetzldbestn2(i) == [shift_bin_hetzldbestn2{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_hetzldbestn2{j}) + position;
        if tmp <= position
            loc_hetzldbestn2(i) = (j-1)/16;
            check = 1;
        end
    end
end
%Rep 3
for i = 1:length(activeM_hetzldbestn3)
    check = 0;
    j = 0;
    position = 0;
    tmp = find(activeM_hetzldbestn3(i) == [shift_bin_hetzldbestn3{:}]);
    while check == 0
        j = j+1;
        position = length(shift_bin_hetzldbestn3{j}) + position;
        if tmp <= position
            loc_hetzldbestn3(i) = (j-1)/16;
            check = 1;
        end
    end
end



%% Bursting Quantification

for i = 1:length(activeM_controln1)
    time_active_cn1{i} = find(M_controln1(:,activeM_controln1(i))>= .8*thresh_base_controln1);
    burst_cn1(i) = length(time_active_cn1{i})/(length(T_controln1) - time_active_cn1{i}(1));
end
for i = 1:length(activeM_controln2)
    time_active_cn2{i} = find(M_controln2(:,activeM_controln2(i))>= .8*thresh_base_controln2);
    burst_cn2(i) = length(time_active_cn2{i})/(length(T_controln2) - time_active_cn2{i}(1));
end
for i = 1:length(activeM_controln4)
    time_active_cn4{i} = find(M_controln4(:,activeM_controln4(i))>= .8*thresh_base_controln4);
    burst_cn4(i) = length(time_active_cn4{i})/(length(T_controln4) - time_active_cn4{i}(1));
end


for i = 1:length(activeM_bestn1)
    time_active_bn1{i} = find(M_bestn1(:,activeM_bestn1(i))>= .8*thresh_base_dl1bestn1);
    burst_bn1(i) = length(time_active_bn1{i})/(length(T_bestn1) - time_active_bn1{i}(1));
end
for i = 1:length(activeM_bestn2)
    time_active_bn2{i} = find(M_bestn2(:,activeM_bestn2(i))>= .8*thresh_base_dl1bestn2);
    burst_bn2(i) = length(time_active_bn2{i})/(length(T_bestn2) - time_active_bn2{i}(1));
end
for i = 1:length(activeM_bestn4)
    time_active_bn4{i} = find(M_bestn4(:,activeM_bestn4(i))>= .8*thresh_base_dl1bestn4);
    burst_bn4(i) = length(time_active_bn4{i})/(length(T_bestn4) - time_active_bn4{i}(1));
end


for i = 1:length(activeM_deln1)
    time_active_dn1{i} = find(M_deln1(:,activeM_deln1(i))>= .8*thresh_base_dl1deln1);
    burst_dn1(i) = length(time_active_dn1{i})/(length(T_deln1) - time_active_dn1{i}(1));
end
for i = 1:length(activeM_deln2)
    time_active_dn2{i} = find(M_deln2(:,activeM_deln2(i))>= .8*thresh_base_dl1deln2);
    burst_dn2(i) = length(time_active_dn2{i})/(length(T_deln2) - time_active_dn2{i}(1));
end
for i = 1:length(activeM_deln4)
    time_active_dn4{i} = find(M_deln4(:,activeM_deln4(i))>= .8*thresh_base_dl1deln4);
    burst_dn4(i) = length(time_active_dn4{i})/(length(T_deln4) - time_active_dn4{i}(1));
end


for i = 1:length(activeM_dl2bestn1)
    time_active_dl2bn1{i} = find(M_dl2bestn1(:,activeM_dl2bestn1(i))>= .8*thresh_base_dl2bestn1);
    burst_dl2bn1(i) = length(time_active_dl2bn1{i})/(length(T_dl2bestn1) - time_active_dl2bn1{i}(1));
    T_start_dl2b{1}(i) = T_dl2bestn1(time_active_dl2bn1{i}(1))/T_dl2bestn1(end);
    if isinf(burst_dl2bn1(i)) == 1
        burst_dl2bn1(i) = NaN;
    end
end
for i = 1:length(activeM_dl2bestn2)
    time_active_dl2bn2{i} = find(M_dl2bestn2(:,activeM_dl2bestn2(i))>= .8*thresh_base_dl2bestn2);
    burst_dl2bn2(i) = length(time_active_dl2bn2{i})/(length(T_dl2bestn2) - time_active_dl2bn2{i}(1));
    T_start_dl2b{2}(i) = T_dl2bestn2(time_active_dl2bn2{i}(1))/T_dl2bestn2(end);
    if isinf(burst_dl2bn2(i)) == 1
        burst_dl2bn2(i) = NaN;
    end
end
for i = 1:length(activeM_dl2bestn5)
    time_active_dl2bn5{i} = find(M_dl2bestn5(:,activeM_dl2bestn5(i))>= .8*thresh_base_dl2bestn5);
    burst_dl2bn5(i) = length(time_active_dl2bn5{i})/(length(T_dl2bestn5) - time_active_dl2bn5{i}(1));
    T_start_dl2b{3}(i) = T_dl2bestn5(time_active_dl2bn5{i}(1))/T_dl2bestn5(end);
end


for i = 1:length(activeM_dl2deln1)
    time_active_dl2dn1{i} = find(M_dl2deln1(:,activeM_dl2deln1(i))>= .8*thresh_base_dl2deln1);
    burst_dl2dn1(i) = length(time_active_dl2dn1{i})/(length(T_dl2deln1) - time_active_dl2dn1{i}(1));
    T_start_dl2d{1}(i) = T_dl2deln1(time_active_dl2dn1{i}(1))/T_dl2deln1(end);
end
for i = 1:length(activeM_dl2deln2)
    time_active_dl2dn2{i} = find(M_dl2deln2(:,activeM_dl2deln2(i))>= .8*thresh_base_dl2deln2);
    burst_dl2dn2(i) = length(time_active_dl2dn2{i})/(length(T_dl2deln2) - time_active_dl2dn2{i}(1));
    T_start_dl2d{2}(i) = T_dl2deln2(time_active_dl2dn2{i}(1))/T_dl2deln2(end);
    if isinf(burst_dl2dn2(i)) == 1
        burst_dl2dn2(i) = NaN;
    end
end
for i = 1:length(activeM_dl2deln4)
    time_active_dl2dn4{i} = find(M_dl2deln4(:,activeM_dl2deln4(i))>= .8*thresh_base_dl2deln4);
    burst_dl2dn4(i) = length(time_active_dl2dn4{i})/(length(T_dl2deln4) - time_active_dl2dn4{i}(1));
    T_start_dl2d{3}(i) = T_dl2deln4(time_active_dl2dn4{i}(1))/T_dl2deln4(end);
    if isinf(burst_dl2dn4(i)) == 1
        burst_dl2dn4(i) = NaN;
    end
end


for i = 1:length(activeM_zldgoodn2)
    time_active_zldgoodn2{i} = find(M_zldgoodn2(:,activeM_zldgoodn2(i))>= .8*thresh_base_zldgoodn2);
    burst_zgn2(i) = length(time_active_zldgoodn2{i})/(length(T_zldgoodn2) - time_active_zldgoodn2{i}(1));
    T_start_zg{1}(i) = T_zldgoodn2(time_active_zldgoodn2{i}(1))/T_zldgoodn2(end);
end
for i = 1:length(activeM_zldgoodn3)
    time_active_zldgoodn3{i} = find(M_zldgoodn3(:,activeM_zldgoodn3(i))>= .8*thresh_base_zldgoodn3);
    burst_zgn3(i) = length(time_active_zldgoodn3{i})/(length(T_zldgoodn3) - time_active_zldgoodn3{i}(1));
    T_start_zg{2}(i) = T_zldgoodn3(time_active_zldgoodn3{i}(1))/T_zldgoodn3(end);
    if isinf(burst_zgn3(i)) == 1
        burst_zgn3(i) = NaN;
    end
end
for i = 1:length(activeM_zldgoodn4)
    time_active_zldgoodn4{i} = find(M_zldgoodn4(:,activeM_zldgoodn4(i))>= .8*thresh_base_zldgoodn4);
    burst_zgn4(i) = length(time_active_zldgoodn4{i})/(length(T_zldgoodn4) - time_active_zldgoodn4{i}(1));
    T_start_zg{3}(i) = T_zldgoodn4(time_active_zldgoodn4{i}(1))/T_zldgoodn4(end);
end


for i = 1:length(activeM_zldweakn1)
    time_active_zldweakn1{i} = find(M_zldweakn1(:,activeM_zldweakn1(i))>= .8*thresh_base_zldweakn1);
    burst_zwn1(i) = length(time_active_zldweakn1{i})/(length(T_zldweakn1) - time_active_zldweakn1{i}(1));
    T_start_zw{1}(i) = T_zldweakn1(time_active_zldweakn1{i}(1))/T_zldweakn1(end);
end
for i = 1:length(activeM_zldweakn2)
    time_active_zldweakn2{i} = find(M_zldweakn2(:,activeM_zldweakn2(i))>= .8*thresh_base_zldweakn2);
    burst_zwn2(i) = length(time_active_zldweakn2{i})/(length(T_zldweakn2) - time_active_zldweakn2{i}(1));
    T_start_zw{2}(i) = T_zldweakn2(time_active_zldweakn2{i}(1))/T_zldweakn2(end);
end
for i = 1:length(activeM_zldweakn5)
    time_active_zldweakn5{i} = find(M_zldweakn5(:,activeM_zldweakn5(i))>= .8*thresh_base_zldweakn5);
    burst_zwn5(i) = length(time_active_zldweakn5{i})/(length(T_zldweakn5) - time_active_zldweakn5{i}(1));
    T_start_zw{3}(i) = T_zldweakn5(time_active_zldweakn5{i}(1))/T_zldweakn5(end);
end
%% Traces for burst
%{
aa = T_zldgoodn2;
ab = M_dl2bestn2;
ac = loc_dl2bestn2;
ad = activeM_dl2bestn2;
ae = thresh_base_dl2bestn2;

k = find(ac == 0.5);
j = find(ac == 5/16 | ac == 11/16);

figure
for i = 1:10
    subplot(4,5,i)
    plot(aa,ab(:,ad(k(i*6)))); hold on
    plot(aa,.8*ae*ones(size(ab(:,ad(k(i*6)))))); hold off
    ylim([0 5000])
end
for i = 11:20
    subplot(4,5,i)
    plot(aa,ab(:,ad(j((i-10)*10)))); hold on
    plot(aa,.8*ae*ones(size(ab(:,ad(j((i-10)*10)))))); hold off
    ylim([0 5000])
end
% title('Nuclei Signal');
% xlabel('Time (min)');
% ylabel('Pixel Intensity');
% legend('Pixel Intensity','Activation threshold');
% ylim([0 6000])
%}
%% Beginning Activation
clear tmp
locations = [0:15]/16;
for i = 1:length(locations)
    tmp{1,i} = T_start_c{1}(loc_controln1 == locations(i));
%     heat_c(1,i) = mean(tmp{i});
    tmp{2,i} = T_start_c{2}(loc_controln2 == locations(i));
%     heat_c(2,i) = mean(tmp{i});
    tmp{3,i} = T_start_c{3}(loc_controln4 == locations(i));
%     heat_c(3,i) = mean(tmp{i});
    tmp{4,i} = T_start_c{4}(loc_controln3 == locations(i));
    tmp{5,i} = T_start_c{5}(loc_controln5 == locations(i));
    tmp{6,i} = T_start_c{6}(loc_controln6 == locations(i));

    heat_c(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i} tmp{4,i} tmp{5,i} tmp{6,i}]);
    std_heat_c(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i} tmp{4,i} tmp{5,i} tmp{6,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i})+length(tmp{4,i})+length(tmp{5,i})+length(tmp{6,i}));
end
clear tmp
for i = 1:length(locations)
    tmp{1,i} = T_start_b{1}(loc_bestn1 == locations(i));
%     heat_b(1,i) = mean(tmp{i});
    tmp{2,i} = T_start_b{2}(loc_bestn2 == locations(i));
%     heat_b(2,i) = mean(tmp{i});
    tmp{3,i} = T_start_b{3}(loc_bestn4 == locations(i));
%     heat_b(3,i) = mean(tmp{i});
    heat_b(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_heat_b(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
end
for i = 1:length(locations)

    tmp{1,i} = T_start_d{1}(loc_deln1 == locations(i));
%     heat_d(1,i) = mean(tmp{i});
    tmp{2,i} = T_start_d{2}(loc_deln2 == locations(i));
%     heat_d(2,i) = mean(tmp{i});
    tmp{3,i} = T_start_d{3}(loc_deln4 == locations(i));
%     heat_d(3,i) = mean(tmp{i});
    heat_d(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_heat_d(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
end
for i = 1:length(locations)

    tmp{1,i} = T_start_dl2b{1}(loc_dl2bestn1 == locations(i));
%     heat_dl2b(1,i) = mean(tmp{i});
    tmp{2,i} = T_start_dl2b{2}(loc_dl2bestn2 == locations(i));
%     heat_dl2b(2,i) = mean(tmp{i});
    tmp{3,i} = T_start_dl2b{3}(loc_dl2bestn5 == locations(i));
%     heat_dl2b(3,i) = mean(tmp{i});
    heat_dl2b(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_heat_dl2b(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
end
for i = 1:length(locations)

    tmp{1,i} = T_start_dl2d{1}(loc_dl2deln1 == locations(i));
%     heat_dl2d(1,i) = mean(tmp{i});
    tmp{2,i} = T_start_dl2d{2}(loc_dl2deln2 == locations(i));
%     heat_dl2d(2,i) = mean(tmp{i});
    tmp{3,i} = T_start_dl2d{3}(loc_dl2deln4 == locations(i));
%     heat_dl2d(3,i) = mean(tmp{i});
    heat_dl2d(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_heat_dl2d(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
end
for i = 1:length(locations)

    tmp{1,i} = T_start_zg{1}(loc_zldgoodn2 == locations(i));
%     heat_zg(1,i) = mean(tmp{i});
    tmp{2,i} = T_start_zg{2}(loc_zldgoodn3 == locations(i));
%     heat_zg(2,i) = mean(tmp{i});
    tmp{3,i} = T_start_zg{3}(loc_zldgoodn4 == locations(i));
%     heat_zg(3,i) = mean(tmp{i});
    heat_zg(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_heat_zg(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
end
    tmp{1,i} = T_start_td{1}(loc_twideln1 == locations(i));
    tmp{2,i} = T_start_td{2}(loc_twideln2 == locations(i));

    heat_td(i) = nanmean([tmp{1,i} tmp{2,i}]);
    std_heat_td(i) = nanstd([tmp{1,i} tmp{2,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i}));
 
%Revision
clear tmp
locations = [0:15]/16;
for i = 1:length(locations)
    tmp{1,i} = T_start_Rc{1}(loc_Rcontroln1 == locations(i));
%     heat_c(1,i) = mean(tmp{i});
    tmp{2,i} = T_start_Rc{2}(loc_Rcontroln2 == locations(i));
%     heat_c(2,i) = mean(tmp{i});
    tmp{3,i} = T_start_Rc{3}(loc_Rcontroln3 == locations(i));

    heat_Rc(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_heat_Rc(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
end
locations = [0:15]/16;
for i = 1:length(locations)
    tmp{1,i} = T_start_Rdz{1}(loc_Rdzn1 == locations(i));
%     heat_c(1,i) = mean(tmp{i});
    tmp{2,i} = T_start_Rdz{2}(loc_Rdzn2 == locations(i));
%     heat_c(2,i) = mean(tmp{i});
    tmp{3,i} = T_start_Rdz{3}(loc_Rdzn3 == locations(i));

    heat_Rdz(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_heat_Rdz(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
end
locations = [0:15]/16;
for i = 1:length(locations)
    tmp{1,i} = T_start_Rd{1}(loc_Rdn1 == locations(i));
%     heat_c(1,i) = mean(tmp{i});
    tmp{2,i} = T_start_Rd{2}(loc_Rdn2 == locations(i));
%     heat_c(2,i) = mean(tmp{i});
    tmp{3,i} = T_start_Rd{3}(loc_Rdn3 == locations(i));

    heat_Rd(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_heat_Rd(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
end
%     tmp{i} = T_start_zw{1}(loc_zldweakn1 == locations(i));
%     heat_zw(1,i) = mean(tmp{i});
%     tmp{i} = T_start_zw{2}(loc_zldweakn2 == locations(i));
%     heat_zw(2,i) = mean(tmp{i});
%     tmp{i} = T_start_zw{3}(loc_zldweakn5 == locations(i));
%     heat_zw(3,i) = mean(tmp{i});

% % heat_c(isnan(heat_c)) = 1;
% std_heat_c = nanstd(heat_c);
% heat_c = nanmean(heat_c);
% % heat_b(isnan(heat_b)) = 1;
% std_heat_b = nanstd(heat_b);
% heat_b = nanmean(heat_b);
% % heat_d(isnan(heat_d)) = 1;
% std_heat_d = nanstd(heat_d);
% heat_d = nanmean(heat_d);
% % heat_dl2b(isnan(heat_dl2b)) = 1;
% std_heat_dl2b = nanstd(heat_dl2b);
% heat_dl2b = nanmean(heat_dl2b);
% % heat_dl2d(isnan(heat_dl2d)) = 1;
% std_heat_dl2d = nanstd(heat_dl2d);
% heat_dl2d = nanmean(heat_dl2d);
% % heat_zg(isnan(heat_zg)) = 1;
% std_heat_zg = std(heat_zg);
% heat_zg = mean(heat_zg);
% % heat_zw(isnan(heat_zw)) = 1;
% % std_heat_zw = std(heat_zw);
% % heat_zw = mean(heat_zw);

%het
%{
for i = 1:length(locations)
    tmp{i} = T_start_hetc{1}(loc_hetcontroln1 == locations(i));
    heat_hetc(1,i) = mean(tmp{i});
    tmp{i} = T_start_hetc{2}(loc_hetcontroln2 == locations(i));
    heat_hetc(2,i) = mean(tmp{i});
    tmp{i} = T_start_hetc{3}(loc_hetcontroln3 == locations(i));
    heat_hetc(3,i) = mean(tmp{i});
    
    tmp{i} = T_start_hetb{1}(loc_hetbestn1 == locations(i));
    heat_hetb(1,i) = mean(tmp{i});
    tmp{i} = T_start_hetb{2}(loc_hetbestn2 == locations(i));
    heat_hetb(2,i) = mean(tmp{i});
    tmp{i} = T_start_hetb{3}(loc_hetbestn3 == locations(i));
    heat_hetb(3,i) = mean(tmp{i});
       
    tmp{i} = T_start_hetdl2b{1}(loc_hetdl2bestn1 == locations(i));
    heat_hetdl2b(1,i) = mean(tmp{i});
    tmp{i} = T_start_hetdl2b{2}(loc_hetdl2bestn2 == locations(i));
    heat_hetdl2b(2,i) = mean(tmp{i});
    tmp{i} = T_start_hetdl2b{3}(loc_hetdl2bestn3 == locations(i));
    heat_hetdl2b(3,i) = mean(tmp{i});
    
    tmp{i} = T_start_hetdl2d{1}(loc_hetdl2deln1 == locations(i));
    heat_hetdl2d(1,i) = mean(tmp{i});
    tmp{i} = T_start_hetdl2d{2}(loc_hetdl2deln2 == locations(i));
    heat_hetdl2d(2,i) = mean(tmp{i});
    tmp{i} = T_start_hetdl2d{3}(loc_hetdl2deln3 == locations(i));
    heat_hetdl2d(3,i) = mean(tmp{i});
    
    tmp{i} = T_start_hetzg{1}(loc_hetzldbestn1 == locations(i));
    heat_hetzg(1,i) = mean(tmp{i});
    tmp{i} = T_start_hetzg{2}(loc_hetzldbestn2 == locations(i));
    heat_hetzg(2,i) = mean(tmp{i});
    tmp{i} = T_start_hetzg{3}(loc_hetzldbestn3 == locations(i));
    heat_hetzg(3,i) = mean(tmp{i});
    
end
std_heat_hetc = nanstd(heat_hetc);
heat_hetc = nanmean(heat_hetc);
std_heat_hetb = nanstd(heat_hetb);
heat_hetb = nanmean(heat_hetb);
std_heat_hetdl2b = nanstd(heat_hetdl2b);
heat_hetdl2b = nanmean(heat_hetdl2b);
std_heat_hetdl2d = nanstd(heat_hetdl2d);
heat_hetdl2d = nanmean(heat_hetdl2d);
std_heat_hetzg = std(heat_hetzg);
heat_hetzg = mean(heat_hetzg);
%}
%go active

figure,errorbar(locations,heat_c,std_heat_c,'k','LineWidth',2);hold on
errorbar(locations,heat_b,std_heat_b,'color',rgb('red'),'LineWidth',2); hold on
errorbar(locations,heat_d,std_heat_d,'--','color',rgb('lightcoral'),'LineWidth',2)
ylim([.3 1])
xlim([0 .9375])
title('Time of Activation Across Embryo')
ylabel('Activation Time')
xlabel('Location')
legend('control','Dl1 strong','Dl1 weak','Location', 'northwest')

figure,errorbar(locations,heat_c,std_heat_c,'k','LineWidth',2);hold on
errorbar(locations,heat_b,std_heat_b,'ro'); hold on
errorbar(locations,heat_d,std_heat_d,'r*')
ylim([.2 1])

figure,errorbar(locations,heat_c,std_heat_c,'k','LineWidth',2);hold on
errorbar(locations,heat_dl2b,std_heat_dl2b,'color',rgb('blue'),'LineWidth',2)
errorbar(locations,heat_dl2d,std_heat_dl2d,'--','color',rgb('cornflowerblue'),'LineWidth',2)
ylim([.3 1])
xlim([0 .9375])
title('Time of Activation Across Embryo')
ylabel('Activation Time')
xlabel('Location')
legend('control','Dl2 strong','Dl2 weak','Location', 'northwest')

figure,errorbar(locations,heat_b,std_heat_b,'color',rgb('red'),'LineWidth',2); hold on
errorbar(locations,heat_dl2b,std_heat_dl2b,'color',rgb('blue'),'LineWidth',2)
errorbar(locations,heat_zg,std_heat_zg,'color',rgb('green'),'LineWidth',2)
ylim([.2 1])
xlim([0 .9375])
title('Time of Activation Across Embryo')
ylabel('Activation Time')
xlabel('Location')
legend('Dl1 Strong','Dl2 Strong','Zld strong','Location', 'northwest')


%Revision
figure,errorbar(locations,heat_c,std_heat_c,'Color',rgb('gray'),'LineWidth',2);hold on
errorbar(locations,heat_Rc,std_heat_Rc,'Color',rgb('black'),'LineWidth',2);hold on
errorbar(locations,heat_Rdz,std_heat_Rdz,'Color',rgb('red'),'LineWidth',2);hold on
errorbar(locations,heat_Rd,std_heat_Rd,'Color',rgb('blue'),'LineWidth',2);hold on

ylim([.3 1])
xlim([0 .9375])
title('Time of Activation Across Embryo')
ylabel('Activation Time')
xlabel('Location')
legend('control','Dl1 strong','Dl1 weak','Location', 'northwest')

%het
%{
figure,errorbar(locations,heat_hetc,std_heat_hetc,'ko');hold on
errorbar(locations,heat_hetb,std_heat_hetb,'o','color',rgb('red')); hold on
ylim([.3 1])
xlim([0 .9375])
title('Time of Activation Across Embryo')
ylabel('Activation Time')
xlabel('Location')
legend('het control','het Dl1 strong','Location', 'northwest')

figure,errorbar(locations,heat_hetc,std_heat_hetc,'ko');hold on
errorbar(locations,heat_hetdl2b,std_heat_hetdl2b,'o','color',rgb('blue'))
errorbar(locations,heat_hetdl2d,std_heat_hetdl2d,'o','color',rgb('cornflowerblue'))
ylim([.3 1])
xlim([0 .9375])
title('Time of Activation Across Embryo')
ylabel('Activation Time')
xlabel('Location')
legend('het control','het Dl2 strong','het Dl2 weak','Location', 'northwest')

figure,errorbar(locations,heat_hetc,std_heat_hetc,'ko');hold on
errorbar(locations,heat_hetzg,std_heat_hetzg,'o','color',rgb('green')); hold on
ylim([.3 1])
xlim([0 .9375])
title('Time of Activation Across Embryo')
ylabel('Activation Time')
xlabel('Location')
legend('het control','het Zld strong','Location', 'northwest')
%}
%% Max Active Signal

%Max signal reached by active nuclei
Max_c{1} = Max_c{1}(activeM_controln1);
Max_c{2} = Max_c{2}(activeM_controln2);
Max_c{3} = Max_c{3}(activeM_controln4);
Max_c{4} = Max_c{4}(activeM_controln3);
Max_c{5} = Max_c{5}(activeM_controln5);
Max_c{6} = Max_c{6}(activeM_controln6);
Max_b{1} = Max_b{1}(activeM_bestn1);
Max_b{2} = Max_b{2}(activeM_bestn2);
Max_b{3} = Max_b{3}(activeM_bestn4);
Max_d{1} = Max_d{1}(activeM_deln1);
Max_d{2} = Max_d{2}(activeM_deln2);
Max_d{3} = Max_d{3}(activeM_deln4);
Max_dl2b{1} = Max_dl2b{1}(activeM_dl2bestn1);
Max_dl2b{2} = Max_dl2b{2}(activeM_dl2bestn2);
Max_dl2b{3} = Max_dl2b{3}(activeM_dl2bestn5);
Max_dl2d{1} = Max_dl2d{1}(activeM_dl2deln1);
Max_dl2d{2} = Max_dl2d{2}(activeM_dl2deln2);
Max_dl2d{3} = Max_dl2d{3}(activeM_dl2deln4);
Max_zg{1} = Max_zg{1}(activeM_zldgoodn2);
Max_zg{2} = Max_zg{2}(activeM_zldgoodn3);
Max_zg{3} = Max_zg{3}(activeM_zldgoodn4);
Max_zw{1} = Max_zw{1}(activeM_zldweakn1);
Max_zw{2} = Max_zw{2}(activeM_zldweakn2);
Max_zw{3} = Max_zw{3}(activeM_zldweakn5);
Max_td{1} = Max_td{1}(activeM_twideln1);
Max_td{2} = Max_td{2}(activeM_twideln2);

%Revision
Max_Rc{1} = Max_Rc{1}(activeM_Rcn1);
Max_Rc{2} = Max_Rc{2}(activeM_Rcn2);
Max_Rc{3} = Max_Rc{3}(activeM_Rcn3);
Max_Rdz{1} = Max_Rdz{1}(activeM_Rdzn1);
Max_Rdz{2} = Max_Rdz{2}(activeM_Rdzn2);
Max_Rdz{3} = Max_Rdz{3}(activeM_Rdzn3);
Max_Rd{1} = Max_Rd{1}(activeM_Rdn1);
Max_Rd{2} = Max_Rd{2}(activeM_Rdn2);
Max_Rd{3} = Max_Rd{3}(activeM_Rdn3);


clear tmp
locations = [0:15]/16;
for i = 1:length(locations)
    tmp{1,i} = Max_c{1}(loc_controln1 == locations(i));
%     max_c(1,i) = mean(tmp{i});
    tmp{2,i} = Max_c{2}(loc_controln2 == locations(i));
%     max_c(2,i) = mean(tmp{i});
    tmp{3,i} = Max_c{3}(loc_controln4 == locations(i));
%     max_c(3,i) = mean(tmp{i});
    tmp{4,i} = Max_c{4}(loc_controln3 == locations(i));
    tmp{5,i} = Max_c{5}(loc_controln5 == locations(i));
    tmp{6,i} = Max_c{6}(loc_controln6 == locations(i));

    max_c(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i} tmp{4,i} tmp{5,i} tmp{6,i}]);
    std_max_c(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i} tmp{4,i} tmp{5,i} tmp{6,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i})+length(tmp{4,i})+length(tmp{5,i})+length(tmp{6,i}));
    clear tmp
    tmp{1,i} = Max_b{1}(loc_bestn1 == locations(i));
%     max_b(1,i) = mean(tmp{i});
    tmp{2,i} = Max_b{2}(loc_bestn2 == locations(i));
%     max_b(2,i) = mean(tmp{i});
    tmp{3,i} = Max_b{3}(loc_bestn4 == locations(i));
%     max_b(3,i) = mean(tmp{i});
    max_b(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_max_b(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
    clear tmp
    tmp{1,i} = Max_d{1}(loc_deln1 == locations(i));
%     max_d(1,i) = mean(tmp{i});
    tmp{2,i} = Max_d{2}(loc_deln2 == locations(i));
%     max_d(2,i) = mean(tmp{i});
    tmp{3,i} = Max_d{3}(loc_deln4 == locations(i));
%     max_d(3,i) = mean(tmp{i});
    max_d(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_max_d(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
    clear tmp
    tmp{1,i} = Max_dl2b{1}(loc_dl2bestn1 == locations(i));
%     max_dl2b(1,i) = mean(tmp{i});
    tmp{2,i} = Max_dl2b{2}(loc_dl2bestn2 == locations(i));
%     max_dl2b(2,i) = mean(tmp{i});
    tmp{3,i} = Max_dl2b{3}(loc_dl2bestn5 == locations(i));
%     max_dl2b(3,i) = mean(tmp{i});
    max_dl2b(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_max_dl2b(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
    clear tmp
    tmp{1,i} = Max_dl2d{1}(loc_dl2deln1 == locations(i));
%     max_dl2d(1,i) = mean(tmp{i});
    tmp{2,i} = Max_dl2d{2}(loc_dl2deln2 == locations(i));
%     max_dl2d(2,i) = mean(tmp{i});
    tmp{3,i} = Max_dl2d{3}(loc_dl2deln4 == locations(i));
%     max_dl2d(3,i) = mean(tmp{i});
    max_dl2d(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_max_dl2d(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
    clear tmp
    tmp{1,i} = Max_zg{1}(loc_zldgoodn2 == locations(i));
%     max_zg(1,i) = mean(tmp{i});
    tmp{2,i} = Max_zg{2}(loc_zldgoodn3 == locations(i));
%     max_zg(2,i) = mean(tmp{i});
    tmp{3,i} = Max_zg{3}(loc_zldgoodn4 == locations(i));
%     max_zg(3,i) = mean(tmp{i});
    max_zg(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_max_zg(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
    clear tmp
    tmp{1,i} = Max_td{1}(loc_twideln1 == locations(i));
    tmp{2,i} = Max_td{2}(loc_twideln2 == locations(i));

    max_td(i) =  nanmean([tmp{1,i} tmp{2,i}]);
    std_max_td(i) = nanstd([tmp{1,i} tmp{2,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i}));
    clear tmp
    %Revision
    tmp{1,i} = Max_Rc{1}(loc_Rcontroln1 == locations(i));
%     max_c(1,i) = mean(tmp{i});
    tmp{2,i} = Max_Rc{2}(loc_Rcontroln2 == locations(i));
%     max_c(2,i) = mean(tmp{i});
    tmp{3,i} = Max_Rc{3}(loc_Rcontroln3 == locations(i));

    max_Rc(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_max_Rc(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
    clear tmp
    tmp{1,i} = Max_Rdz{1}(loc_Rdzn1 == locations(i));
    tmp{2,i} = Max_Rdz{2}(loc_Rdzn2 == locations(i));
    tmp{3,i} = Max_Rdz{3}(loc_Rdzn3 == locations(i));

    max_Rdz(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_max_Rdz(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
    clear tmp
    tmp{1,i} = Max_Rd{1}(loc_Rdn1 == locations(i));
    tmp{2,i} = Max_Rd{2}(loc_Rdn2 == locations(i));
    tmp{3,i} = Max_Rd{3}(loc_Rdn3 == locations(i));

    max_Rd(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_max_Rd(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
    clear tmp
  
    
end

% std_max_c = nanstd(max_c);
% max_c = nanmean(max_c);
% std_max_b = nanstd(max_b);
% max_b = nanmean(max_b);
% std_max_d = nanstd(max_d);
% max_d = nanmean(max_d);
% std_max_dl2b = nanstd(max_dl2b);
% max_dl2b = nanmean(max_dl2b);
% std_max_dl2d = nanstd(max_dl2d);
% max_dl2d = nanmean(max_dl2d);
% std_max_zg = nanstd(max_zg);
% max_zg = nanmean(max_zg);
% 

%het
Max_hetc{1} = Max_hetc{1}(activeM_hetcontroln1);
Max_hetc{2} = Max_hetc{2}(activeM_hetcontroln2);
Max_hetc{3} = Max_hetc{3}(activeM_hetcontroln3);
Max_hetb{1} = Max_hetb{1}(activeM_hetbestn1);
Max_hetb{2} = Max_hetb{2}(activeM_hetbestn2);
Max_hetb{3} = Max_hetb{3}(activeM_hetbestn3);
Max_hetdl2b{1} = Max_hetdl2b{1}(activeM_hetdl2bestn1);
Max_hetdl2b{2} = Max_hetdl2b{2}(activeM_hetdl2bestn2);
Max_hetdl2b{3} = Max_hetdl2b{3}(activeM_hetdl2bestn3);
Max_hetdl2d{1} = Max_hetdl2d{1}(activeM_hetdl2deln1);
Max_hetdl2d{2} = Max_hetdl2d{2}(activeM_hetdl2deln2);
Max_hetdl2d{3} = Max_hetdl2d{3}(activeM_hetdl2deln3);
Max_hetzg{1} = Max_hetzg{1}(activeM_hetzldbestn1);
Max_hetzg{2} = Max_hetzg{2}(activeM_hetzldbestn2);
Max_hetzg{3} = Max_hetzg{3}(activeM_hetzldbestn3);

clear tmp
%{
locations = [0:15]/16;
for i = 1:length(locations)
    tmp{i} = Max_hetc{1}(loc_hetcontroln1 == locations(i));
    max_hetc(1,i) = mean(tmp{i});
    tmp{i} = Max_hetc{2}(loc_hetcontroln2 == locations(i));
    max_hetc(2,i) = mean(tmp{i});
    tmp{i} = Max_hetc{3}(loc_hetcontroln3 == locations(i));
    max_hetc(3,i) = mean(tmp{i});
    
    tmp{i} = Max_hetb{1}(loc_hetbestn1 == locations(i));
    max_hetb(1,i) = mean(tmp{i});
    tmp{i} = Max_hetb{2}(loc_hetbestn2 == locations(i));
    max_hetb(2,i) = mean(tmp{i});
    tmp{i} = Max_hetb{3}(loc_hetbestn3 == locations(i));
    max_hetb(3,i) = mean(tmp{i});
        
    tmp{i} = Max_hetdl2b{1}(loc_hetdl2bestn1 == locations(i));
    max_hetdl2b(1,i) = mean(tmp{i});
    tmp{i} = Max_hetdl2b{2}(loc_hetdl2bestn2 == locations(i));
    max_hetdl2b(2,i) = mean(tmp{i});
    tmp{i} = Max_hetdl2b{3}(loc_hetdl2bestn3 == locations(i));
    max_hetdl2b(3,i) = mean(tmp{i});
    
    tmp{i} = Max_hetdl2d{1}(loc_hetdl2deln1 == locations(i));
    max_hetdl2d(1,i) = mean(tmp{i});
    tmp{i} = Max_hetdl2d{2}(loc_hetdl2deln2 == locations(i));
    max_hetdl2d(2,i) = mean(tmp{i});
    tmp{i} = Max_hetdl2d{3}(loc_hetdl2deln3 == locations(i));
    max_hetdl2d(3,i) = mean(tmp{i});
    
    tmp{i} = Max_hetzg{1}(loc_hetzldbestn1 == locations(i));
    max_hetzg(1,i) = mean(tmp{i});
    tmp{i} = Max_hetzg{2}(loc_hetzldbestn2 == locations(i));
    max_hetzg(2,i) = mean(tmp{i});
    tmp{i} = Max_hetzg{3}(loc_hetzldbestn3 == locations(i));
    max_hetzg(3,i) = mean(tmp{i});
    
end

std_max_hetc = nanstd(max_hetc);
max_hetc = nanmean(max_hetc);
std_max_hetb = nanstd(max_hetb);
max_hetb = nanmean(max_hetb);
std_max_hetdl2b = nanstd(max_hetdl2b);
max_hetdl2b = nanmean(max_hetdl2b);
std_max_hetdl2d = nanstd(max_hetdl2d);
max_hetdl2d = nanmean(max_hetdl2d);
std_max_hetzg = nanstd(max_hetzg);
max_hetzg = nanmean(max_hetzg);
%}
%go max

figure,errorbar(locations,max_c,std_max_c,'k','LineWidth',2);hold on
errorbar(locations,max_b,std_max_b,'color',rgb('red'),'LineWidth',2);hold on
errorbar(locations,max_d,std_max_d,'--','color',rgb('lightcoral'),'LineWidth',2)
ylim([2000 5000])
xlim([0 .9375])
title('Maximum Signal Intensity')
ylabel('Max Intensity')
xlabel('location')
legend('control','Dl1 strong','Dl1 weak','Location','northwest')

figure,errorbar(locations,max_c,std_max_c,'k','LineWidth',2);hold on
errorbar(locations,max_dl2b,std_max_dl2b,'b','LineWidth',2)
errorbar(locations,max_dl2d,std_max_dl2d,'--','color',rgb('cornflowerblue'),'LineWidth',2)
ylim([2000 5000])
xlim([0 .9375])
title('Maximum Signal Intensity')
ylabel('Max Intensity')
xlabel('location')
legend('control','Dl2 strong','Dl2 weak','Location','northwest')

figure,errorbar(locations,max_c,std_max_c,'k','LineWidth',2);hold on
errorbar(locations,max_b,std_max_b,'color',rgb('red'),'LineWidth',2);hold on
errorbar(locations,max_dl2b,std_max_dl2b,'b','LineWidth',2)
errorbar(locations,max_zg,std_max_zg,'color',rgb('green'),'LineWidth',2)
ylim([2000 5000])
xlim([0 .9375])
title('Maximum Signal Intensity')
ylabel('Max Intensity')
xlabel('location')
legend('control','Dl1 strong','Dl2 strong','Zld strong','Location','northwest')

%het
%{
figure,errorbar(locations,max_hetc,std_max_hetc,'k');hold on
errorbar(locations,max_hetb,std_max_hetb,'color',rgb('red'));hold on
ylim([0 8000])
xlim([0 .9375])
title('Maximum Signal Intensity')
ylabel('Max Intensity')
xlabel('location')
legend('het control','het Dl1 strong','Location','northwest')

figure,errorbar(locations,max_hetc,std_max_hetc,'k');hold on
errorbar(locations,max_hetdl2b,std_max_hetdl2b,'b')
errorbar(locations,max_hetdl2d,std_max_hetdl2d,'--','color',rgb('cornflowerblue'))
ylim([0 8000])
xlim([0 .9375])
title('Maximum Signal Intensity')
ylabel('Max Intensity')
xlabel('location')
legend('het control','het Dl2 strong','hetDl2 weak','Location','northwest')

figure,errorbar(locations,max_hetc,std_max_hetc,'k');hold on
errorbar(locations,max_hetzg,std_max_hetzg,'color',rgb('green'));hold on
ylim([0 8000])
xlim([0 .9375])
title('Maximum Signal Intensity')
ylabel('Max Intensity')
xlabel('location')
legend('het control','het Zld strong','Location','northwest')

figure,errorbar(locations,max_hetb,std_max_hetb,'color',rgb('red'));hold on
errorbar(locations,max_hetdl2b,std_max_hetdl2b,'b')
errorbar(locations,max_hetzg,std_max_hetzg,'color',rgb('green'))
ylim([0 8000])
xlim([0 .9375])
title('Maximum Signal Intensity')
ylabel('Max Intensity')
xlabel('location')
legend('het Dl1 strong','het Dl2 strong','het Zld strong','Location','northwest')
%}

%% Time to Max signal
%Control
for i = 1:16
    for j = 1:length(intersect(shift_bin_controln1{i},activeM_controln1))
        tmp = intersect(shift_bin_controln1{i},activeM_controln1);
        [~,maxtimec1{i}] = max(M_controln1(:,tmp));
        Tmaxc{1,i}(j) = T_controln1(maxtimec1{i}(j))/T_controln1(end);
    end
end
if length(Tmaxc) ~= 16
    for i = length(Tmaxc)+1:16
        Tmaxc{1,i} = nan;
    end
end

for i = 1:16
    for j = 1:length(intersect(shift_bin_controln2{i},activeM_controln2))
        tmp = intersect(shift_bin_controln2{i},activeM_controln2);
        [~,maxtimec2{i}] = max(M_controln2(:,tmp));
        Tmaxc{2,i}(j) = T_controln2(maxtimec2{i}(j))/T_controln2(end);
    end
end


for i = 1:16
    for j = 1:length(intersect(shift_bin_controln4{i},activeM_controln4))
        tmp = intersect(shift_bin_controln4{i},activeM_controln4);
        [~,maxtimec4{i}] = max(M_controln4(:,tmp));
        Tmaxc{3,i}(j) = T_controln4(maxtimec4{i}(j))/T_controln4(end);
    end
end

for i = 1:16
    for j = 1:length(intersect(shift_bin_controln3{i},activeM_controln3))
        tmp = intersect(shift_bin_controln3{i},activeM_controln3);
        [~,maxtimec3{i}] = max(M_controln3(:,tmp));
        Tmaxc{4,i}(j) = T_controln3(maxtimec3{i}(j))/T_controln3(end);
    end
end

for i = 1:16
    for j = 1:length(intersect(shift_bin_controln5{i},activeM_controln5))
        tmp = intersect(shift_bin_controln5{i},activeM_controln5);
        [~,maxtimec5{i}] = max(M_controln5(:,tmp));
        Tmaxc{5,i}(j) = T_controln5(maxtimec5{i}(j))/T_controln5(end);
    end
end

for i = 1:16
    for j = 1:length(intersect(shift_bin_controln6{i},activeM_controln6))
        tmp = intersect(shift_bin_controln6{i},activeM_controln6);
        [~,maxtimec6{i}] = max(M_controln6(:,tmp));
        Tmaxc{6,i}(j) = T_controln6(maxtimec6{i}(j))/T_controln6(end);
    end
end

%Best
for i = 1:16
    for j = 1:length(intersect(shift_bin_bestn1{i},activeM_bestn1))
        tmp = intersect(shift_bin_bestn1{i},activeM_bestn1);
        [~,maxtimeb1{i}] = max(M_bestn1(:,tmp));
        Tmaxb{1,i}(j) = T_bestn1(maxtimeb1{i}(j))/T_bestn1(end);
    end
end
if length(Tmaxb) ~= 16
    for i = length(Tmaxb)+1:16
        Tmaxb{1,i} = nan;
    end
end

for i = 1:16
    for j = 1:length(intersect(shift_bin_bestn2{i},activeM_bestn2))
        tmp = intersect(shift_bin_bestn2{i},activeM_bestn2);
        [~,maxtimeb2{i}] = max(M_bestn2(:,tmp));
        Tmaxb{2,i}(j) = T_bestn2(maxtimeb2{i}(j))/T_bestn2(end);
    end
end

for i = 1:16
    for j = 1:length(intersect(shift_bin_bestn4{i},activeM_bestn4))
        tmp = intersect(shift_bin_bestn4{i},activeM_bestn4);
        [~,maxtimeb4{i}] = max(M_bestn4(:,tmp));
        Tmaxb{3,i}(j) = T_bestn4(maxtimeb4{i}(j))/T_bestn4(end);
    end
end

%Del
for i = 1:16
    for j = 1:length(intersect(shift_bin_deln1{i},activeM_deln1))
        tmp = intersect(shift_bin_deln1{i},activeM_deln1);
        [~,maxtimed1{i}] = max(M_deln1(:,tmp));
        Tmaxd{1,i}(j) = T_deln1(maxtimed1{i}(j))/T_deln1(end);
    end
end
if length(Tmaxd) ~= 16
    for i = length(Tmaxd)+1:16
        Tmaxd{1,i} = nan;
    end
end

for i = 1:16
    for j = 1:length(intersect(shift_bin_deln2{i},activeM_deln2))
        tmp = intersect(shift_bin_deln2{i},activeM_deln2);
        [~,maxtimed2{i}] = max(M_deln2(:,tmp));
        Tmaxd{2,i}(j) = T_deln2(maxtimed2{i}(j))/T_deln2(end);
    end
end

for i = 1:16
    for j = 1:length(intersect(shift_bin_deln4{i},activeM_deln4))
        tmp = intersect(shift_bin_deln4{i},activeM_deln4);
        [~,maxtimed4{i}] = max(M_deln4(:,tmp));
        Tmaxd{3,i}(j) = T_deln4(maxtimed4{i}(j))/T_deln4(end);
    end
end

%Dl2 Best
for i = 1:16
    for j = 1:length(intersect(shift_bin_dl2bestn1{i},activeM_dl2bestn1))
        tmp = intersect(shift_bin_dl2bestn1{i},activeM_dl2bestn1);
        [~,maxtimedl2b1{i}] = max(M_dl2bestn1(:,tmp));
        Tmaxdl2b{1,i}(j) = T_dl2bestn1(maxtimedl2b1{i}(j))/T_dl2bestn1(end);
    end
end
if length(Tmaxdl2b) ~= 16
    for i = length(Tmaxdl2b)+1:16
        Tmaxdl2b{1,i} = nan;
    end
end

for i = 1:16
    for j = 1:length(intersect(shift_bin_dl2bestn2{i},activeM_dl2bestn2))
        tmp = intersect(shift_bin_dl2bestn2{i},activeM_dl2bestn2);
        [~,maxtimedl2b2{i}] = max(M_dl2bestn2(:,tmp));
        Tmaxdl2b{2,i}(j) = T_dl2bestn2(maxtimedl2b2{i}(j))/T_dl2bestn2(end);
    end
end

for i = 1:16
    for j = 1:length(intersect(shift_bin_dl2bestn5{i},activeM_dl2bestn5))
        tmp = intersect(shift_bin_dl2bestn5{i},activeM_dl2bestn5);
        [~,maxtimedl2b5{i}] = max(M_dl2bestn5(:,tmp));
        Tmaxdl2b{3,i}(j) = T_dl2bestn5(maxtimedl2b5{i}(j))/T_dl2bestn5(end);
    end
end

%Dl2 Del
for i = 1:16
    for j = 1:length(intersect(shift_bin_dl2deln1{i},activeM_dl2deln1))
        tmp = intersect(shift_bin_dl2deln1{i},activeM_dl2deln1);
        [~,maxtimedl2d1{i}] = max(M_dl2deln1(:,tmp));
        Tmaxdl2d{1,i}(j) = T_dl2deln1(maxtimedl2d1{i}(j))/T_dl2deln1(end);
    end
end
if length(Tmaxdl2d) ~= 16
    for i = length(Tmaxdl2d)+1:16
        Tmaxdl2d{1,i} = nan;
    end
end

for i = 1:16
    for j = 1:length(intersect(shift_bin_dl2deln2{i},activeM_dl2deln2))
        tmp = intersect(shift_bin_dl2deln2{i},activeM_dl2deln2);
        [~,maxtimedl2d2{i}] = max(M_dl2deln2(:,tmp));
        Tmaxdl2d{2,i}(j) = T_dl2deln2(maxtimedl2d2{i}(j))/T_dl2deln2(end);
    end
end

for i = 1:16
    for j = 1:length(intersect(shift_bin_dl2deln4{i},activeM_dl2deln4))
        tmp = intersect(shift_bin_dl2deln4{i},activeM_dl2deln4);
        [~,maxtimedl2d4{i}] = max(M_dl2deln4(:,tmp));
        Tmaxdl2d{3,i}(j) = T_dl2deln4(maxtimedl2d4{i}(j))/T_dl2deln4(end);
    end
end


%Zld good
for i = 1:16
    for j = 1:length(intersect(shift_bin_zldgoodn2{i},activeM_zldgoodn2))
        tmp = intersect(shift_bin_zldgoodn2{i},activeM_zldgoodn2);
        [~,maxtimezg1{i}] = max(M_zldgoodn2(:,tmp));
        Tmaxzg{1,i}(j) = T_zldgoodn2(maxtimezg1{i}(j))/T_zldgoodn2(end);
    end
end
if length(Tmaxzg) ~= 16
    for i = length(Tmaxzg)+1:16
        Tmaxzg{1,i} = nan;
    end
end

for i = 1:16
    for j = 1:length(intersect(shift_bin_zldgoodn3{i},activeM_zldgoodn3))
        tmp = intersect(shift_bin_zldgoodn3{i},activeM_zldgoodn3);
        [~,maxtimezg{i}] = max(M_zldgoodn3(:,tmp));
        Tmaxzg{2,i}(j) = T_zldgoodn3(maxtimezg{i}(j))/T_zldgoodn3(end);
    end
end

for i = 1:16
    for j = 1:length(intersect(shift_bin_zldgoodn4{i},activeM_zldgoodn4))
        tmp = intersect(shift_bin_zldgoodn4{i},activeM_zldgoodn4);
        [~,maxtimezg4{i}] = max(M_zldgoodn4(:,tmp));
        Tmaxzg{3,i}(j) = T_zldgoodn4(maxtimezg4{i}(j))/T_zldgoodn4(end);
    end
end


for i = 1:16
    meanTc(i) = nanmean([Tmaxc{1,i} Tmaxc{2,i} Tmaxc{3,i} Tmaxc{4,i} Tmaxc{5,i} Tmaxc{5,i}]);
    stdTc(i) = nanstd([Tmaxc{1,i} Tmaxc{2,i} Tmaxc{3,i} Tmaxc{4,i} Tmaxc{5,i} Tmaxc{5,i}])./sqrt(length(Tmaxc{1,i})+length(Tmaxc{2,i})+length(Tmaxc{3,i})+length(Tmaxc{4,i})+length(Tmaxc{5,i})+length(Tmaxc{6,i}));
    meanTb(i) = nanmean([Tmaxb{1,i} Tmaxb{2,i} Tmaxb{3,i}]);
    stdTb(i) = nanstd([Tmaxb{1,i} Tmaxb{2,i} Tmaxb{3,i}])./sqrt(length(Tmaxb{1,i})+length(Tmaxb{2,i})+length(Tmaxb{3,i}));
    meanTd(i) = nanmean([Tmaxd{1,i} Tmaxd{2,i} Tmaxd{3,i}]);
    stdTd(i) = nanstd([Tmaxd{1,i} Tmaxd{2,i} Tmaxd{3,i}])./sqrt(length(Tmaxd{1,i})+length(Tmaxd{2,i})+length(Tmaxd{3,i}));
    meanTdl2b(i) = nanmean([Tmaxdl2b{1,i} Tmaxdl2b{2,i} Tmaxdl2b{3,i}]);
    stdTdl2b(i) = nanstd([Tmaxdl2b{1,i} Tmaxdl2b{2,i} Tmaxdl2b{3,i}])./sqrt(length(Tmaxdl2b{1,i})+length(Tmaxdl2b{2,i})+length(Tmaxdl2b{3,i}));
    meanTdl2d(i) = nanmean([Tmaxdl2d{1,i} Tmaxdl2d{2,i} Tmaxdl2d{3,i}]);
    stdTdl2d(i) = nanstd([Tmaxdl2d{1,i} Tmaxdl2d{2,i} Tmaxdl2d{3,i}])./sqrt(length(Tmaxdl2d{1,i})+length(Tmaxdl2d{2,i})+length(Tmaxdl2d{3,i}));
    meanTzg(i) = nanmean([Tmaxzg{1,i} Tmaxzg{2,i} Tmaxzg{3,i}]);
    stdTzg(i) = nanstd([Tmaxzg{1,i} Tmaxzg{2,i} Tmaxzg{3,i}])./sqrt(length(Tmaxzg{1,i})+length(Tmaxzg{2,i})+length(Tmaxzg{3,i}));

end
meanTc(meanTc == 0) = nan;
stdTc(stdTc == 0) = nan;
meanTb(meanTb == 0) = nan;
stdTb(stdTb == 0) = nan;
meanTd(meanTd == 0) = nan;
stdTd(stdTd == 0) = nan;
meanTdl2b(meanTdl2b == 0) = nan;
stdTdl2b(stdTdl2b == 0) = nan;
meanTdl2d(meanTdl2d == 0) = nan;
stdTdl2d(stdTdl2d == 0) = nan;
meanTzg(meanTzg == 0) = nan;
stdTzg(stdTzg == 0) = nan;

figure;
errorbar([0:15]/16,meanTc,stdTc,'o','color', rgb('black')); hold on
errorbar([0:15]/16,meanTb,stdTb,'o','color', rgb('red')); 
errorbar([0:15]/16,meanTd,stdTd,'o','color', rgb('lightcoral')); hold off
title('Time to Maximum Signal')
ylabel('% Through NC14')
xlabel('location')
legend('control','Dl1 strong','Dl1 weak','Location','northwest')
ylim([.5 1])
xlim([0 .9375])

figure;
errorbar([0:15]/16,meanTc,stdTc,'o','color', rgb('black')); hold on
errorbar([0:15]/16,meanTdl2b,stdTdl2b,'o','color', rgb('blue')); 
errorbar([0:15]/16,meanTdl2d,stdTdl2d,'o','color', rgb('cornflowerblue')); hold off
title('Time to Maximum Signal')
ylabel('% Through NC14')
xlabel('location')
legend('control','Dl2 strong','Dl2 weak','Location','northwest')
ylim([.5 1])
xlim([0 .9375])

figure;
errorbar([0:15]/16,meanTb,stdTb,'o','color', rgb('red')); hold on
errorbar([0:15]/16,meanTdl2b,stdTdl2b,'o','color', rgb('blue')); 
errorbar([0:15]/16,meanTzg,stdTzg,'o','color', rgb('green')); hold off
title('Time to Maximum Signal')
ylabel('% Through NC14')
xlabel('location')
legend('Dl1 strong','Dl2 strong','zld good','Location','northwest')
ylim([.5 1])
xlim([0 .9375])

%% Variance plots
%{
for i = 1:16
        c(1,i) = length(bin_controln1{1,i});
        c(2,i) = length(bin_controln2{1,i});
        c(3,i) = length(bin_controln4{1,i});
        b(1,i) = length(bin_bestn1{1,i});
        b(2,i) = length(bin_bestn2{1,i});
        b(3,i) = length(bin_bestn4{1,i});
        d(1,i) = length(bin_deln1{1,i});
        d(2,i) = length(bin_deln2{1,i});
        d(3,i) = length(bin_deln4{1,i});
        dl2b(1,i) = length(bin_dl2bestn1{1,i});
        dl2b(2,i) = length(bin_dl2bestn2{1,i});
        dl2b(3,i) = length(bin_dl2bestn5{1,i});
        dl2d(1,i) = length(bin_dl2deln1{1,i});
        dl2d(2,i) = length(bin_dl2deln2{1,i});
        dl2d(3,i) = length(bin_dl2deln4{1,i});
end
stdArea_bin_controln1 = sqrt(c(1,:)).*stdArea_bin_controln1;
stdArea_bin_controln2 = sqrt(c(1,:)).*stdArea_bin_controln2;
stdArea_bin_controln4 = sqrt(c(1,:)).*stdArea_bin_controln4;
stdArea_bin_bestn1 = sqrt(b(1,:)).*stdArea_bin_bestn1;
stdArea_bin_bestn2 = sqrt(b(1,:)).*stdArea_bin_bestn2;
stdArea_bin_bestn4 = sqrt(b(1,:)).*stdArea_bin_bestn4;
stdArea_bin_deln1 = sqrt(d(1,:)).*stdArea_bin_deln1;
stdArea_bin_deln2 = sqrt(d(1,:)).*stdArea_bin_deln2;
stdArea_bin_deln4 = sqrt(d(1,:)).*stdArea_bin_deln4;
stdArea_bin_dl2bestn1 = sqrt(dl2b(1,:)).*stdArea_bin_dl2bestn1;
stdArea_bin_dl2bestn2 = sqrt(dl2b(1,:)).*stdArea_bin_dl2bestn2;
stdArea_bin_dl2bestn5 = sqrt(dl2b(1,:)).*stdArea_bin_dl2bestn5;
stdArea_bin_dl2deln1 = sqrt(dl2d(1,:)).*stdArea_bin_dl2deln1;
stdArea_bin_dl2deln2 = sqrt(dl2d(1,:)).*stdArea_bin_dl2deln2;
stdArea_bin_dl2deln4 = sqrt(dl2d(1,:)).*stdArea_bin_dl2deln4;
%}
%% Intensity distribution
thing = [.05:.05:1];
p = .7;
% for q = 1:20
%     p = p + .05;
t = round(p*length(T_controln1));
for i = 1:16
    if isempty(shift_bin_controln1)
        sigc{1,i} = 0;
    else
        for j = 1:length(shift_bin_controln1{i})
            sigc{1,i}(j) = M_controln1(t,shift_bin_controln1{i}(j));
        end
    end
end
t = round(p*length(T_controln2));
for i = 1:16
    if isempty(shift_bin_controln2)
        sigc{2,i} = 0;
    else
        for j = 1:length(shift_bin_controln2{i})
            sigc{2,i}(j) = M_controln2(t,shift_bin_controln2{i}(j));
        end
    end
end
t = round(p*length(T_controln4));
for i = 1:16
    if isempty(shift_bin_controln4)
        sigc{3,i} = 0;
    else
        for j = 1:length(shift_bin_controln4{i})
            sigc{3,i}(j) = M_controln4(t,shift_bin_controln4{i}(j));
        end
    end
end
t = round(p*length(T_controln3));
for i = 1:16
    if isempty(shift_bin_controln3)
        sigc{3,i} = 0;
    else
        for j = 1:length(shift_bin_controln3{i})
            sigc{3,i}(j) = M_controln3(t,shift_bin_controln3{i}(j));
        end
    end
end
t = round(p*length(T_controln5));
for i = 1:16
    if isempty(shift_bin_controln5)
        sigc{5,i} = 0;
    else
        for j = 1:length(shift_bin_controln5{i})
            sigc{5,i}(j) = M_controln5(t,shift_bin_controln5{i}(j));
        end
    end
end
t = round(p*length(T_controln6));
for i = 1:16
    if isempty(shift_bin_controln6)
        sigc{6,i} = 0;
    else
        for j = 1:length(shift_bin_controln6{i})
            sigc{6,i}(j) = M_controln6(t,shift_bin_controln6{i}(j));
        end
    end
end
%Dl1 best
t = round(p*length(T_bestn1));
for i = 1:16
    if isempty(shift_bin_bestn1)
        sigb{1,i} = 0;
    else
        for j = 1:length(shift_bin_bestn1{i})
            sigb{1,i}(j) = M_bestn1(t,shift_bin_bestn1{i}(j));
        end
    end
end
t = round(p*length(T_bestn2));
for i = 1:16
    if isempty(shift_bin_bestn2)
        sigb{2,i} = 0;
    else
        for j = 1:length(shift_bin_bestn2{i})
            sigb{2,i}(j) = M_bestn2(t,shift_bin_bestn2{i}(j));
        end
    end
end
t = round(p*length(T_bestn4));
for i = 1:16
    if isempty(shift_bin_bestn4)
        sigb{3,i} = 0;
    else
        for j = 1:length(shift_bin_bestn4{i})
            sigb{3,i}(j) = M_bestn4(t,shift_bin_bestn4{i}(j));
        end
    end
end

%Dl1 del
t = round(p*length(T_deln1));
for i = 1:16
    if isempty(shift_bin_deln1)
        sigd{1,i} = 0;
    else
        for j = 1:length(shift_bin_deln1{i})
            sigd{1,i}(j) = M_deln1(t,shift_bin_deln1{i}(j));
        end
    end
end
t = round(p*length(T_deln2));
for i = 1:16
    if isempty(shift_bin_deln2)
        sigd{2,i} = 0;
    else
        for j = 1:length(shift_bin_deln2{i})
            sigd{2,i}(j) = M_deln2(t,shift_bin_deln2{i}(j));
        end
    end
end
t = round(p*length(T_deln4));
for i = 1:16
    if isempty(shift_bin_deln4)
        sigd{3,i} = 0;
    else
        for j = 1:length(shift_bin_deln4{i})
            sigd{3,i}(j) = M_deln4(t,shift_bin_deln4{i}(j));
        end
    end
end

%Dl2 best
t = round(p*length(T_dl2bestn1));
for i = 1:16
    if isempty(shift_bin_dl2bestn1)
        sigdl2b{1,i} = 0;
    else
        for j = 1:length(shift_bin_dl2bestn1{i})
            sigdl2b{1,i}(j) = M_dl2bestn1(t,shift_bin_dl2bestn1{i}(j));
        end
    end
end
t = round(p*length(T_dl2bestn2));
for i = 1:16
    if isempty(shift_bin_dl2bestn2)
        sigdl2b{2,i} = 0;
    else
        for j = 1:length(shift_bin_dl2bestn2{i})
            sigdl2b{2,i}(j) = M_dl2bestn2(t,shift_bin_dl2bestn2{i}(j));
        end
    end
end
t = round(p*length(T_dl2bestn5));
for i = 1:16
    if isempty(shift_bin_dl2bestn5)
        sigdl2b{3,i} = 0;
    else
        for j = 1:length(shift_bin_dl2bestn5{i})
            sigdl2b{3,i}(j) = M_dl2bestn5(t,shift_bin_dl2bestn5{i}(j));
        end
    end
end

%Dl2 del
t = round(p*length(T_dl2deln1));
for i = 1:16
    if isempty(shift_bin_dl2deln1)
        sigdl2d{1,i} = 0;
    else
        for j = 1:length(shift_bin_dl2deln1{i})
            sigdl2d{1,i}(j) = M_dl2deln1(t,shift_bin_dl2deln1{i}(j));
        end
    end
end
t = round(p*length(T_dl2deln2));
for i = 1:16
    if isempty(shift_bin_dl2deln2)
        sigdl2d{2,i} = 0;
    else
        for j = 1:length(shift_bin_dl2deln2{i})
            sigdl2d{2,i}(j) = M_dl2deln2(t,shift_bin_dl2deln2{i}(j));
        end
    end
end
t = round(p*length(T_dl2deln4));
for i = 1:16
    if isempty(shift_bin_dl2deln4)
        sigdl2d{3,i} = 0;
    else
        for j = 1:length(shift_bin_dl2deln4{i})
            sigdl2d{3,i}(j) = M_dl2deln4(t,shift_bin_dl2deln4{i}(j));
        end
    end
end

%Zld Good
t = round(p*length(T_zldgoodn2));
for i = 1:16
    if isempty(shift_bin_zldgoodn2)
        sigzg{1,i} = 0;
    else
        for j = 1:length(shift_bin_zldgoodn2{i})
            sigzg{1,i}(j) = M_zldgoodn2(t,shift_bin_zldgoodn2{i}(j));
        end
    end
end
t = round(p*length(T_zldgoodn3));
for i = 1:16
    if isempty(shift_bin_zldgoodn3)
        sigzg{2,i} = 0;
    else
        for j = 1:length(shift_bin_zldgoodn3{i})
            sigzg{2,i}(j) = M_zldgoodn3(t,shift_bin_zldgoodn3{i}(j));
        end
    end
end
t = round(p*length(T_zldgoodn4));
for i = 1:16
    if isempty(shift_bin_zldgoodn4)
        sigzg{3,i} = 0;
    else
        for j = 1:length(shift_bin_zldgoodn4{i})
            sigzg{3,i}(j) = M_zldgoodn4(t,shift_bin_zldgoodn4{i}(j));
        end
    end
end

for i = 1:3
    sigc{i,16} = [];
    sigb{i,16} = [];
    sigd{i,16} = [];
    sigdl2b{i,15} = [];
    sigdl2d{i,15} = [];
    sigzg{i,15} = [];
    sigdl2b{i,16} = [];
    sigdl2d{i,16} = [];
    sigzg{i,16} = [];
end

for i = 1:16
    tmp1{i} = [sigc{1,i} sigc{2,i} sigc{3,i} sigc{4,i} sigc{5,i} sigc{6,i}];
    tmp2{i} = [sigb{1,i} sigb{2,i} sigb{3,i}];
    tmp3{i} = [sigd{1,i} sigd{2,i} sigd{3,i}];
    tmp4{i} = [sigdl2b{1,i} sigdl2b{2,i} sigdl2b{3,i}];
    tmp5{i} = [sigdl2d{1,i} sigdl2d{2,i} sigdl2d{3,i}];
    tmp6{i} = [sigzg{1,i} sigzg{2,i} sigzg{3,i}];
end

for i = 1:16
    ave_sigc(i) = mean(tmp1{i});
    ave_sigb(i) = mean(tmp2{i});
    ave_sigd(i) = mean(tmp3{i});
    ave_sigdl2b(i) = mean(tmp4{i});
    ave_sigdl2d(i) = mean(tmp5{i});
    ave_sigzg(i) = mean(tmp6{i});

    std_sigc(i) = std(tmp1{i})./sqrt(length(tmp1{i}));
    std_sigb(i) = std(tmp2{i})./sqrt(length(tmp2{i}));
    std_sigd(i) = std(tmp3{i})./sqrt(length(tmp3{i}));
    std_sigdl2b(i) = std(tmp4{i})./sqrt(length(tmp4{i}));
    std_sigdl2d(i) = std(tmp5{i})./sqrt(length(tmp5{i}));
    std_sigzg(i) = std(tmp6{i})./sqrt(length(tmp6{i}));
end


%     h = figure;
%     errorbar([0:bins_bestn1-1]/bins_bestn1,ave_sigc,std_sigc,'color',rgb('black')); hold on
%     errorbar([0:bins_bestn1-1]/bins_bestn1,ave_sigb,std_sigb,'color',rgb('red')); 
%     errorbar([0:bins_bestn1-1]/bins_bestn1,ave_sigd,std_sigd,'--','color',rgb('lightcoral'));
%     errorbar([0:bins_bestn1-1]/bins_bestn1,ave_sigdl2b,std_sigdl2b,'color',rgb('blue')); 
%     errorbar([0:bins_bestn1-1]/bins_bestn1,ave_sigdl2d,std_sigdl2d,'--','color',rgb('cornflowerblue'));
%     title('Signal Across Embryo')
%     xlabel('Location')
%     ylabel('Nuclei Signal')
%     legend('control', 'dl1 strong','dl1 weak', 'dl2 strong', 'dl2 weak','Location', 'northwest')
%     ylim([0 3000])
%     xlim([0 .9375])
%     a = annotation('textbox',[.155 .10 .4 .57],'String',['% NC14 = ',num2str(thing(q)*100)],'FitBoxToText','on');
%     a.EdgeColor = 'none';
%     F(q) = getframe(gcf);
%     close(h)


clear tmp1 tmp2 tmp3 tmp4 tmp5 tmp6
% end

% video = VideoWriter('G:\My Drive\Research\Affinity Project\Signal_over_time.mp4');
% video.FrameRate = 1;
% video.Quality = 100;
% open(video)
% writeVideo(video,F);
% close(video)

figure;
errorbar([0:15]./16,ave_sigc,std_sigc,'Color',rgb('black')); hold on
errorbar([0:15]./16,ave_sigb,std_sigb,'Color',rgb('red'))
errorbar([0:15]./16,ave_sigd,std_sigd,'--','Color',rgb('lightcoral'))
title('Nuclei Signal Across Embryo at 70% NC14')
xlabel('Location')
ylabel('Nuclei Signal (A.U.)')
legend('control', 'Dl1 strong','Dl1 weak','Location', 'northwest')
xlim([0 .9375])
ylim([0 2.6e3])

figure;
errorbar([0:15]./16,ave_sigc,std_sigc,'Color',rgb('black')); hold on
errorbar([0:15]./16,ave_sigdl2b,std_sigdl2b,'Color',rgb('blue'))
errorbar([0:15]./16,ave_sigdl2d,std_sigdl2d,'--','Color',rgb('cornflowerblue'))
title('Nuclei Signal Across Embryo at 70% NC14')
xlabel('Location')
ylabel('Nuclei Signal (A.U.)')
legend('control', 'Dl2 strong','Dl2 weak','Location', 'northwest')
xlim([0 .9375])
ylim([0 2.6e3])

figure;
errorbar([0:15]./16,ave_sigc,std_sigc,'Color',rgb('black')); hold on
errorbar([0:15]./16,ave_sigb,std_sigb,'Color',rgb('red'))
errorbar([0:15]./16,ave_sigdl2b,std_sigdl2b,'Color',rgb('blue'))
errorbar([0:15]./16,ave_sigzg,std_sigzg,'Color',rgb('green'))
title('Nuclei Signal Across Embryo')
xlabel('Location')
ylabel('Nuclei Signal (A.U.)')
legend('control', 'Dl1 strong','Dl1 weak','Zld strong','Location', 'northwest')
xlim([0 .9375])
ylim([0 2.6e3])

%% Kinetics Rates
load('G:\Shared drives\Lim_Lab\Sam\Affinity_Project\kinetics_data.mat')

%Active Nuclei Only
for i = 1:16
    %For getting active nuclei
    [~,ic1{i}] = intersect(shift_bin_controln1{i},activeM_controln1);
    [~,ic2{i}] = intersect(shift_bin_controln2{i},activeM_controln2);
    [~,ic4{i}] = intersect(shift_bin_controln4{i},activeM_controln4);
    [~,ic3{i}] = intersect(shift_bin_controln3{i},activeM_controln3);
    [~,ic5{i}] = intersect(shift_bin_controln5{i},activeM_controln5);
    [~,ic6{i}] = intersect(shift_bin_controln6{i},activeM_controln6);
    [~,ib1{i}] = intersect(shift_bin_bestn1{i},activeM_bestn1);
    [~,ib2{i}] = intersect(shift_bin_bestn2{i},activeM_bestn2);
    [~,ib4{i}] = intersect(shift_bin_bestn4{i},activeM_bestn4);
    [~,id1{i}] = intersect(shift_bin_deln1{i},activeM_deln1);
    [~,id2{i}] = intersect(shift_bin_deln2{i},activeM_deln2);
    [~,id4{i}] = intersect(shift_bin_deln4{i},activeM_deln4);
    [~,idl2b1{i}] = intersect(shift_bin_dl2bestn1{i},activeM_dl2bestn1);
    [~,idl2b2{i}] = intersect(shift_bin_dl2bestn2{i},activeM_dl2bestn2);
    [~,idl2b5{i}] = intersect(shift_bin_dl2bestn5{i},activeM_dl2bestn5);
    [~,idl2d1{i}] = intersect(shift_bin_dl2deln1{i},activeM_dl2deln1);
    [~,idl2d2{i}] = intersect(shift_bin_dl2deln2{i},activeM_dl2deln2);
    [~,idl2d4{i}] = intersect(shift_bin_dl2deln4{i},activeM_dl2deln4);
    [~,izg2{i}] = intersect(shift_bin_zldgoodn2{i},activeM_zldgoodn2);
    [~,izg3{i}] = intersect(shift_bin_zldgoodn3{i},activeM_zldgoodn3);
    [~,izg4{i}] = intersect(shift_bin_zldgoodn4{i},activeM_zldgoodn4);
end

locations = [0:15]/16;
for i = 1:16
    [control.kon{1,i}] = controlreps.rep1.kons(shift_bin_controln1{i}(ic1{i}));
    [control.kon{2,i}] = controlreps.rep2.kons(shift_bin_controln2{i}(ic2{i}));
    [control.kon{3,i}] = controlreps.rep3.kons(shift_bin_controln4{i}(ic4{i}));
    [control.kon{4,i}] = controlreps.rep4.kons(shift_bin_controln3{i}(ic3{i}));
    [control.kon{5,i}] = controlreps.rep5.kons(shift_bin_controln5{i}(ic5{i}));
    [control.kon{6,i}] = controlreps.rep6.kons(shift_bin_controln6{i}(ic6{i}));
    [control.koff{1,i}] = controlreps.rep1.koffs(shift_bin_controln1{i}(ic1{i}));
    [control.koff{2,i}] = controlreps.rep2.koffs(shift_bin_controln2{i}(ic2{i}));
    [control.koff{3,i}] = controlreps.rep3.koffs(shift_bin_controln4{i}(ic4{i}));
    [control.koff{4,i}] = controlreps.rep4.koffs(shift_bin_controln3{i}(ic3{i}));
    [control.koff{5,i}] = controlreps.rep5.koffs(shift_bin_controln5{i}(ic5{i}));
    [control.koff{6,i}] = controlreps.rep6.koffs(shift_bin_controln6{i}(ic6{i}));
    [control.kinis{1,i}] = controlreps.rep1.kinis(shift_bin_controln1{i}(ic1{i}));
    [control.kinis{2,i}] = controlreps.rep2.kinis(shift_bin_controln2{i}(ic2{i}));
    [control.kinis{3,i}] = controlreps.rep3.kinis(shift_bin_controln4{i}(ic4{i}));
    [control.kinis{4,i}] = controlreps.rep4.kinis(shift_bin_controln3{i}(ic3{i}));
    [control.kinis{5,i}] = controlreps.rep5.kinis(shift_bin_controln5{i}(ic5{i}));
    [control.kinis{6,i}] = controlreps.rep6.kinis(shift_bin_controln6{i}(ic6{i}));
    [control.totmrna{1,i}] = controlreps.rep1.totmrna(shift_bin_controln1{i}(ic1{i}));
    [control.totmrna{2,i}] = controlreps.rep2.totmrna(shift_bin_controln2{i}(ic2{i}));
    [control.totmrna{3,i}] = controlreps.rep3.totmrna(shift_bin_controln4{i}(ic4{i}));
    [control.totmrna{4,i}] = controlreps.rep4.totmrna(shift_bin_controln3{i}(ic3{i}));
    [control.totmrna{5,i}] = controlreps.rep5.totmrna(shift_bin_controln5{i}(ic5{i}));
    [control.totmrna{6,i}] = controlreps.rep6.totmrna(shift_bin_controln6{i}(ic6{i}));
    
    [control.konmean(i)] = nanmean([control.kon{1,i} control.kon{2,i} control.kon{3,i} control.kon{4,i} control.kon{5,i} control.kon{6,i}]);
    [control.konstd(i)] = nanstd([control.kon{1,i} control.kon{2,i} control.kon{3,i} control.kon{4,i} control.kon{5,i} control.kon{6,i}])./sqrt(length(control.kon{1,i})+length(control.kon{2,i})+length(control.kon{3,i})+length(control.kon{4,i})+length(control.kon{5,i})+length(control.kon{6,i}));
    [control.koffmean(i)] = nanmean([control.koff{1,i} control.koff{2,i} control.koff{3,i} control.koff{4,i} control.koff{5,i} control.koff{6,i}]);
    [control.koffstd(i)] = nanstd([control.koff{1,i} control.koff{2,i} control.koff{3,i}])./sqrt(length(control.kon{1,i})+length(control.kon{2,i})+length(control.kon{3,i})+length(control.kon{4,i})+length(control.kon{5,i})+length(control.kon{6,i}));
    [control.kinismean(i)] = nanmean([control.kinis{1,i} control.kinis{2,i} control.kinis{3,i} control.kinis{4,i} control.kinis{5,i} control.kinis{6,i}]);
    [control.kinisstd(i)] = nanstd([control.kinis{1,i} control.kinis{2,i} control.kinis{3,i} control.kinis{4,i} control.kinis{5,i} control.kinis{6,i}])./sqrt(length(control.kon{1,i})+length(control.kon{2,i})+length(control.kon{3,i})+length(control.kon{4,i})+length(control.kon{5,i})+length(control.kon{6,i}));
    [control.totmrnamean(i)] = nanmean([control.totmrna{1,i} control.totmrna{2,i} control.totmrna{3,i} control.totmrna{4,i} control.totmrna{5,i} control.totmrna{6,i}]);
    [control.totmrnastd(i)] = nanstd([control.totmrna{1,i} control.totmrna{2,i} control.totmrna{3,i} control.totmrna{4,i} control.totmrna{5,i} control.totmrna{6,i}])./sqrt(length(control.kon{1,i})+length(control.kon{2,i})+length(control.kon{3,i})+length(control.kon{4,i})+length(control.kon{5,i})+length(control.kon{6,i}));
    
    %Dl1 best
    [dl1best.kon{1,i}] = dl1bestreps.rep1.kons(shift_bin_bestn1{i}(ib1{i}));
    [dl1best.kon{2,i}] = dl1bestreps.rep2.kons(shift_bin_bestn2{i}(ib2{i}));
    [dl1best.kon{3,i}] = dl1bestreps.rep3.kons(shift_bin_bestn4{i}(ib4{i}));
    [dl1best.koff{1,i}] = dl1bestreps.rep1.koffs(shift_bin_bestn1{i}(ib1{i}));
    [dl1best.koff{2,i}] = dl1bestreps.rep2.koffs(shift_bin_bestn2{i}(ib2{i}));
    [dl1best.koff{3,i}] = dl1bestreps.rep3.koffs(shift_bin_bestn4{i}(ib4{i}));
    [dl1best.kinis{1,i}] = dl1bestreps.rep1.kinis(shift_bin_bestn1{i}(ib1{i}));
    [dl1best.kinis{2,i}] = dl1bestreps.rep2.kinis(shift_bin_bestn2{i}(ib2{i}));
    [dl1best.kinis{3,i}] = dl1bestreps.rep3.kinis(shift_bin_bestn4{i}(ib4{i}));
    [dl1best.totmrna{1,i}] = dl1bestreps.rep1.totmrna(shift_bin_bestn1{i}(ib1{i}));
    [dl1best.totmrna{2,i}] = dl1bestreps.rep2.totmrna(shift_bin_bestn2{i}(ib2{i}));
    [dl1best.totmrna{3,i}] = dl1bestreps.rep3.totmrna(shift_bin_bestn4{i}(ib4{i}));
    
    [dl1best.konmean(i)] = nanmean([dl1best.kon{1,i} dl1best.kon{2,i} dl1best.kon{3,i}]);
    [dl1best.konstd(i)] = nanstd([dl1best.kon{1,i} dl1best.kon{2,i} dl1best.kon{3,i}])./sqrt(length(dl1best.kon{1,i})+length(dl1best.kon{2,i})+length(dl1best.kon{3,i}));
    [dl1best.koffmean(i)] = nanmean([dl1best.koff{1,i} dl1best.koff{2,i} dl1best.koff{3,i}]);
    [dl1best.koffstd(i)] = nanstd([dl1best.koff{1,i} dl1best.koff{2,i} dl1best.koff{3,i}])./sqrt(length(dl1best.kon{1,i})+length(dl1best.kon{2,i})+length(dl1best.kon{3,i}));
    [dl1best.kinismean(i)] = nanmean([dl1best.kinis{1,i} dl1best.kinis{2,i} dl1best.kinis{3,i}]);
    [dl1best.kinisstd(i)] = nanstd([dl1best.kinis{1,i} dl1best.kinis{2,i} dl1best.kinis{3,i}])./sqrt(length(dl1best.kon{1,i})+length(dl1best.kon{2,i})+length(dl1best.kon{3,i}));
    [dl1best.totmrnamean(i)] = nanmean([dl1best.totmrna{1,i} dl1best.totmrna{2,i} dl1best.totmrna{3,i}]);
    [dl1best.totmrnastd(i)] = nanstd([dl1best.totmrna{1,i} dl1best.totmrna{2,i} dl1best.totmrna{3,i}])./sqrt(length(dl1best.kon{1,i})+length(dl1best.kon{2,i})+length(dl1best.kon{3,i}));
    
    %Dl1 del
    [dl1del.kon{1,i}] = dl1delreps.rep1.kons(shift_bin_deln1{i}(id1{i}));
    [dl1del.kon{2,i}] = dl1delreps.rep2.kons(shift_bin_deln2{i}(id2{i}));
    [dl1del.kon{3,i}] = dl1delreps.rep3.kons(shift_bin_deln4{i}(id4{i}));
    [dl1del.koff{1,i}] = dl1delreps.rep1.koffs(shift_bin_deln1{i}(id1{i}));
    [dl1del.koff{2,i}] = dl1delreps.rep2.koffs(shift_bin_deln2{i}(id2{i}));
    [dl1del.koff{3,i}] = dl1delreps.rep3.koffs(shift_bin_deln4{i}(id4{i}));
    [dl1del.kinis{1,i}] = dl1delreps.rep1.kinis(shift_bin_deln1{i}(id1{i}));
    [dl1del.kinis{2,i}] = dl1delreps.rep2.kinis(shift_bin_deln2{i}(id2{i}));
    [dl1del.kinis{3,i}] = dl1delreps.rep3.kinis(shift_bin_deln4{i}(id4{i}));
    [dl1del.totmrna{1,i}] = dl1delreps.rep1.totmrna(shift_bin_deln1{i}(id1{i}));
    [dl1del.totmrna{2,i}] = dl1delreps.rep2.totmrna(shift_bin_deln2{i}(id2{i}));
    [dl1del.totmrna{3,i}] = dl1delreps.rep3.totmrna(shift_bin_deln4{i}(id4{i}));
    
    [dl1del.konmean(i)] = nanmean([dl1del.kon{1,i} dl1del.kon{2,i} dl1del.kon{3,i}]);
    [dl1del.konstd(i)] = nanstd([dl1del.kon{1,i} dl1del.kon{2,i} dl1del.kon{3,i}])./sqrt(length(dl1del.kon{1,i})+length(dl1del.kon{2,i})+length(dl1del.kon{3,i}));
    [dl1del.koffmean(i)] = nanmean([dl1del.koff{1,i} dl1del.koff{2,i} dl1del.koff{3,i}]);
    [dl1del.koffstd(i)] = nanstd([dl1del.koff{1,i} dl1del.koff{2,i} dl1del.koff{3,i}])./sqrt(length(dl1del.kon{1,i})+length(dl1del.kon{2,i})+length(dl1del.kon{3,i}));
    [dl1del.kinismean(i)] = nanmean([dl1del.kinis{1,i} dl1del.kinis{2,i} dl1del.kinis{3,i}]);
    [dl1del.kinisstd(i)] = nanstd([dl1del.kinis{1,i} dl1del.kinis{2,i} dl1del.kinis{3,i}])./sqrt(length(dl1del.kon{1,i})+length(dl1del.kon{2,i})+length(dl1del.kon{3,i}));
    [dl1del.totmrnamean(i)] = nanmean([dl1del.totmrna{1,i} dl1del.totmrna{2,i} dl1del.totmrna{3,i}]);
    [dl1del.totmrnastd(i)] = nanstd([dl1del.totmrna{1,i} dl1del.totmrna{2,i} dl1del.totmrna{3,i}])./sqrt(length(dl1del.kon{1,i})+length(dl1del.kon{2,i})+length(dl1del.kon{3,i}));
    
    %Dl2 best
    [dl2best.kon{1,i}] = dl2bestreps.rep1.kons(shift_bin_dl2bestn1{i}(idl2b1{i}));
    [dl2best.kon{2,i}] = dl2bestreps.rep2.kons(shift_bin_dl2bestn2{i}(idl2b2{i}));
    [dl2best.kon{3,i}] = dl2bestreps.rep3.kons(shift_bin_dl2bestn5{i}(idl2b5{i}));
    [dl2best.koff{1,i}] = dl2bestreps.rep1.koffs(shift_bin_dl2bestn1{i}(idl2b1{i}));
    [dl2best.koff{2,i}] = dl2bestreps.rep2.koffs(shift_bin_dl2bestn2{i}(idl2b2{i}));
    [dl2best.koff{3,i}] = dl2bestreps.rep3.koffs(shift_bin_dl2bestn5{i}(idl2b5{i}));
    [dl2best.kinis{1,i}] = dl2bestreps.rep1.kinis(shift_bin_dl2bestn1{i}(idl2b1{i}));
    [dl2best.kinis{2,i}] = dl2bestreps.rep2.kinis(shift_bin_dl2bestn2{i}(idl2b2{i}));
    [dl2best.kinis{3,i}] = dl2bestreps.rep3.kinis(shift_bin_dl2bestn5{i}(idl2b5{i}));
    [dl2best.totmrna{1,i}] = dl2bestreps.rep1.totmrna(shift_bin_dl2bestn1{i}(idl2b1{i}));
    [dl2best.totmrna{2,i}] = dl2bestreps.rep2.totmrna(shift_bin_dl2bestn2{i}(idl2b2{i}));
    [dl2best.totmrna{3,i}] = dl2bestreps.rep3.totmrna(shift_bin_dl2bestn5{i}(idl2b5{i}));
    
    [dl2best.konmean(i)] = nanmean([dl2best.kon{1,i} dl2best.kon{2,i} dl2best.kon{3,i}]);
    [dl2best.konstd(i)] = nanstd([dl2best.kon{1,i} dl2best.kon{2,i} dl2best.kon{3,i}])./sqrt(length(dl2best.kon{1,i})+length(dl2best.kon{2,i})+length(dl2best.kon{3,i}));
    [dl2best.koffmean(i)] = nanmean([dl2best.koff{1,i} dl2best.koff{2,i} dl2best.koff{3,i}]);
    [dl2best.koffstd(i)] = nanstd([dl2best.koff{1,i} dl2best.koff{2,i} dl2best.koff{3,i}])./sqrt(length(dl2best.kon{1,i})+length(dl2best.kon{2,i})+length(dl2best.kon{3,i}));
    [dl2best.kinismean(i)] = nanmean([dl2best.kinis{1,i} dl2best.kinis{2,i} dl2best.kinis{3,i}]);
    [dl2best.kinisstd(i)] = nanstd([dl2best.kinis{1,i} dl2best.kinis{2,i} dl2best.kinis{3,i}])./sqrt(length(dl2best.kon{1,i})+length(dl2best.kon{2,i})+length(dl2best.kon{3,i}));
    [dl2best.totmrnamean(i)] = nanmean([dl2best.totmrna{1,i} dl2best.totmrna{2,i} dl2best.totmrna{3,i}]);
    [dl2best.totmrnastd(i)] = nanstd([dl2best.totmrna{1,i} dl2best.totmrna{2,i} dl2best.totmrna{3,i}])./sqrt(length(dl2best.kon{1,i})+length(dl2best.kon{2,i})+length(dl2best.kon{3,i}));
    
    %Dl2 del
    [dl2del.kon{1,i}] = dl2delreps.rep1.kons(shift_bin_dl2deln1{i}(idl2d1{i}));
    [dl2del.kon{2,i}] = dl2delreps.rep2.kons(shift_bin_dl2deln2{i}(idl2d2{i}));
    [dl2del.kon{3,i}] = dl2delreps.rep3.kons(shift_bin_dl2deln4{i}(idl2d4{i}));
    [dl2del.koff{1,i}] = dl2delreps.rep1.koffs(shift_bin_dl2deln1{i}(idl2d1{i}));
    [dl2del.koff{2,i}] = dl2delreps.rep2.koffs(shift_bin_dl2deln2{i}(idl2d2{i}));
    [dl2del.koff{3,i}] = dl2delreps.rep3.koffs(shift_bin_dl2deln4{i}(idl2d4{i}));
    [dl2del.kinis{1,i}] = dl2delreps.rep1.kinis(shift_bin_dl2deln1{i}(idl2d1{i}));
    [dl2del.kinis{2,i}] = dl2delreps.rep2.kinis(shift_bin_dl2deln2{i}(idl2d2{i}));
    [dl2del.kinis{3,i}] = dl2delreps.rep3.kinis(shift_bin_dl2deln4{i}(idl2d4{i}));
    [dl2del.totmrna{1,i}] = dl2delreps.rep1.totmrna(shift_bin_dl2deln1{i}(idl2d1{i}));
    [dl2del.totmrna{2,i}] = dl2delreps.rep2.totmrna(shift_bin_dl2deln2{i}(idl2d2{i}));
    [dl2del.totmrna{3,i}] = dl2delreps.rep3.totmrna(shift_bin_dl2deln4{i}(idl2d4{i}));
    
    [dl2del.konmean(i)] = nanmean([dl2del.kon{1,i} dl2del.kon{2,i} dl2del.kon{3,i}]);
    [dl2del.konstd(i)] = nanstd([dl2del.kon{1,i} dl2del.kon{2,i} dl2del.kon{3,i}])./sqrt(length(dl2del.kon{1,i})+length(dl2del.kon{2,i})+length(dl2del.kon{3,i}));
    [dl2del.koffmean(i)] = nanmean([dl2del.koff{1,i} dl2del.koff{2,i} dl2del.koff{3,i}]);
    [dl2del.koffstd(i)] = nanstd([dl2del.koff{1,i} dl2del.koff{2,i} dl2del.koff{3,i}])./sqrt(length(dl2del.kon{1,i})+length(dl2del.kon{2,i})+length(dl2del.kon{3,i}));
    [dl2del.kinismean(i)] = nanmean([dl2del.kinis{1,i} dl2del.kinis{2,i} dl2del.kinis{3,i}]);
    [dl2del.kinisstd(i)] = nanstd([dl2del.kinis{1,i} dl2del.kinis{2,i} dl2del.kinis{3,i}])./sqrt(length(dl2del.kon{1,i})+length(dl2del.kon{2,i})+length(dl2del.kon{3,i}));
    [dl2del.totmrnamean(i)] = nanmean([dl2del.totmrna{1,i} dl2del.totmrna{2,i} dl2del.totmrna{3,i}]);
    [dl2del.totmrnastd(i)] = nanstd([dl2del.totmrna{1,i} dl2del.totmrna{2,i} dl2del.totmrna{3,i}])./sqrt(length(dl2del.kon{1,i})+length(dl2del.kon{2,i})+length(dl2del.kon{3,i}));
    
    %Zld good
    [zld.kon{1,i}] = zldbestreps.rep1.kons(shift_bin_zldgoodn2{i}(izg2{i}));
    [zld.kon{2,i}] = zldbestreps.rep2.kons(shift_bin_zldgoodn3{i}(izg3{i}));
    [zld.kon{3,i}] = zldbestreps.rep3.kons(shift_bin_zldgoodn4{i}(izg4{i}));
    [zld.koff{1,i}] = zldbestreps.rep1.koffs(shift_bin_zldgoodn2{i}(izg2{i}));
    [zld.koff{2,i}] = zldbestreps.rep2.koffs(shift_bin_zldgoodn3{i}(izg3{i}));
    [zld.koff{3,i}] = zldbestreps.rep3.koffs(shift_bin_zldgoodn4{i}(izg4{i}));
    [zld.kinis{1,i}] = zldbestreps.rep1.kinis(shift_bin_zldgoodn2{i}(izg2{i}));
    [zld.kinis{2,i}] = zldbestreps.rep2.kinis(shift_bin_zldgoodn3{i}(izg3{i}));
    [zld.kinis{3,i}] = zldbestreps.rep3.kinis(shift_bin_zldgoodn4{i}(izg4{i}));
    [zld.totmrna{1,i}] = zldbestreps.rep1.totmrna(shift_bin_zldgoodn2{i}(izg2{i}));
    [zld.totmrna{2,i}] = zldbestreps.rep2.totmrna(shift_bin_zldgoodn3{i}(izg3{i}));
    [zld.totmrna{3,i}] = zldbestreps.rep3.totmrna(shift_bin_zldgoodn4{i}(izg4{i}));
    
    [zld.konmean(i)] = nanmean([zld.kon{1,i} zld.kon{2,i} zld.kon{3,i}]);
    [zld.konstd(i)] = nanstd([zld.kon{1,i} zld.kon{2,i} zld.kon{3,i}])./sqrt(length(zld.kon{1,i})+length(zld.kon{2,i})+length(zld.kon{3,i}));
    [zld.koffmean(i)] = nanmean([zld.koff{1,i} zld.koff{2,i} zld.koff{3,i}]);
    [zld.koffstd(i)] = nanstd([zld.koff{1,i} zld.koff{2,i} zld.koff{3,i}])./sqrt(length(zld.kon{1,i})+length(zld.kon{2,i})+length(zld.kon{3,i}));
    [zld.kinismean(i)] = nanmean([zld.kinis{1,i} zld.kinis{2,i} zld.kinis{3,i}]);
    [zld.kinisstd(i)] = nanstd([zld.kinis{1,i} zld.kinis{2,i} zld.kinis{3,i}])./sqrt(length(zld.kon{1,i})+length(zld.kon{2,i})+length(zld.kon{3,i}));
    [zld.totmrnamean(i)] = nanmean([zld.totmrna{1,i} zld.totmrna{2,i} zld.totmrna{3,i}]);
    [zld.totmrnastd(i)] = nanstd([zld.totmrna{1,i} zld.totmrna{2,i} zld.totmrna{3,i}])./sqrt(length(zld.kon{1,i})+length(zld.kon{2,i})+length(zld.kon{3,i}));
    
end


%Kon
figure;
errorbar([0:15]./16,control.konmean,control.konstd,'Color',rgb('black'),'LineWidth',2); hold on
errorbar([0:15]./16,dl1best.konmean,dl1best.konstd,'Color',rgb('red'),'LineWidth',2)
errorbar([0:15]./16,dl1del.konmean,dl1del.konstd,'--','Color',rgb('lightcoral'),'LineWidth',2)
title('Kon')
xlabel('Location')
ylabel('Kon (A.U.)')
legend('control', 'Dl1 strong','Dl1 weak','Location', 'northwest')
xlim([0 .9375])
ylim([0.3 1])

figure;
errorbar([0:15]./16,control.konmean,control.konstd,'Color',rgb('black'),'LineWidth',2); hold on
errorbar([0:15]./16,dl2best.konmean,dl2best.konstd,'Color',rgb('blue'),'LineWidth',2)
errorbar([0:15]./16,dl2del.konmean,dl2del.konstd,'--','Color',rgb('cornflowerblue'),'LineWidth',2)
title('Kon')
xlabel('Location')
ylabel('Kon (A.U.)')
legend('control', 'Dl2 strong','Dl2 weak','Zld strong','Location', 'northwest')
xlim([0 .9375])
ylim([0.3 1])

figure;
errorbar([0:15]./16,dl1best.konmean,dl1best.konstd,'Color',rgb('red'),'LineWidth',2); hold on
errorbar([0:15]./16,dl2best.konmean,dl2best.konstd,'Color',rgb('blue'),'LineWidth',2)
errorbar([0:15]./16,zld.konmean,zld.konstd,'Color',rgb('green'),'LineWidth',2); hold on
title('Kon')
xlabel('Location')
ylabel('Kon (A.U.)')
legend('Dl1 strong', 'Dl2 strong','Zld strong','Location', 'northwest')
xlim([0 .9375])
ylim([0.3 1])


%Koff
figure;
errorbar([0:15]./16,control.koffmean,control.koffstd,'Color',rgb('black'),'LineWidth',2); hold on
errorbar([0:15]./16,dl1best.koffmean,dl1best.koffstd,'Color',rgb('red'),'LineWidth',2)
errorbar([0:15]./16,dl1del.koffmean,dl1del.koffstd,'--','Color',rgb('lightcoral'),'LineWidth',2)
title('Koff')
xlabel('Location')
ylabel('Koff (A.U.)')
legend('control', 'Dl1 strong','Dl1 weak','Location', 'northwest')
xlim([0 .9375])
ylim([0 .9])

figure;
errorbar([0:15]./16,control.koffmean,control.koffstd,'Color',rgb('black'),'LineWidth',2); hold on
errorbar([0:15]./16,dl2best.koffmean,dl2best.koffstd,'Color',rgb('blue'),'LineWidth',2)
errorbar([0:15]./16,dl2del.koffmean,dl2del.koffstd,'--','Color',rgb('cornflowerblue'),'LineWidth',2)
title('Koff')
xlabel('Location')
ylabel('Koff (A.U.)')
legend('control', 'Dl2 strong','Dl2 weak','Location', 'northwest')
xlim([0 .9375])
ylim([0 .9])

figure;
errorbar([0:15]./16,dl1best.koffmean,dl1best.koffstd,'Color',rgb('red'),'LineWidth',2); hold on
errorbar([0:15]./16,dl2best.koffmean,dl2best.koffstd,'Color',rgb('blue'),'LineWidth',2)
errorbar([0:15]./16,zld.koffmean,zld.koffstd,'Color',rgb('green'),'LineWidth',2); hold on
title('Koff')
xlabel('Location')
ylabel('Koff (A.U.)')
legend('Dl1 strong', 'Dl2 strong','Zld strong','Location', 'northwest')
xlim([0 .9375])
ylim([0 .9])

%Kinis
figure;
errorbar([0:15]./16,control.kinismean,control.kinisstd,'Color',rgb('black'),'LineWidth',2); hold on
errorbar([0:15]./16,dl1best.kinismean,dl1best.kinisstd,'Color',rgb('red'),'LineWidth',2)
errorbar([0:15]./16,dl1del.kinismean,dl1del.kinisstd,'--','Color',rgb('lightcoral'),'LineWidth',2)
title('Kini')
xlabel('Location')
ylabel('Kini (A.U.)')
legend('control', 'Dl1 strong','Dl1 weak','Location', 'northwest')
xlim([0 .9375])
ylim([0 2500])
 
figure;
errorbar([0:15]./16,control.kinismean,control.kinisstd,'Color',rgb('black'),'LineWidth',2); hold on
errorbar([0:15]./16,dl2best.kinismean,dl2best.kinisstd,'Color',rgb('blue'),'LineWidth',2)
errorbar([0:15]./16,dl2del.kinismean,dl2del.kinisstd,'--','Color',rgb('cornflowerblue'),'LineWidth',2)
title('Kini')
xlabel('Location')
ylabel('Kini (A.U.)')
legend('control', 'Dl2 strong','Dl2 weak','Location', 'northwest')
xlim([0 .9375])
ylim([0 2400])

figure;
errorbar([0:15]./16,dl1best.kinismean,dl1best.kinisstd,'Color',rgb('red'),'LineWidth',2); hold on
errorbar([0:15]./16,dl2best.kinismean,dl2best.kinisstd,'Color',rgb('blue'),'LineWidth',2)
errorbar([0:15]./16,zld.kinismean,zld.kinisstd,'Color',rgb('green'),'LineWidth',2); hold on
title('Kini')
xlabel('Location')
ylabel('Kini (A.U.)')
legend('Dl1 strong', 'Dl2 strong','Zld strong','Location', 'northwest')
xlim([0 .9375])
ylim([0 2500])


%Total mRNA
figure;
errorbar([0:15]./16,dl1best.totmrnamean,dl1best.totmrnastd,'Color',rgb('red')); hold on
errorbar([0:15]./16,dl2best.totmrnamean,dl2best.totmrnastd,'Color',rgb('blue'))
errorbar([0:15]./16,zld.totmrnamean,zld.totmrnastd,'Color',rgb('green')); hold on
title('Kini')
xlabel('Location')
ylabel('Kini (A.U.)')
legend('Dl1 strong', 'Dl2 strong','Zld strong','Location', 'northwest')
xlim([0 .9375])


subplot(1,3,1)
scatter(linxc1(80, :), linyc1(80, :), 40, controlreps.rep1.kons./(max(controlreps.rep1.kons)), 'filled')
title('Kon')
subplot(1,3,2)
scatter(linxc1(80, :), linyc1(80, :), 40, controlreps.rep1.koffs./(max(controlreps.rep1.koffs)), 'filled')
title('Koff')
subplot(1,3,3)
scatter(linxc1(80, :), linyc1(80, :), 40, controlreps.rep1.kinis./(max(controlreps.rep1.kinis)), 'filled')
title('Kini')

%Only Active Nuclei
subplot(1,3,1)
scatter(linxc1(80, activeM_controln1), linyc1(80, activeM_controln1), 40, controlreps.rep1.kons(activeM_controln1)./(max(controlreps.rep1.kons)), 'filled')
title('Kon')
subplot(1,3,2)
scatter(linxc1(80, activeM_controln1), linyc1(80, activeM_controln1), 40, controlreps.rep1.koffs(activeM_controln1)./(max(controlreps.rep1.koffs)), 'filled')
title('Koff')
subplot(1,3,3)
scatter(linxc1(80, activeM_controln1), linyc1(80, activeM_controln1), 40, controlreps.rep1.kinis(activeM_controln1)./(max(controlreps.rep1.kinis)), 'filled')
title('Kini')

%% Average Trace
normtracec1 = interp1(time_controln1,M_controln1,norm_time);
normtracec2 = interp1(time_controln2,M_controln2,norm_time);
normtracec4 = interp1(time_controln4,M_controln4,norm_time);
normtracec3 = interp1(time_controln3,M_controln3,norm_time);
normtracec5 = interp1(time_controln5,M_controln5,norm_time);
normtracec6 = interp1(time_controln6,M_controln6,norm_time);

normtraceb1 = interp1(time_bestn1,M_bestn1,norm_time);
normtraceb2 = interp1(time_bestn2,M_bestn2,norm_time);
normtraceb4 = interp1(time_bestn4,M_bestn4,norm_time);

normtraced1 = interp1(time_deln1,M_deln1,norm_time);
normtraced2 = interp1(time_deln2,M_deln2,norm_time);
normtraced4 = interp1(time_deln4,M_deln4,norm_time);

normtracedl2b1 = interp1(time_dl2bestn1,M_dl2bestn1,norm_time);
normtracedl2b2 = interp1(time_dl2bestn2,M_dl2bestn2,norm_time);
normtracedl2b5 = interp1(time_dl2bestn5,M_dl2bestn5,norm_time);

normtracedl2d1 = interp1(time_dl2deln1,M_dl2deln1,norm_time);
normtracedl2d2 = interp1(time_dl2deln2,M_dl2deln2,norm_time);
normtracedl2d4 = interp1(time_dl2deln4,M_dl2deln4,norm_time);

normtracezg2 = interp1(time_zldgoodn2,M_zldgoodn2,norm_time);
normtracezg3 = interp1(time_zldgoodn3,M_zldgoodn3,norm_time);
normtracezg4 = interp1(time_zldgoodn4,M_zldgoodn4,norm_time);

%Revision
normtraceRc1 = interp1(time_Rcontroln1,M_Rcn1,norm_time);
normtraceRc2 = interp1(time_Rcontroln2,M_Rcn2,norm_time);
normtraceRc3 = interp1(time_Rcontroln3,M_Rcn3,norm_time);
normtraceRdz1 = interp1(time_Rdzn1,M_Rdzn1,norm_time);
normtraceRdz2 = interp1(time_Rdzn2,M_Rdzn2,norm_time);
normtraceRdz3 = interp1(time_Rdzn3,M_Rdzn3,norm_time);
normtraceRd1 = interp1(time_Rdn1,M_Rdn1,norm_time);
normtraceRd2 = interp1(time_Rdn2,M_Rdn2,norm_time);
normtraceRd3 = interp1(time_Rdn3,M_Rdn3,norm_time);

%Active Nuclei Only
for i = 1:16
    %For getting active nuclei
    [ic1{i}] = intersect(shift_bin_controln1{i},activeM_controln1);
    [ic2{i}] = intersect(shift_bin_controln2{i},activeM_controln2);
    [ic4{i}] = intersect(shift_bin_controln4{i},activeM_controln4);
    [ic3{i}] = intersect(shift_bin_controln3{i},activeM_controln3);
    [ic5{i}] = intersect(shift_bin_controln5{i},activeM_controln5);
    [ic6{i}] = intersect(shift_bin_controln6{i},activeM_controln6);
    [ib1{i}] = intersect(shift_bin_bestn1{i},activeM_bestn1);
    [ib2{i}] = intersect(shift_bin_bestn2{i},activeM_bestn2);
    [ib4{i}] = intersect(shift_bin_bestn4{i},activeM_bestn4);
    [id1{i}] = intersect(shift_bin_deln1{i},activeM_deln1);
    [id2{i}] = intersect(shift_bin_deln2{i},activeM_deln2);
    [id4{i}] = intersect(shift_bin_deln4{i},activeM_deln4);
    [idl2b1{i}] = intersect(shift_bin_dl2bestn1{i},activeM_dl2bestn1);
    [idl2b2{i}] = intersect(shift_bin_dl2bestn2{i},activeM_dl2bestn2);
    [idl2b5{i}] = intersect(shift_bin_dl2bestn5{i},activeM_dl2bestn5);
    [idl2d1{i}] = intersect(shift_bin_dl2deln1{i},activeM_dl2deln1);
    [idl2d2{i}] = intersect(shift_bin_dl2deln2{i},activeM_dl2deln2);
    [idl2d4{i}] = intersect(shift_bin_dl2deln4{i},activeM_dl2deln4);
    [izg2{i}] = intersect(shift_bin_zldgoodn2{i},activeM_zldgoodn2);
    [izg3{i}] = intersect(shift_bin_zldgoodn3{i},activeM_zldgoodn3);
    [izg4{i}] = intersect(shift_bin_zldgoodn4{i},activeM_zldgoodn4);
    
%     %Revision
    [iRc1{i}] = intersect(shift_bin_Rcn1{i},activeM_Rcn1);
    [iRc2{i}] = intersect(shift_bin_Rcn2{i},activeM_Rcn2);
    [iRc3{i}] = intersect(shift_bin_Rcn3{i},activeM_Rcn3);
    [iRdz1{i}] = intersect(shift_bin_Rdzn1{i},activeM_Rdzn1);
    [iRdz2{i}] = intersect(shift_bin_Rdzn2{i},activeM_Rdzn2);
    [iRdz3{i}] = intersect(shift_bin_Rdzn3{i},activeM_Rdzn3);
    [iRd1{i}] = intersect(shift_bin_Rdn1{i},activeM_Rdn1);
    [iRd2{i}] = intersect(shift_bin_Rdn2{i},activeM_Rdn2);
    [iRd3{i}] = intersect(shift_bin_Rdn3{i},activeM_Rdn3);

end
for i = 1:8
    

    for j = 1:181
        trace_c{i}(j) = mean([normtracec1(j,ic1{i}) normtracec1(j,ic1{17-i}) normtracec2(j,ic2{i}) normtracec2(j,ic2{17-i}) normtracec4(j,ic4{i}) normtracec4(j,ic4{17-i}) normtracec3(j,ic3{i}) normtracec3(j,ic3{17-i}) normtracec5(j,ic5{i}) normtracec5(j,ic5{17-i}) normtracec6(j,ic6{i}) normtracec6(j,ic6{17-i})]);
        trace_b{i}(j) = mean([normtraceb1(j,ib1{i}) normtraceb1(j,ib1{17-i}) normtraceb2(j,ib2{i}) normtraceb2(j,ib2{17-i}) normtraceb4(j,ib4{i}) normtraceb4(j,ib4{17-i})]);
        trace_d{i}(j) = mean([normtraced1(j,id1{i}) normtraced1(j,id1{17-i}) normtraced2(j,id2{i}) normtraced2(j,id2{17-i}) normtraced4(j,id4{i}) normtraced4(j,id4{17-i})]);
        trace_dl2b{i}(j) = mean([normtracedl2b1(j,idl2b1{i}) normtracedl2b1(j,idl2b1{17-i}) normtracedl2b2(j,idl2b2{i}) normtracedl2b2(j,idl2b2{17-i}) normtracedl2b5(j,idl2b5{i}) normtracedl2b5(j,idl2b5{17-i})]);
        trace_dl2d{i}(j) = mean([normtracedl2d1(j,idl2d1{i}) normtracedl2d1(j,idl2d1{17-i}) normtracedl2d2(j,idl2d2{i}) normtracedl2d2(j,idl2d2{17-i}) normtracedl2d4(j,idl2d4{i}) normtracedl2d4(j,idl2d4{17-i})]);
        trace_zg{i}(j) = mean([normtracezg2(j,izg2{i}) normtracezg2(j,izg2{17-i}) normtracezg3(j,izg3{i}) normtracezg3(j,izg3{17-i}) normtracezg4(j,izg4{i}) normtracezg4(j,izg4{17-i})]);

        stdtrace_c{i}(j) = std([normtracec1(j,ic1{i}) normtracec1(j,ic1{17-i}) normtracec2(j,ic2{i}) normtracec2(j,ic2{17-i}) normtracec4(j,ic4{i}) normtracec4(j,ic4{17-i}) normtracec3(j,ic3{i}) normtracec3(j,ic3{17-i}) normtracec5(j,ic5{i}) normtracec5(j,ic5{17-i}) normtracec6(j,ic6{i}) normtracec6(j,ic6{17-i})])./sqrt(length(ic1{i})+length(ic1{17-i})+length(ic2{i})+length(ic2{17-i})+length(ic4{i})+length(ic4{17-i})+length(ic3{i})+length(ic3{17-i})+length(ic5{i})+length(ic5{17-i})+length(ic6{i})+length(ic6{17-i}));
        stdtrace_b{i}(j) = std([normtraceb1(j,ib1{i}) normtraceb1(j,ib1{17-i}) normtraceb2(j,ib2{i}) normtraceb2(j,ib2{17-i}) normtraceb4(j,ib4{i}) normtraceb4(j,ib4{17-i})])./sqrt(length(ib1{i})+length(ib1{17-i})+length(ib2{i})+length(ib2{17-i})+length(ib4{i})+length(ib4{17-i}));
        stdtrace_d{i}(j) = std([normtraced1(j,id1{i}) normtraced1(j,id1{17-i}) normtraced2(j,id2{i}) normtraced2(j,id2{17-i}) normtraced4(j,id4{i}) normtraced4(j,id4{17-i})])./sqrt(length(id1{i})+length(id1{17-i})+length(id2{i})+length(id2{17-i})+length(id4{i})+length(id4{17-i}));
        stdtrace_dl2b{i}(j) = std([normtracedl2b1(j,idl2b1{i}) normtracedl2b1(j,idl2b1{17-i}) normtracedl2b2(j,idl2b2{i}) normtracedl2b2(j,idl2b2{17-i}) normtracedl2b5(j,idl2b5{i}) normtracedl2b5(j,idl2b5{17-i})])./sqrt(length(idl2b1{i})+length(idl2b1{17-i})+length(idl2b2{i})+length(idl2b2{17-i})+length(idl2b5{i})+length(idl2b5{17-i}));
        stdtrace_dl2d{i}(j) = std([normtracedl2d1(j,idl2d1{i}) normtracedl2d1(j,idl2d1{17-i}) normtracedl2d2(j,idl2d2{i}) normtracedl2d2(j,idl2d2{17-i}) normtracedl2d4(j,idl2d4{i}) normtracedl2d4(j,idl2d4{17-i})])./sqrt(length(idl2d1{i})+length(idl2d1{17-i})+length(idl2d2{i})+length(idl2d2{17-i})+length(idl2d4{i})+length(idl2d4{17-i}));
        stdtrace_zg{i}(j) = std([normtracezg2(j,izg2{i}) normtracezg2(j,izg2{17-i}) normtracezg3(j,izg3{i}) normtracezg3(j,izg3{17-i}) normtracezg4(j,izg4{i}) normtracezg4(j,izg4{17-i})])./sqrt(length(izg2{i})+length(izg2{17-2})+length(izg3{i})+length(izg3{17-2})+length(izg4{i})+length(izg4{17-2}));

        
        %Revision
        trace_Rc{i}(j) = mean([normtraceRc1(j,iRc1{i}) normtraceRc1(j,iRc1{17-i}) normtraceRc2(j,iRc2{i}) normtraceRc2(j,iRc2{17-i}) normtraceRc3(j,iRc3{i}) normtraceRc3(j,iRc3{17-i})]);
        stdtrace_Rc{i}(j) = std([normtraceRc1(j,iRc1{i}) normtraceRc1(j,iRc1{17-i}) normtraceRc2(j,iRc2{i}) normtraceRc2(j,iRc2{17-i}) normtraceRc3(j,iRc3{i}) normtraceRc3(j,iRc3{17-i})])./sqrt(length(iRc1{i})+length(iRc1{17-i})+length(iRc2{i})+length(iRc2{17-i})+length(iRc3{i})+length(iRc3{17-i}));
        trace_Rdz{i}(j) = mean([normtraceRdz1(j,iRdz1{i}) normtraceRdz1(j,iRdz1{17-i}) normtraceRdz2(j,iRdz2{i}) normtraceRdz2(j,iRdz2{17-i}) normtraceRdz3(j,iRdz3{i}) normtraceRdz3(j,iRdz3{17-i})]);
        stdtrace_Rdz{i}(j) = std([normtraceRdz1(j,iRdz1{i}) normtraceRdz1(j,iRdz1{17-i}) normtraceRdz2(j,iRdz2{i}) normtraceRdz2(j,iRdz2{17-i}) normtraceRdz3(j,iRdz3{i}) normtraceRdz3(j,iRdz3{17-i})])./sqrt(length(iRdz1{i})+length(iRdz1{17-i})+length(iRdz2{i})+length(iRdz2{17-i})+length(iRdz3{i})+length(iRdz3{17-i}));
        trace_Rd{i}(j) = mean([normtraceRd1(j,iRd1{i}) normtraceRd1(j,iRd1{17-i}) normtraceRd2(j,iRd2{i}) normtraceRd2(j,iRd2{17-i}) normtraceRd3(j,iRd3{i}) normtraceRd3(j,iRd3{17-i})]);
        stdtrace_Rd{i}(j) = std([normtraceRd1(j,iRd1{i}) normtraceRd1(j,iRd1{17-i}) normtraceRd2(j,iRd2{i}) normtraceRd2(j,iRd2{17-i}) normtraceRd3(j,iRd3{i}) normtraceRd3(j,iRd3{17-i})])./sqrt(length(iRd1{i})+length(iRd1{17-i})+length(iRd2{i})+length(iRd2{17-i})+length(iRd3{i})+length(iRd3{17-i}));

    end
end

figure; 
for i = 1:8
subplot(2,4,i)
shadedErrorBar(norm_time,trace_c{i},stdtrace_c{i},{'k','LineWidth',2}); ylim([0 2000]); hold on
shadedErrorBar(norm_time,trace_b{i},stdtrace_b{i},{'r','LineWidth',2}); ylim([0 2000])
% shadedErrorBar(norm_time,trace_d{i},stdtrace_d{i},{'color', rgb('lightcoral'),'LineWidth',2}); ylim([0 2000])
shadedErrorBar(norm_time,trace_dl2b{i},stdtrace_dl2b{i},{'b','LineWidth',2}); ylim([0 3000])
% shadedErrorBar(norm_time,trace_dl2d{i},stdtrace_dl2d{i},{'color', rgb('cornflowerblue'),'LineWidth',2}); ylim([0 3000])
shadedErrorBar(norm_time,trace_zg{i},stdtrace_zg{i},{'color',rgb('green'),'LineWidth',2}); ylim([0 3000])

end


%Revision
figure; 
for i = 1:8
subplot(2,4,i)
shadedErrorBar(norm_time,trace_c{i},stdtrace_c{i},{'Color',rgb('grey'),'LineWidth',2}); ylim([0 2000]); hold on
shadedErrorBar(norm_time,trace_Rc{i},stdtrace_Rc{i},{'Color',rgb('black'),'LineWidth',2}); ylim([0 2000]); hold on
shadedErrorBar(norm_time,trace_Rdz{i},stdtrace_Rdz{i},{'Color',rgb('red'),'LineWidth',2}); ylim([0 2000]); hold on
shadedErrorBar(norm_time,trace_Rd{i},stdtrace_Rd{i},{'Color',rgb('blue'),'LineWidth',2}); ylim([0 2000]); hold on

end
% for i = 1:8
%     for j = 1:181
%         if isnan(trace_c{i}(j))
%             trace_c{i}(j) = 0;
%             trace_b{i}(j) = 0;
%             trace_dl2b{i}(j) = 0;
%             trace_zg{i}(j) = 0;
%         end
%     end
%     out_c(i) = trapz(norm_time,[trace_c{i}]);
%     out_b(i) = trapz(norm_time,[trace_b{i}]);
%     out_dl2b(i) = trapz(norm_time,[trace_dl2b{i}]);
%     out_zg(i) = trapz(norm_time,[trace_zg{i}]);
% end
% figure;plot(1:8,out_c,'ko'); hold on
% plot(1:8,out_b,'ro')
%% Active Nuclei per bin
for i = 1:16
    activebinc(i) = (length(ic1{i})+ length(ic2{i})+length(ic3{i})+length(ic4{i})+length(ic5{i})+length(ic6{i}))/6;
    activebinb(i) = (length(ib1{i})+ length(ib2{i})+length(ib4{i}))/3;
    activebind(i) = (length(id1{i})+ length(id2{i})+length(id4{i}))/3;
    activebindl2b(i) = (length(idl2b1{i})+ length(idl2b2{i})+length(idl2b5{i}))/3;
    activebindl2d(i) = (length(idl2d1{i})+ length(idl2d2{i})+length(idl2d4{i}))/3;
    activebinzg(i) = (length(izg2{i})+ length(izg3{i})+length(izg4{i}))/3;
    
    stdbinc(i) = std([length(ic1{i}) length(ic2{i}) length(ic3{i}) length(ic4{i}) length(ic5{i}) length(ic6{i})]);
    stdbinb(i) = std([length(ib1{i}) length(ib2{i}) length(ib4{i})]);
    stdbind(i) = std([length(id1{i}) length(id2{i}) length(id4{i})]);
    stdbindl2b(i) = std([length(idl2b1{i}) length(idl2b2{i}) length(idl2b5{i})]);
    stdbindl2d(i) = std([length(idl2d1{i}) length(idl2d2{i}) length(idl2d4{i})]);
    stdbinzg(i) = std([length(izg2{i}) length(izg3{i}) length(izg4{i})]);

end

figure;
errorbar([0:15]/16, activebinc,stdbinc,'Color',rgb('black'),'LineWidth',2); hold on
errorbar([0:15]/16, activebinb,stdbinb,'Color',rgb('red'),'LineWidth',2);
errorbar([0:15]/16, activebind,stdbind,'--','Color',rgb('lightcoral'),'LineWidth',2);
title('Average Active Nuclei Per Bin')

figure;
errorbar([0:15]/16, activebinc,stdbinc,'Color',rgb('black'),'LineWidth',2); hold on
errorbar([0:15]/16, activebindl2b,stdbindl2b,'Color',rgb('blue'),'LineWidth',2);
errorbar([0:15]/16, activebindl2d,stdbindl2d,'--','Color',rgb('cornflowerblue'),'LineWidth',2);

figure;
errorbar([0:15]/16, activebinc,stdbinc,'Color',rgb('black'),'LineWidth',2); hold on
errorbar([0:15]/16, activebinzg,stdbinzg,'Color',rgb('green'),'LineWidth',2);

%% Trace Heatmap
%{
heatmap_c = [];
for i = 1:16
    for j = 1:length(ic1{i})
        heatmap_c = [heatmap_c M_controln1(:,ic1{i}(j))];
    end
end
heatmap_c = heatmap_c';

for i = 1:16
    clock(i,:) = trace_c{i};
end
figure;heatmap(heatmap_c,'Colormap',jet)
figure;heatmap(clock,'Colormap',jet)
%}
%% Mean Active Signal
for i = 1:length(T_start_c{1})
    [~,medidx_c{1}(i)] = intersect(time_controln1,T_start_c{1}(i));
    activesignal_c{1}{i} = M_controln1(medidx_c{1}(i):end,activeM_controln1(i));
end
for i = 1:length(T_start_c{2})
    [~,medidx_c{2}(i)] = intersect(time_controln2,T_start_c{2}(i));
    activesignal_c{2}{i} = M_controln2(medidx_c{2}(i):end,activeM_controln2(i));
end
for i = 1:length(T_start_c{3})
    [~,medidx_c{3}(i)] = intersect(time_controln4,T_start_c{3}(i));
    activesignal_c{3}{i} = M_controln4(medidx_c{3}(i):end,activeM_controln4(i));
end
for i = 1:length(T_start_c{4})
    [~,medidx_c{4}(i)] = intersect(time_controln3,T_start_c{4}(i));
    activesignal_c{4}{i} = M_controln3(medidx_c{4}(i):end,activeM_controln3(i));
end
for i = 1:length(T_start_c{5})
    [~,medidx_c{5}(i)] = intersect(time_controln5,T_start_c{5}(i));
    activesignal_c{5}{i} = M_controln5(medidx_c{5}(i):end,activeM_controln5(i));
end
for i = 1:length(T_start_c{6})
    [~,medidx_c{6}(i)] = intersect(time_controln6,T_start_c{6}(i));
    activesignal_c{6}{i} = M_controln6(medidx_c{6}(i):end,activeM_controln6(i));
end

for i = 1:length(T_start_b{1})
    [~,medidx_b{1}(i)] = intersect(time_bestn1,T_start_b{1}(i));
    activesignal_b{1}{i} = M_bestn1(medidx_b{1}(i):end,activeM_bestn1(i));
end
for i = 1:length(T_start_b{2})
    [~,medidx_b{2}(i)] = intersect(time_bestn2,T_start_b{2}(i));
    activesignal_b{2}{i} = M_bestn2(medidx_b{2}(i):end,activeM_bestn2(i));
end
for i = 1:length(T_start_b{3})
    [~,medidx_b{3}(i)] = intersect(time_bestn4,T_start_b{3}(i));
    activesignal_b{3}{i} = M_bestn4(medidx_b{3}(i):end,activeM_bestn4(i));
end

for i = 1:length(T_start_d{1})
    [~,medidx_d{1}(i)] = intersect(time_deln1,T_start_d{1}(i));
    activesignal_d{1}{i} = M_deln1(medidx_d{1}(i):end,activeM_deln1(i));
end
for i = 1:length(T_start_d{2})
    [~,medidx_d{2}(i)] = intersect(time_deln2,T_start_d{2}(i));
    activesignal_d{2}{i} = M_deln2(medidx_d{2}(i):end,activeM_deln2(i));
end
for i = 1:length(T_start_d{3})
    [~,medidx_d{3}(i)] = intersect(time_deln4,T_start_d{3}(i));
    activesignal_d{3}{i} = M_deln4(medidx_d{3}(i):end,activeM_deln4(i));
end

for i = 1:length(T_start_dl2b{1})
    [~,medidx_dl2b{1}(i)] = intersect(time_dl2bestn1,T_start_dl2b{1}(i));
    activesignal_dl2b{1}{i} = M_dl2bestn1(medidx_dl2b{1}(i):end,activeM_dl2bestn1(i));
end
for i = 1:length(T_start_dl2b{2})
    [~,medidx_dl2b{2}(i)] = intersect(time_dl2bestn2,T_start_dl2b{2}(i));
    activesignal_dl2b{2}{i} = M_dl2bestn2(medidx_dl2b{2}(i):end,activeM_dl2bestn2(i));
end
for i = 1:length(T_start_dl2b{3})
    [~,medidx_dl2b{3}(i)] = intersect(time_dl2bestn5,T_start_dl2b{3}(i));
    activesignal_dl2b{3}{i} = M_dl2bestn5(medidx_dl2b{3}(i):end,activeM_dl2bestn5(i));
end

for i = 1:length(T_start_dl2d{1})
    [~,medidx_dl2d{1}(i)] = intersect(time_dl2deln1,T_start_dl2d{1}(i));
    activesignal_dl2d{1}{i} = M_dl2deln1(medidx_dl2d{1}(i):end,activeM_dl2deln1(i));
end
for i = 1:length(T_start_dl2d{2})
    [~,medidx_dl2d{2}(i)] = intersect(time_dl2deln2,T_start_dl2d{2}(i));
    activesignal_dl2d{2}{i} = M_dl2deln2(medidx_dl2d{2}(i):end,activeM_dl2deln2(i));
end
for i = 1:length(T_start_dl2d{3})
    [~,medidx_dl2d{3}(i)] = intersect(time_dl2deln4,T_start_dl2d{3}(i));
    activesignal_dl2d{3}{i} = M_dl2deln4(medidx_dl2d{3}(i):end,activeM_dl2deln4(i));
end

for i = 1:length(T_start_zg{1})
    [~,medidx_zg{1}(i)] = intersect(time_zldgoodn2,T_start_zg{1}(i));
    activesignal_zg{1}{i} = M_zldgoodn2(medidx_zg{1}(i):end,activeM_zldgoodn2(i));
end
for i = 1:length(T_start_zg{2})
    [~,medidx_zg{2}(i)] = intersect(time_zldgoodn3,T_start_zg{2}(i));
    activesignal_zg{2}{i} = M_zldgoodn3(medidx_zg{2}(i):end,activeM_zldgoodn3(i));
end
for i = 1:length(T_start_zg{3})
    [~,medidx_zg{3}(i)] = intersect(time_zldgoodn4,T_start_zg{3}(i));
    activesignal_zg{3}{i} = M_zldgoodn4(medidx_zg{3}(i):end,activeM_zldgoodn4(i));
end


%Revision
for i = 1:length(T_start_Rc{1})
    [~,medidx_Rc{1}(i)] = intersect(time_Rcontroln1,T_start_Rc{1}(i));
    activesignal_Rc{1}{i} = M_Rcn1(medidx_Rc{1}(i):end,activeM_Rcn1(i));
end
for i = 1:length(T_start_Rc{2})
    [~,medidx_Rc{2}(i)] = intersect(time_Rcontroln2,T_start_Rc{2}(i));
    activesignal_Rc{2}{i} = M_Rcn2(medidx_Rc{2}(i):end,activeM_Rcn2(i));
end
for i = 1:length(T_start_Rc{3})
    [~,medidx_Rc{3}(i)] = intersect(time_Rcontroln3,T_start_Rc{3}(i));
    activesignal_Rc{3}{i} = M_Rcn3(medidx_Rc{3}(i):end,activeM_Rcn3(i));
end
for i = 1:length(T_start_Rdz{1})
    [~,medidx_Rdz{1}(i)] = intersect(time_Rdzn1,T_start_Rdz{1}(i));
    activesignal_Rdz{1}{i} = M_Rdzn1(medidx_Rdz{1}(i):end,activeM_Rdzn1(i));
end
for i = 1:length(T_start_Rdz{2})
    [~,medidx_Rdz{2}(i)] = intersect(time_Rdzn2,T_start_Rdz{2}(i));
    activesignal_Rdz{2}{i} = M_Rdzn2(medidx_Rdz{2}(i):end,activeM_Rdzn2(i));
end
for i = 1:length(T_start_Rdz{3})
    [~,medidx_Rdz{3}(i)] = intersect(time_Rdzn3,T_start_Rdz{3}(i));
    activesignal_Rdz{3}{i} = M_Rdzn3(medidx_Rdz{3}(i):end,activeM_Rdzn3(i));
end
for i = 1:length(T_start_Rd{1})
    [~,medidx_Rd{1}(i)] = intersect(time_Rdn1,T_start_Rd{1}(i));
    activesignal_Rd{1}{i} = M_Rdn1(medidx_Rd{1}(i):end,activeM_Rdn1(i));
end
for i = 1:length(T_start_Rd{2})
    [~,medidx_Rd{2}(i)] = intersect(time_Rdn2,T_start_Rd{2}(i));
    activesignal_Rd{2}{i} = M_Rdn2(medidx_Rd{2}(i):end,activeM_Rdn2(i));
end
for i = 1:length(T_start_Rd{3})
    [~,medidx_Rd{3}(i)] = intersect(time_Rdn3,T_start_Rd{3}(i));
    activesignal_Rd{3}{i} = M_Rdn3(medidx_Rd{3}(i):end,activeM_Rdn3(i));
end
for i = 1:6
    for j = 1:length(activesignal_c{i})
        ave_activesignal_c{i}(j) = mean(activesignal_c{i}{j});
    end
end

for i = 1:3
    for j = 1:length(activesignal_b{i})
        ave_activesignal_b{i}(j) = mean(activesignal_b{i}{j});
    end
    for j = 1:length(activesignal_d{i})
        ave_activesignal_d{i}(j) = mean(activesignal_d{i}{j});
    end
    for j = 1:length(activesignal_dl2b{i})
        ave_activesignal_dl2b{i}(j) = mean(activesignal_dl2b{i}{j});
    end
    for j = 1:length(activesignal_dl2d{i})
        ave_activesignal_dl2d{i}(j) = mean(activesignal_dl2d{i}{j});
    end
    for j = 1:length(activesignal_zg{i})
        ave_activesignal_zg{i}(j) = mean(activesignal_zg{i}{j});
    end
    
    
    
    for j = 1:length(activesignal_Rc{i})
        ave_activesignal_Rc{i}(j) = mean(activesignal_Rc{i}{j});
    end
    for j = 1:length(activesignal_Rdz{i})
        ave_activesignal_Rdz{i}(j) = mean(activesignal_Rdz{i}{j});
    end
    for j = 1:length(activesignal_Rd{i})
        ave_activesignal_Rd{i}(j) = mean(activesignal_Rd{i}{j});
    end
end

clear tmp
locations = [0:15]/16;
for i = 1:length(locations)
    tmp{1,i} = ave_activesignal_c{1}(loc_controln1 == locations(i));
%     max_c(1,i) = mean(tmp{i});
    tmp{2,i} = ave_activesignal_c{2}(loc_controln2 == locations(i));
%     max_c(2,i) = mean(tmp{i});
    tmp{3,i} = ave_activesignal_c{3}(loc_controln4 == locations(i));
%     max_c(3,i) = mean(tmp{i});
    tmp{4,i} = ave_activesignal_c{4}(loc_controln3 == locations(i));
    tmp{5,i} = ave_activesignal_c{5}(loc_controln5 == locations(i));
    tmp{6,i} = ave_activesignal_c{6}(loc_controln6 == locations(i));
    
    med_c(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i} tmp{4,i} tmp{5,i} tmp{6,i}]);
    std_med_c(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i} tmp{4,i} tmp{5,i} tmp{6,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i})+length(tmp{4,i})+length(tmp{5,i})+length(tmp{6,i}));
end
control = tmp;
clear tmp
for i = 1:length(locations)
    tmp{1,i} = ave_activesignal_b{1}(loc_bestn1 == locations(i));
%     max_b(1,i) = mean(tmp{i});
    tmp{2,i} = ave_activesignal_b{2}(loc_bestn2 == locations(i));
%     max_b(2,i) = mean(tmp{i});
    tmp{3,i} = ave_activesignal_b{3}(loc_bestn4 == locations(i));
%     max_b(3,i) = mean(tmp{i});
    med_b(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_med_b(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
end
best = tmp;
clear tmp
for i = 1:length(locations)
    tmp{1,i} = ave_activesignal_d{1}(loc_deln1 == locations(i));
%     max_d(1,i) = mean(tmp{i});
    tmp{2,i} = ave_activesignal_d{2}(loc_deln2 == locations(i));
%     max_d(2,i) = mean(tmp{i});
    tmp{3,i} = ave_activesignal_d{3}(loc_deln4 == locations(i));
%     max_d(3,i) = mean(tmp{i});
    med_d(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_med_d(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
end
del = tmp;
clear tmp
for i = 1:length(locations)
    tmp{1,i} = ave_activesignal_dl2b{1}(loc_dl2bestn1 == locations(i));
%     max_dl2b(1,i) = mean(tmp{i});
    tmp{2,i} = ave_activesignal_dl2b{2}(loc_dl2bestn2 == locations(i));
%     max_dl2b(2,i) = mean(tmp{i});
    tmp{3,i} = ave_activesignal_dl2b{3}(loc_dl2bestn5 == locations(i));
%     max_dl2b(3,i) = mean(tmp{i});
    med_dl2b(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_med_dl2b(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
end
dl2best = tmp;
clear tmp
for i = 1:length(locations)
    tmp{1,i} = ave_activesignal_dl2d{1}(loc_dl2deln1 == locations(i));
%     max_dl2d(1,i) = mean(tmp{i});
    tmp{2,i} = ave_activesignal_dl2d{2}(loc_dl2deln2 == locations(i));
%     max_dl2d(2,i) = mean(tmp{i});
    tmp{3,i} = ave_activesignal_dl2d{3}(loc_dl2deln4 == locations(i));
%     max_dl2d(3,i) = mean(tmp{i});
    med_dl2d(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_med_dl2d(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
end
dl2del = tmp;
clear tmp
for i = 1:length(locations)
    tmp{1,i} = ave_activesignal_zg{1}(loc_zldgoodn2 == locations(i));
%     max_zg(1,i) = mean(tmp{i});
    tmp{2,i} = ave_activesignal_zg{2}(loc_zldgoodn3 == locations(i));
%     max_zg(2,i) = mean(tmp{i});
    tmp{3,i} = ave_activesignal_zg{3}(loc_zldgoodn4 == locations(i));
%     max_zg(3,i) = mean(tmp{i});
    med_zg(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_med_zg(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
end
zld = tmp;
clear tmp
for i = 1:length(locations) 
    %Revision
    tmp{1,i} = ave_activesignal_Rc{1}(loc_Rcontroln1 == locations(i));
%     max_c(1,i) = mean(tmp{i});
    tmp{2,i} = ave_activesignal_Rc{2}(loc_Rcontroln2 == locations(i));
%     max_c(2,i) = mean(tmp{i});
    tmp{3,i} = ave_activesignal_Rc{3}(loc_Rcontroln3 == locations(i));
    
    med_Rc(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_med_Rc(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
end
clear tmp
for i = 1:length(locations)
    tmp{1,i} = ave_activesignal_Rdz{1}(loc_Rdzn1 == locations(i));
    tmp{2,i} = ave_activesignal_Rdz{2}(loc_Rdzn2 == locations(i));
    tmp{3,i} = ave_activesignal_Rdz{3}(loc_Rdzn3 == locations(i));
    
    med_Rdz(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_med_Rdz(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));
end
clear tmp
for i = 1:length(locations)
    tmp{1,i} = ave_activesignal_Rd{1}(loc_Rdn1 == locations(i));
    tmp{2,i} = ave_activesignal_Rd{2}(loc_Rdn2 == locations(i));
    tmp{3,i} = ave_activesignal_Rd{3}(loc_Rdn3 == locations(i));
    
    med_Rd(i) = nanmean([tmp{1,i} tmp{2,i} tmp{3,i}]);
    std_med_Rd(i) = nanstd([tmp{1,i} tmp{2,i} tmp{3,i}])./sqrt(length(tmp{1,i})+length(tmp{2,i})+length(tmp{3,i}));

end

figure,errorbar(locations,med_c,std_med_c,'k','LineWidth',2);hold on
errorbar(locations,med_b,std_med_b,'color',rgb('red'),'LineWidth',2);hold on
errorbar(locations,med_d,std_med_d,'--','color',rgb('lightcoral'),'LineWidth',2)
ylim([0 3000])
xlim([0 .9375])
title('Median Active Signal Intensity')
ylabel('Median Intensity')
xlabel('location')
legend('control','Dl1 strong','Dl1 weak','Location','northwest')

figure,errorbar(locations,med_c,std_med_c,'k','LineWidth',2);hold on
errorbar(locations,med_dl2b,std_med_dl2b,'b','LineWidth',2)
errorbar(locations,med_dl2d,std_med_dl2d,'--','color',rgb('cornflowerblue'),'LineWidth',2)
ylim([0 3000])
xlim([0 .9375])
title('Median Active Signal Intensity')
ylabel('Median Intensity')
xlabel('location')
legend('control','Dl2 strong','Dl2 weak','Location','northwest')

figure,errorbar(locations,med_c,std_med_c,'k','LineWidth',2);hold on
errorbar(locations,med_b,std_med_b,'color',rgb('red'),'LineWidth',2);hold on
errorbar(locations,med_dl2b,std_med_dl2b,'b','LineWidth',2)
errorbar(locations,med_zg,std_med_zg,'color',rgb('green'),'LineWidth',2)
ylim([0 3000])
xlim([0 .9375])
title('Median Active Signal Intensity')
ylabel('Median Intensity')
xlabel('location')
legend('control','Dl1 strong','Dl2 strong','Zld strong','Location','northwest')


%Revsion 

figure,errorbar(locations,med_c,std_med_c,'Color',rgb('grey'),'LineWidth',2);hold on
errorbar(locations,med_Rc,std_med_Rc,'Color',rgb('black'),'LineWidth',2);hold on
errorbar(locations,med_Rdz,std_med_Rdz,'Color',rgb('red'),'LineWidth',2);hold on
errorbar(locations,med_Rd,std_med_Rd,'Color',rgb('blue'),'LineWidth',2);hold on

ylim([0 3000])
xlim([0 .9375])
title('Median Active Signal Intensity')
ylabel('Median Intensity')
xlabel('location')
legend('control','Dl1 strong','Dl1 weak','Location','northwest')

%Percent Difference
diffmax(1,:) = (max_b - max_c)./max_c;
diffmax(2,:) = (max_d - max_c)./max_c;
diffmax(3,:) = (max_dl2b - max_c)./max_c;
diffmax(4,:) = (max_dl2d - max_c)./max_c;
diffmax(5,:) = (max_zg - max_c)./max_c;

diffmed(1,:) = (med_b - med_c)./med_c;
diffmed(2,:) = (med_d - med_c)./med_c;
diffmed(3,:) = (med_dl2b - med_c)./med_c;
diffmed(4,:) = (med_dl2d - med_c)./med_c;
diffmed(5,:) = (med_zg - med_c)./med_c;




%% Welch's ttest
for i = 1:16
    [meanb(i),pb(i)] = ttest2([[control{1,i}] [control{2,i}] [control{3,i}] [control{4,i}] [control{5,i}] [control{6,i}]], [[best{1,i}] [best{2,i}] [best{3,i}]],'Vartype','unequal');
    [meand(i),pd(i)] = ttest2([[control{1,i}] [control{2,i}] [control{3,i}] [control{4,i}] [control{5,i}] [control{6,i}]], [[del{1,i}] [del{2,i}] [del{3,i}]],'Vartype','unequal');
    [meandl2b(i),pdl2b(i)] = ttest2([[control{1,i}] [control{2,i}] [control{3,i}] [control{4,i}] [control{5,i}] [control{6,i}]], [[dl2best{1,i}] [dl2best{2,i}] [dl2best{3,i}]],'Vartype','unequal');
    [meandl2d(i),pdl2d(i)] = ttest2([[control{1,i}] [control{2,i}] [control{3,i}] [control{4,i}] [control{5,i}] [control{6,i}]], [[dl2del{1,i}] [dl2del{2,i}] [dl2del{3,i}]],'Vartype','unequal');
    [meanz(i),pz(i)] = ttest2([[control{1,i}] [control{2,i}] [control{3,i}] [control{4,i}] [control{5,i}] [control{6,i}]], [[zld{1,i}] [zld{2,i}] [zld{3,i}]],'Vartype','unequal');   
    

    
%     activateb(i) = ttest2([[actc{1,i}] [actc{2,i}] [actc{3,i}] [actc{4,i}] [actc{5,i}] [actc{6,i}]], [[actb{1,i}] [actb{2,i}] [actb{3,i}]],'Vartype','unequal');
%     activated(i) = ttest2([[actc{1,i}] [actc{2,i}] [actc{3,i}] [actc{4,i}] [actc{5,i}] [actc{6,i}]], [[actd{1,i}] [actd{2,i}] [actd{3,i}]],'Vartype','unequal');
%     activatedl2b(i) = ttest2([[actc{1,i}] [actc{2,i}] [actc{3,i}] [actc{4,i}] [actc{5,i}] [actc{6,i}]], [[actdl2b{1,i}] [actdl2b{2,i}] [actdl2b{3,i}]],'Vartype','unequal');
%     activatedl2d(i) = ttest2([[actc{1,i}] [actc{2,i}] [actc{3,i}] [actc{4,i}] [actc{5,i}] [actc{6,i}]], [[actdl2d{1,i}] [actdl2d{2,i}] [actdl2d{3,i}]],'Vartype','unequal');
%     activatez(i) = ttest2([[actc{1,i}] [actc{2,i}] [actc{3,i}] [actc{4,i}] [actc{5,i}] [actc{6,i}]], [[actz{1,i}] [actz{2,i}] [actz{3,i}]],'Vartype','unequal');   
%     
end
% 
% for i = 3:15
%     outb(i) = ttest2([Area_c{i}], [Area_b{i}],'Vartype','unequal');
%     outdl2b(i) = ttest2([Area_c{i}], [Area_dl2b{i}],'Vartype','unequal');
%     outdl2d(i) = ttest2([Area_c{i}], [Area_dl2d{i}],'Vartype','unequal');
% end
% for i = 3:14
%     outd(i) = ttest2([Area_c{i}], [Area_d{i}],'Vartype','unequal');
% end
% for i = 4:15
%     outz(i) = ttest2([Area_c{i}], [Area_zg{i}],'Vartype','unequal');
% end    

