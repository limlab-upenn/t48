% [disc, bins] = spat_bin_embryo('t48_control', 1, [], 16);
% [outspatial, difftotvecbintrans, bins] = get_HMM_bursts_v2('t48_control', 1, [], 8);
% 
% %statevec = hmmviterbi(seq,TRANS,EMIS)
% ons = [];
% offs = [];
% for i = 1:length(outspatial)
%     temp = outspatial(i).Trans;
%     ons = [ons temp(1, 2)];
%     offs = [offs temp(2, 1)];
% end
% 
% figure
% yyaxis left
% plot(ons)
% yyaxis right
% plot(offs)
% 
% loc = load('G:\Shared drives\Lim_Lab\Sam\Affinity_Project\control_best_deletion_post_analysis.mat','loc_controln1');
% loc = cell2mat(struct2cell(loc))*16;
% 
% figure;
% nucleitomap = randi(length(difftotvecbintrans), 6, 1);
% for i = 1:6
%     subplot(2, 3, i)
%     whichguy = nucleitomap(i);
%     statevec = hmmviterbi(difftotvecbintrans(whichguy, :),outspatial(loc(whichguy)).Trans,outspatial(loc(whichguy)).Emmis);
%     yyaxis left
%     plot(difftotvecbintrans(whichguy, :))
%     yyaxis right
%     plot(statevec)
% end

%%
clear
% close all
% bins = [1237,1797.73333333333,2358.46666666667,2919.20000000000,3479.93333333333,4040.66666666667,4601.40000000000,5162.13333333333,5722.86666666667,6283.60000000000,6844.33333333333,7405.06666666667,7965.80000000000,8526.53333333333,9087.26666666667,9648];
% %sidu
% bins = [0,522.166666666667,1044.33333333333,1566.50000000000,2088.66666666667,2610.83333333333,3133,3655.16666666667,4177.33333333333,4699.50000000000,5221.66666666667,5743.83333333333,6266,6788.16666666667,7310.33333333333,7832.50000000000]; %only control
% bins = [0,643.200000000000,1286.40000000000,1929.60000000000,2572.80000000000,3216,3859.20000000000,4502.40000000000,5145.60000000000,5788.80000000000,6432,7075.20000000000,7718.40000000000,8361.60000000000,9004.80000000000,9648]; %sidu starting at 0
% bins = [984.466591117588,1562.03548504308,2139.60437896858,2717.17327289407,3294.74216681956,3872.31106074506,4449.87995467055,5027.44884859605,5605.01774252154,6182.58663644704,6760.15553037253,7337.72442429802,7915.29331822352,8492.86221214901,9070.43110607451,9648]; %from lowest active threshold value
bins = [984.466591117588,1420.82881837642,1857.19104563524,2293.55327289407,2729.91550015290,3166.27772741173,3602.63995467055,4039.00218192938,4475.36440918821,4911.72663644703,5348.08886370586,5784.45109096469,6220.81331822352,6657.17554548235,7093.53777274117,7529.90000000000]; %from lowest active threshold to smoothed max
% for i = 1
% [~, ~,bins] = get_HMM_bursts_v2('t48_control', i, [], 8);
% end

q = 0;
for i = [1 2 4 3 5 6]
    q = q +1;
    [thingc{q}] = get_HMM_bursts_v4('t48_control', i, bins, 16);
end
q = 0;
for i = [1 2 4]
    q = q +1;
    [thingb{q}] = get_HMM_bursts_v4('t48_dl1_best', i, bins, 16);
    [thingd{q}] = get_HMM_bursts_v4('t48_dl1_deletion', i, bins, 16);
    [thingdl2d{q}] = get_HMM_bursts_v4('t48_dl2_deletion', i, bins, 16);
end
q = 0;
for i = [1 2 5]
    q = q +1;
    [thingdl2b{q}] = get_HMM_bursts_v4('t48_dl2_best', i, bins, 16);
end
q = 0;
for i = [2 3 4]
    q = q +1;
    [thingzg{q}] = get_HMM_bursts_v4('t48_zld_best', i, bins, 16);
end

%Mean and std
for i = 1:16
    for j = 1:6
        Onc(j,i) = [thingc{j}(i).Pon]*.3;
        Offc(j,i) = [thingc{j}(i).Poff]*.3;
    end
    for j = 1:3
        Onb(j,i) = [thingb{j}(i).Pon]*.3;
        Offb(j,i) = [thingb{j}(i).Poff]*.3;
        Ond(j,i) = [thingd{j}(i).Pon]*.3;
        Offd(j,i) = [thingd{j}(i).Poff]*.3;
        Ondl2b(j,i) = [thingdl2b{j}(i).Pon]*.3;
        Offdl2b(j,i) = [thingdl2b{j}(i).Poff]*.3;
        Ondl2d(j,i) = [thingdl2d{j}(i).Pon]*.3;
        Offdl2d(j,i) = [thingdl2d{j}(i).Poff]*.3;
        Onzg(j,i) = [thingzg{j}(i).Pon]*.3;
        Offzg(j,i) = [thingzg{j}(i).Poff]*.3;
    end
end

Ponc = mean(Onc);
stdPonc = std(Onc);
Ponb = mean(Onb);
stdPonb = std(Onb);
Pond = mean(Ond);
stdPond = std(Ond);
Pondl2b = mean(Ondl2b);
stdPondl2b = std(Ondl2b);
Pondl2d = mean(Ondl2d);
stdPondl2d = std(Ondl2d);
Ponzg = mean(Onzg);
stdPonzg = std(Onzg);

Offc(Offc == 0) = NaN;
Offb(Offb == 0) = NaN;
Offd(Offd == 0) = NaN;
Offdl2b(Offdl2b == 0) = NaN;
Offdl2d(Offdl2d == 0) = NaN;
Offzg(Offzg == 0) = NaN;

Poffc = nanmean(1./Offc);
stdPoffc = nanstd(1./Offc);
Poffb = nanmean(1./Offb);
stdPoffb = nanstd(1./Offb);
Poffd = nanmean(1./Offd);
stdPoffd = nanstd(1./Offd);
Poffdl2b = nanmean(1./Offdl2b);
stdPoffdl2b = nanstd(1./Offdl2b);
Poffdl2d = nanmean(1./Offdl2d);
stdPoffdl2d = nanstd(1./Offdl2d);
Poffzg = nanmean(1./Offzg);
stdPoffzg = nanstd(1./Offzg);


figure; errorbar(([1:16]-1)./16,Ponc,stdPonc,'Color',rgb('black'),'LineWidth',2); hold on
errorbar(([1:16]-1)./16,Ponb,stdPonb,'Color',rgb('red'),'LineWidth',2);
errorbar(([1:16]-1)./16,Pondl2b,stdPondl2b,'Color',rgb('blue'),'LineWidth',2);
errorbar(([1:16]-1)./16,Ponzg,stdPonzg,'Color',rgb('green'),'LineWidth',2);
ylim([0 .07])
xlim([0 .9375])
legend('Control','Dl1 optimized','Dl2 optimized','Location','Northwest')
title('Kon')

figure; errorbar(([1:16]-1)./16,Poffc,stdPoffc,'Color',rgb('black'),'LineWidth',2); hold on
errorbar(([1:16]-1)./16,Poffb,stdPoffb,'Color',rgb('red'),'LineWidth',2);
errorbar(([1:16]-1)./16,Poffdl2b,stdPoffdl2b,'Color',rgb('blue'),'LineWidth',2);
errorbar(([1:16]-1)./16,Poffzg,stdPoffzg,'Color',rgb('green'),'LineWidth',2);
ylim([0 .4])
xlim([0 .9375])
legend('Control','Dl1 optimized','Dl2 optimized','Location','Northwest')
title('Koff')

figure; errorbar(([1:16]-1)./16,Ponc,stdPonc,'Color',rgb('black'),'LineWidth',2); hold on
errorbar(([1:16]-1)./16,Pond,stdPond,'--','Color',rgb('lightcoral'),'LineWidth',2);
errorbar(([1:16]-1)./16,Pondl2d,stdPondl2d,'--','Color',rgb('cornflowerblue'),'LineWidth',2);
ylim([0 .07])
xlim([0 .9375])
legend('Control','Dl1 deletion','Dl2 deletion','Location','Northwest')
title('Kon')

figure; errorbar(([1:16]-1)./16,Poffc,stdPoffc,'Color',rgb('black'),'LineWidth',2); hold on
errorbar(([1:16]-1)./16,Poffd,stdPoffd,'--','Color',rgb('lightcoral'),'LineWidth',2);
errorbar(([1:16]-1)./16,Poffdl2d,stdPoffdl2d,'--','Color',rgb('cornflowerblue'),'LineWidth',2);
ylim([0 .4])
xlim([0 .9375])
legend('Control','Dl1 deletion','Dl2 deletion','Location','Northwest')
title('Koff')


%Burst frequency
for i = 1:16
    for j = 1:6
        Freqc(j,i) = Onc(j,i)*Offc(j,i)/(Onc(j,i)+Offc(j,i));
    end
    for j = 1:3
        Freqb(j,i) = Onb(j,i)*Offb(j,i)/(Onb(j,i)+Offb(j,i));
        Freqd(j,i) = Ond(j,i)*Offd(j,i)/(Ond(j,i)+Offd(j,i));
        Freqdl2b(j,i) = Ondl2b(j,i)*Offdl2b(j,i)/(Ondl2b(j,i)+Offdl2b(j,i));
        Freqdl2d(j,i) = Ondl2d(j,i)*Offdl2d(j,i)/(Ondl2d(j,i)+Offdl2d(j,i));
        Freqzg(j,i) = Onzg(j,i)*Offzg(j,i)/(Onzg(j,i)+Offzg(j,i));
    end
end

for i = 1:16
    for j = 1:6
        Occ(j,i) = Onc(j,i)/(Onc(j,i)+Offc(j,i));
    end
    for j = 1:3
        Ocb(j,i) = Onb(j,i)/(Onb(j,i)+Offb(j,i));
        Ocd(j,i) = Ond(j,i)/(Ond(j,i)+Offd(j,i));
        Ocdl2b(j,i) = Ondl2b(j,i)/(Ondl2b(j,i)+Offdl2b(j,i));
        Ocdl2d(j,i) = Ondl2d(j,i)/(Ondl2d(j,i)+Offdl2d(j,i));
        Oczg(j,i) = Onzg(j,i)/(Onzg(j,i)+Offzg(j,i));
    end
end

Fc = nanmean(Freqc);
stdFc = nanstd(Freqc);
Fb = nanmean(Freqb);
stdFb = nanstd(Freqb);
Fd = nanmean(Freqd);
stdFd = nanstd(Freqd);
Fdl2b = nanmean(Freqdl2b);
stdFdl2b = nanstd(Freqdl2b);
Fdl2d = nanmean(Freqdl2d);
stdFdl2d = nanstd(Freqdl2d);
Fzg = nanmean(Freqzg);
stdFzg = nanstd(Freqzg);

Oc = mean(Occ);
stdOc = std(Occ);
Ob = mean(Ocb);
stdOb = std(Ocb);
Od = mean(Ocd);
stdOd = std(Ocd);
Odl2b = mean(Ocdl2b);
stdOdl2b = std(Ocdl2b);
Odl2d = mean(Ocdl2d);
stdOdl2d = std(Ocdl2d);
Ozg = mean(Oczg);
stdOzg = std(Oczg);

figure,errorbar(([1:16]-1)./16,Fc,stdFc,'Color',rgb('black'),'LineWidth',2); hold on
errorbar(([1:16]-1)./16,Fb,stdFb,'Color',rgb('red'),'LineWidth',2); 
errorbar(([1:16]-1)./16,Fdl2b,stdFdl2b,'Color',rgb('blue'),'LineWidth',2); 
errorbar(([1:16]-1)./16,Fzg,stdFzg,'Color',rgb('green'),'LineWidth',2); 
ylim([0 .05])
xlim([0 .9375])
legend('Control','Dl1 Opt','Dl2 Opt','Zld Opt','Location','Northwest')
title('Burst Frequency')

figure,errorbar(([1:16]-1)./16,Fc,stdFc,'Color',rgb('black'),'LineWidth',2); hold on
errorbar(([1:16]-1)./16,Fd,stdFd,'--','Color',rgb('lightcoral'),'LineWidth',2); 
errorbar(([1:16]-1)./16,Fdl2d,stdFdl2d,'--','Color',rgb('cornflowerblue'),'LineWidth',2); 
ylim([0 .05])
xlim([0 .9375])
legend('Control','Dl1 Del','Dl2 Del','Location','Northwest')
title('Burst Frequency')

figure,errorbar(([1:16]-1)./16,Oc,stdOc,'Color',rgb('black'),'LineWidth',2); hold on
errorbar(([1:16]-1)./16,Ob,stdOb,'Color',rgb('red'),'LineWidth',2); 
errorbar(([1:16]-1)./16,Odl2b,stdOdl2b,'Color',rgb('blue'),'LineWidth',2); 
errorbar(([1:16]-1)./16,Ozg,stdOzg,'Color',rgb('green'),'LineWidth',2); 
ylim([0 .25])
xlim([0 .9375])
legend('Control','Dl1 Opt','Dl2 Opt','Zld Opt','Location','Northwest')
title('Occupancy')

figure,errorbar(([1:16]-1)./16,Oc,stdOc,'Color',rgb('black'),'LineWidth',2); hold on
errorbar(([1:16]-1)./16,Od,stdOd,'--','Color',rgb('lightcoral'),'LineWidth',2); 
errorbar(([1:16]-1)./16,Odl2d,stdOdl2d,'--','Color',rgb('cornflowerblue'),'LineWidth',2); 
ylim([0 .25])
xlim([0 .9375])
legend('Control','Dl1 Del','Dl2 Del','Location','Northwest')
title('Occupancy')
%% Kini
for i = 1:6
    for j = 1:16
        binstemp = mean(thingc{i}(j).bintrajs);
        binstemp1 = diff(binstemp');
        binstemp1 = binstemp1';
        binstemp1 = [0 binstemp1];
        
        statestemp = mean(thingc{i}(j).statesmat-1);
        
        kinic{i,j}(:) = binstemp1./statestemp;
        clear binstemp binstemp1 statestemp
    end
end

for i = 1:3
    for j = 1:16
        binstemp = mean(thingb{i}(j).bintrajs);
        binstemp1 = diff(binstemp');
        binstemp1 = binstemp1';
        binstemp1 = [0 binstemp1];
        
        statestemp = mean(thingb{i}(j).statesmat-1);
        
        kinib{i,j}(:) = binstemp1./statestemp;
        clear binstemp binstemp1 statestemp
    end
end
for i = 1:3
    for j = 1:16
        binstemp = mean(thingd{i}(j).bintrajs);
        binstemp1 = diff(binstemp');
        binstemp1 = binstemp1';
        binstemp1 = [0 binstemp1];
        
        statestemp = mean(thingd{i}(j).statesmat-1);
        
        kinid{i,j}(:) = binstemp1./statestemp;
        clear binstemp binstemp1 statestemp
    end
end
for i = 1:3
    for j = 1:16
        binstemp = mean(thingdl2b{i}(j).bintrajs);
        binstemp1 = diff(binstemp');
        binstemp1 = binstemp1';
        binstemp1 = [0 binstemp1];
        
        statestemp = mean(thingdl2b{i}(j).statesmat-1);
        
        kinidl2b{i,j}(:) = binstemp1./statestemp;
        clear binstemp binstemp1 statestemp
    end
end
for i = 1:3
    for j = 1:16
        binstemp = mean(thingdl2d{i}(j).bintrajs);
        binstemp1 = diff(binstemp');
        binstemp1 = binstemp1';
        binstemp1 = [0 binstemp1];
        
        statestemp = mean(thingdl2d{i}(j).statesmat-1);
        
        kinidl2d{i,j}(:) = binstemp1./statestemp;
        clear binstemp binstemp1 statestemp
    end
end
for i = 1:3
    for j = 1:16
        binstemp = mean(thingzg{i}(j).bintrajs);
        binstemp1 = diff(binstemp');
        binstemp1 = binstemp1';
        binstemp1 = [0 binstemp1];
        
        statestemp = mean(thingzg{i}(j).statesmat-1);
        
        kinizg{i,j}(:) = binstemp1./statestemp;
        clear binstemp binstemp1 statestemp
    end
end
for i = 1:3
    subplot(3,6,1+6*(i-1))
    plot(kinic{i,8});
    ylim([-15 5])
    subplot(3,6,2+6*(i-1))
    plot(kinib{i,8});
    ylim([-15 5])
    subplot(3,6,3+6*(i-1))
    plot(kinidl2b{i,8});
    ylim([-15 5])
    subplot(3,6,4+6*(i-1))
    plot(kinid{i,8});ylim([-15 5])
    subplot(3,6,5+6*(i-1))
    plot(kinidl2d{i,8});ylim([-15 5])
    subplot(3,6,6+6*(i-1))
    plot(kinizg{i,8});ylim([-15 5])
end

subplot(3,6,1)
title('Control')
subplot(3,6,2)
title('Dl1 Opt')
subplot(3,6,3)
title('Dl2 Opt')
subplot(3,6,4)
title('Dl1 Del')
subplot(3,6,5)
title('Dl2 Del')
subplot(3,6,6)
title('Zld Opt')

load('G:\Shared drives\Lim_Lab\Sam\Affinity_Project\control_best_deletion_post_analysis.mat')

%% Promoter State

m = [7 7 7 7 8 8 8 8 9 9 9 9 9 10 10 10 11 11 11 11]; %DV bin number
n = 1:2:40; %nuclei in that bin
for p = 1:20
i = m(p);
j = n(p)+3;

figure;
subplot(2,1,1)
plot(time_zldgoodn4,smooth(M_zldgoodn4(:,izg4{i}(j)),4),'k','LineWidth',2)
xlim([0 1])
subplot(2,1,2)
plot(time_zldgoodn4,[0 thingzg{1,3}(i).promoter(j,:)],'k','LineWidth',2)
xlim([0 1])
ylim([0 1.25])
end