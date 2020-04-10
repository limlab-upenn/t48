clear;
[filename, pathname]= uigetfile({'*.mat','Data files(*.mat)'},'Choose');
load([pathname filename]);

%Smoothing
%{
% smovingM=smooth(M_raw(:,100));
% slowessM=smooth(M_raw(:,100),'lowess');
% subplot(1,2,1), plot(T,smovingM,'b'); hold on;
% title('Standard Smooth','FontSize', 16);
% ylabel('Pixel Intensity','FontWeight','bold');
% xlabel('Time (min)','FontWeight','bold');
% subplot(1,2,2), plot(T,slowessM,'b'); hold on;
% title('lowess Smooth','FontSize', 16);
% ylabel('Pixel Intensity','FontWeight','bold');
% xlabel('Time (min)','FontWeight','bold');
% legend('Original', 'lowess smoothing');
% subplot(1,2,1), plot(T,M_raw(:,100),'r');
% legend('Original', 'Normal smoothing');
% subplot(1,2,2), plot(T,M_raw(:,100),'r');
% legend('Original', 'lowess smoothing');
%}
%Normalize by subtracting M by Mm and then by the smallest value in M
M_raw = M;
Mm = double(Mm);


for i=1:size(M,2)
    M(:,i) = M(:,i)-Mm(:,i);
    M(:,i) = M(:,i)-min(M(:,i));
    
    %sM(:,i) = smooth(M1(:,i));
end

%active nuclei
for i = 1:size(M,2)
    check(i) = max(M(:,i));
    tmp(:,i) = sort(M(:,i));   
    maxM = tmp(end-1,:);
end

thresh = .8;
thresh_base = mean(mean(Mm));
activeM = find(maxM > thresh*thresh_base);
    
%time of activation
for j = 1:length(activeM)
    T_M = find(M(1:end,activeM(j)) > thresh*thresh_base);
%     T_activation(j) = T(T_M(1)-1);  %many cases activation starts time point before due to steep slope
    T_activation(j) = T(T_M(2));   %use second time point to ensure nuclei is active
end
aveT_act = mean(T_activation);
stdT_act = std(T_activation);

%Number of nuclei active vs time
nuclei_by_T = 0;
for i = 1:length(T)
    nuclei_active = find(T_activation == T(i));
    if i ==1
        nuclei_by_T(i) = nuclei_by_T + size(nuclei_active,2);
    else
        nuclei_by_T(i) = nuclei_by_T(i-1) + size(nuclei_active,2);
    end
end

T_start = T_activation./T(end);
Max = max(M);

figure;
subplot(4,5,5)
plot(T,nuclei_by_T/length(activeM))
title('Number Active Nuclei vs Time')
ylabel('Number of Active Nuclei')
xlabel('Time (min)')
xlim([0 60])

%area under curve
% for j = 1:length(activeM)
%     Area(j) = trapz(T,M(:,activeM(j)));
% end

 %maximum amplitude
for i = 1:length(activeM)
    maxAmp(i) = max(M(:,activeM(i)));
    aveAmp = mean(maxAmp);
    stdAmp = std(maxAmp);
end
   
% figure, histogram(T_activation,20);
% title('Activation Time');
% ylabel('Frequency');
% xlabel('Time (min)');
% figure, histogram(maxAmp,20);
% title('Maximum Amplitude per Active Nucleus');
% ylabel('Frequency');
% xlabel('Pixel Intensity');
% figure, histogram(Area,16);   
% title('Area of Active Nuclei Signal');
% ylabel('Frequency');
% xlabel('Area Under Signal Curve');

% Max_zw{3} = check;
% clearvars -except Max_zw
%%
% figure,
k=74;
j = 1;
for i = 1:20
    j = j+18;
    subplot(4,5,i)
    plot(T,M_raw(:,j)); hold on
%     plot(T,.8*thresh_base*ones(size(M(:,k))));
    title('Nuclei Signal');
    ylim([0 10000])
    xlabel('Time (min)');
    ylabel('Pixel Intensity');
%     legend('Pixel Intensity','Activatin threshold');
end
%%
% plot(T,M_raw(:,activeM(138))); hold on
% plot(T,Mm(activeM(138))*ones(size(M(:,100))));
% title('Median and Nuclei Signal');
% ylabel('Pixel Intensity');
% xlabel('Time (min)');
% legend('Signal','Median');


%% Image Bins
% imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\t48_control\t48controln1\';
% I1 = imread([imagepath 'segmentation\070.tif']);
% height = 490; %pixel height of embryo
% time = 70;
% x = [3;248];
% y = [236;510];
% cond = 1;

% imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\t48_control\t48controln2\';
% I1 = imread([imagepath 'segmentation\053.tif']);
% height = 535; %pixel height of embryo
% time = 53;
% x = [3;222];
% y = [262;507];
% cond = 2;

% imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\t48_control\t48controln4\';
% I1 = imread([imagepath 'segmentation\051.tif']);
% height = 530; %pixel height of embryo
% time = 51;
% x = [3;223];
% y = [285;508];
% cond = 7;

% imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\t48_control\t48controln3\';
% I1 = imread([imagepath 'segmentation\059.tif']);
% height = 530; %pixel height of embryo
% time = 59;
% x = [3;238];
% y = [239;509];
% cond = 46;

% imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\t48_control\t48controln5\';
% I1 = imread([imagepath 'segmentation\057.tif']);
% height = 470; %pixel height of embryo
% time = 57;
% x = [3;293];
% y = [194;508];
% cond = 47;

% imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\t48_control\t48controln6\';
% I1 = imread([imagepath 'segmentation\055.tif']);
% height = 500; %pixel height of embryo
% time = 55;
% x = [4;259];
% y = [225;509];
% cond = 48;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_dl1_deletion\t48dl1deln1\';
% I1 = imread([imagepath 'segmentation\053.tif']);
% height = 530;
% time = 53;
% x = [4;240];
% y = [257;510];
% cond = 3;

% imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\t48_dl1_deletion\t48dl1deln2\';
% I1 = imread([imagepath 'segmentation\054.tif']);
% height = 555;
% time = 54;
% x = [5;221];
% y = [283;509];
% cond = 4;

% imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\t48_dl1_deletion\t48dl1deln4\';
% I1 = imread([imagepath 'segmentation\054.tif']);
% height = 420;
% time = 54;
% x = [2;285];
% y = [206;509];
% cond = 8;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_dl1_best\t48dl1bestn1\';
% I1 = imread([imagepath 'segmentation\058.tif']);
% height = 460;
% time = 58;
% x = [6;259];
% y = [236;507];
% cond = 5;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_dl1_best\t48dl1bestn2\';
% I1 = imread([imagepath 'segmentation\057.tif']);
% height = 530;
% time = 57;
% x = [4;229];
% y = [263;507];
% cond = 6;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_dl1_best\t48dl1bestn4\';
% I1 = imread([imagepath 'segmentation\057.tif']);
% height = 530;
% time = 57;
% x = [2;222];
% y = [268;510];
% cond = 9;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_dl2_best\t48dl2bestn1\';
% I1 = imread([imagepath 'segmentation\058.tif']);
% height = 430;
% time = 58;
% x = [2;274];
% y = [229;506];
% cond = 10;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_dl2_best\t48dl2bestn2\';
% I1 = imread([imagepath 'segmentation\056.tif']);
% height = 520;
% time = 56;
% x = [3;248];
% y = [231;510];
% cond = 11;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_dl2_best\t48dl2bestn5\';
% I1 = imread([imagepath 'segmentation\050.tif']);
% height = 500;
% time = 50;
% x = [2;243];
% y = [237;502];
% cond = 12;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_dl2_deletion\t48dl2deln1\';
% I1 = imread([imagepath 'segmentation\052.tif']);
% height = 530;
% time = 52;
% x = [3;227];
% y = [233;509];
% cond = 13;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_dl2_deletion\t48dl2deln2\';
% I1 = imread([imagepath 'segmentation\055.tif']);
% height = 445;
% time = 55;
% x = [3;305];
% y = [214;507];
% cond = 14;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_dl2_deletion\t48dl2deln4\';
% I1 = imread([imagepath 'segmentation\049.tif']);
% height = 500;
% time = 49;
% x = [3;255];
% y = [228;507];
% cond = 15;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_zld_best\t48zldgoodn2\';
% I1 = imread([imagepath 'segmentation\053.tif']);
% height = 510;
% time = 53;
% x = [1;266];
% y = [190;509];
% cond = 16;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_zld_best\t48zldgoodn3\';
% I1 = imread([imagepath 'segmentation\054.tif']);
% height = 465;
% time = 54;
% x = [2;275];
% y = [218;506];
% cond = 17;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_zld_best\t48zldgoodn4\';
% I1 = imread([imagepath 'segmentation\052.tif']);
% height = 545;
% time = 52;
% x = [5;215];
% y = [268;513];
% cond = 18;

% imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\t48_twi4_deletion\twi4deln1\';
% I1 = imread([imagepath 'segmentation\042.tif']);
% height = 545;
% time = 42;
% x = [3;277];
% y = [241;509];
% cond = 49;

% imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\t48_twi4_deletion\twi4deln2\';
% I1 = imread([imagepath 'segmentation\059.tif']);
% height = 545;
% time = 42;
% x = [3;199];
% y = [292;502];
% cond = 50;

%Revision

% imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\Revision\t48_control\t48Rcontroln1\';
% I1 = imread([imagepath 'segmentation\049.tif']);
% height = 510; %pixel height of embryo
% time = 49;
% x = [4;252];
% y = [261;507];
% cond = 51;

% imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\Revision\t48_control\t48Rcontroln2\';
% I1 = imread([imagepath 'segmentation\043.tif']);
% height = 550; %pixel height of embryo
% time = 43;
% x = [5;219];
% y = [302;506];
% cond = 52;

% imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\Revision\t48_control\t48Rcontroln3\';
% I1 = imread([imagepath 'segmentation\045.tif']);
% height = 590; %pixel height of embryo
% time = 45;
% x = [4;186];
% y = [341;503];
% cond = 53;

% imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\Revision\t48_Dl_del_Zld_opt\t48Rdldelzldoptn1\';
% I1 = imread([imagepath 'segmentation\051.tif']);
% height = 565; %pixel height of embryo
% time = 51;
% x = [7;210];
% y = [302;507];
% cond = 54;

% imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\Revision\t48_Dl_del_Zld_opt\t48Rdldelzldoptn2\';
% I1 = imread([imagepath 'segmentation\051.tif']);
% height = 590; %pixel height of embryo
% time = 51;
% x = [5;192];
% y = [312;510];
% cond = 55;

% imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\Revision\t48_Dl_del_Zld_opt\t48Rdldelzldoptn3\';
% I1 = imread([imagepath 'segmentation\056.tif']);
% height = 550; %pixel height of embryo
% time = 56;
% x = [5;218];
% y = [307;509];
% cond = 56;

% imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\Revision\t48_Dl_del\t48Rdldeln1\';
% I1 = imread([imagepath 'segmentation\050.tif']);
% height = 570; %pixel height of embryo
% time = 50;
% x = [2;195];
% y = [292;508];
% cond = 57;

% imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\Revision\t48_Dl_del\t48Rdldeln2\';
% I1 = imread([imagepath 'segmentation\050.tif']);
% height = 610; %pixel height of embryo
% time = 50;
% x = [5;179];
% y = [309;505];
% cond = 58;

imagepath = 'G:\Shared drives\Lim_Lab\Sam\Affinity_Project\Revision\t48_Dl_del\t48Rdldeln3\';
I1 = imread([imagepath 'segmentation\050.tif']);
height = 550; %pixel height of embryo
time = 50;
x = [1;202];
y = [284;508];
cond = 59;
%% Dorsal het

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\MCPhet_t48_control\mcphett48controln1\';
% I1 = imread([imagepath 'segmentation\050.tif']);
% height = 500;
% time = 50;
% x = [0;281];
% y = [211;512];
% cond = 19;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\MCPhet_t48_control\mcphett48controln2\';
% I1 = imread([imagepath 'segmentation\051.tif']);
% height = 580;
% time = 51;
% x = [1;215];
% y = [261;511];
% cond = 20;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\MCPhet_t48_control\mcphett48controln3\';
% I1 = imread([imagepath 'segmentation\051.tif']);
% height = 545;
% time = 51;
% x = [5;215];
% y = [268;513];
% cond = 21;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\Dlhet_t48_control\dlhett48controln1\';
% I1 = imread([imagepath 'segmentation\059.tif']);
% height = 510;
% time = 59;
% x = [4;255];
% y = [240;509];
% cond = 22;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\Dlhet_t48_control\dlhett48controln3\';
% I1 = imread([imagepath 'segmentation\056.tif']);
% height = 540;
% time = 56;
% x = [4;241];
% y = [258;507];
% cond = 23;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\Dlhet_t48_control\dlhett48controln4\';
% I1 = imread([imagepath 'segmentation\056.tif']);
% height = 580;
% time = 56;
% x = [1;227];
% y = [268;509];
% cond = 24;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\Dlhet_t48_dl2best\dlhett48dl2bestn2\';
% I1 = imread([imagepath 'segmentation\052.tif']);
% height = 570;
% time = 52;
% x = [4;220];
% y = [237;509];
% cond = 25;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\Dlhet_t48_dl2best\dlhett48dl2bestn3\';
% I1 = imread([imagepath 'segmentation\054.tif']);
% height = 500;
% time = 54;
% x = [5;268];
% y = [207;512];
% cond = 26;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\Dlhet_t48_dl2best\dlhett48dl2bestn5\';
% I1 = imread([imagepath 'segmentation\054.tif']);
% height = 530;
% time = 54;
% x = [3;235];
% y = [240;508];
% cond = 27;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_zld_weak\t48zldcagn1\';
% I1 = imread([imagepath 'segmentation\063.tif']);
% height = 510;
% time = 63;
% x = [4;249];
% y = [245;508];
% cond = 28;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_zld_weak\t48zldcagn2\';
% I1 = imread([imagepath 'segmentation\055.tif']);
% height = 495;
% time = 55;
% x = [7;257];
% y = [257;511];
% cond = 29;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_zld_weak\t48zldcagn5\';
% I1 = imread([imagepath 'segmentation\051.tif']);
% height = 560;
% time = 51;
% x = [5;214];
% y = [269;505];
% cond = 30;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48singleMCP\t48singleMCPn1\';
% I1 = imread([imagepath 'segmentation\051.tif']);
% height = 400;
% time = 79;
%% Het

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_het_control\t48hetcontroln1\';
% I1 = imread([imagepath 'segmentation\084.tif']);
% height = 560;
% time = 84;
% x = [3;215];
% y = [297;507];
% cond = 31;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_het_control\t48hetcontroln2\';
% I1 = imread([imagepath 'segmentation\081.tif']);
% height = 620;
% time = 81;
% x = [9;196];
% y = [329;507];
% cond = 32;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_het_control\t48hetcontroln3\';
% I1 = imread([imagepath 'segmentation\087.tif']);
% height = 650;
% time = 87;
% x = [9;210];
% y = [288;509];
% cond = 33;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_het_dl1_best\t48hetdl1bestn1\';
% I1 = imread([imagepath 'segmentation\093.tif']);
% height = 560;
% time = 93;
% x = [9;239];
% y = [290;507];
% cond = 34;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_het_dl1_best\t48hetdl1bestn2\';
% I1 = imread([imagepath 'segmentation\090.tif']);
% height = 550;
% time = 90;
% x = [3;216];
% y = [296;508];
% cond = 35;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_het_dl1_best\t48hetdl1bestn3\';
% I1 = imread([imagepath 'segmentation\076.tif']);
% height = 670;
% time = 76;
% x = [2;191];
% y = [325;510];
% cond = 36;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_het_dl2_best\t48hetdl2bestn1\';
% I1 = imread([imagepath 'segmentation\080.tif']);
% height = 550;
% time = 80;
% x = [3;222];
% y = [271;510];
% cond = 37;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_het_dl2_best\t48hetdl2bestn2\';
% I1 = imread([imagepath 'segmentation\077.tif']);
% height = 610;
% time = 77;
% x = [4;208];
% y = [284;506];
% cond = 38;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_het_dl2_best\t48hetdl2bestn3\';
% I1 = imread([imagepath 'segmentation\082.tif']);
% height = 430;
% time = 82;
% x = [2;295];
% y = [233;508];
% cond = 39;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_het_dl2_deletion\t48hetdl2deln1\';
% I1 = imread([imagepath 'segmentation\086.tif']);
% height = 440;
% time = 86;
% x = [1;311];
% y = [199;510];
% cond = 40;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_het_dl2_deletion\t48hetdl2deln2\';
% I1 = imread([imagepath 'segmentation\077.tif']);
% height = 630;
% time = 77;
% x = [4;182];
% y = [331;507];
% cond = 41;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_het_dl2_deletion\t48hetdl2deln3\';
% I1 = imread([imagepath 'segmentation\078.tif']);
% height = 610;
% time = 78;
% x = [5;231];
% y = [277;509];
% cond = 42;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_het_zld_best\t48hetzldbestn1\';
% I1 = imread([imagepath 'segmentation\081.tif']);
% height = 625;
% time = 81;
% x = [5;214];
% y = [286;507];
% cond = 43;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_het_zld_best\t48hetzldbestn2\';
% I1 = imread([imagepath 'segmentation\076.tif']);
% height = 400;
% time = 76;
% x = [6;339];
% y = [205;510];
% cond = 44;

% imagepath = 'G:\Team Drives\Lim_Lab\Sam\Affinity_Project\t48_het_zld_best\t48hetzldbestn3\';
% I1 = imread([imagepath 'segmentation\082.tif']);
% height = 695;
% time = 82;
% x = [1;163];
% y = [364;508];
% cond = 45;

%% Un-comment when doing new run
% figure(2); imshow(I1); hold on; 
% [x,y] = ginput(2); % pick two points from the image

%%
p1 = polyfit(x,y,1); % get a linear fit to the points

plot(x,y,'wo-');
yint1 = y(1)-p1(1)*x(1);


bins = 16;
 
signalT = 1:length(T);%[find(T == 9.998999999999999) find(T == 24.997500000000000) find(T == 40.995899999999990)];

for i = 1:bins
    if i ~= 1
        yint2(i) = yint2(i-1) - height/bins; %subtract since y-axis is flipped
        yint1(i) = yint1(i-1) - height/bins;
    else
        yint2(i) = yint1(i) - height/bins;
    end
    
    y1 = p1(1)*lineage_cx(time,:)+yint1(i);
    y2 = p1(1)*lineage_cx(time,:)+yint2(i);

    
    bin{i} = find(y2 < lineage_cy(time,:) & y1 > lineage_cy(time,:));
    liney1(i,:) = p1(1).*[0:512]+yint1(i);
    liney2(i,:) = p1(1).*[0:512]+yint2(i);
  
end

%Segement into A-P regions
y_shift = 241*sqrt(2);

y_AP1 = -1*lineage_cx(time,:)+yint1(1)+100;
y_AP2 = -1*lineage_cx(time,:)+yint1(1)+100+y_shift;
bin_A = find(lineage_cy(time,:) < y_AP1);
bin_mid = find(lineage_cy(time,:) > y_AP1 & lineage_cy(time,:) < y_AP2);
bin_P = find(lineage_cy(time,:) > y_AP2);

liney3 = -1.*[0:512]+yint1(1)+100;
liney4 = -1.*[0:512]+yint1(1)+100+y_shift;

plot([0:512],liney3,'w',[0:512],liney4,'w')
plot([0:512],liney1,'w',[0:512],liney2,'w',[0:512],liney3,'w',[0:512],liney4,'r'); hold off

%AP bins for output
for i = 1:length(bin)
    for j = 1:length(intersect(bin{1,i},bin_A))
        tmp = intersect(bin{1,i},bin_A);
        Area_bin_A{1,i}(j) = trapz(T,M(:,tmp(j)));
    end
    for j = 1:length(intersect(bin{1,i},bin_mid))
        tmp = intersect(bin{1,i},bin_mid);
        Area_bin_mid{1,i}(j) = trapz(T,M(:,tmp(j)));
    end
    for j = 1:length(intersect(bin{1,i},bin_P))
        tmp = intersect(bin{1,i},bin_P);
        Area_bin_P{1,i}(j) = trapz(T,M(:,tmp(j)));
    end
end


%removes all inactive nuclei 
for k=1:bins
    active_bin{1,k}=intersect(bin{1,k},activeM);
end


%bin output/signal
for k = 2:length(signalT)
    for i = 1:size(bin,2)
        for j = 1:size(bin{1,i},2)
             Area_bin{k-1,i}(j) = trapz(T(1:k),M((1:k),bin{1,i}(j)));
        end
    end
end
for k = 1:length(signalT)
    for i = 1:size(bin,2)
        for j = 1:size(bin{1,i},2)
%             Area_bin{k,i}(j) = trapz(T(k),M((k),bin{1,i}(j)));
            signal_bin{k,i}(j) = M(signalT(k),bin{1,i}(j));
        end
        for q = 1:length(intersect(bin{1,i},bin_A))
            tmp = intersect(bin{1,i},bin_A);
            signal_bin_A{k,i}(q) = M(signalT(k),tmp(q));
        end
        for q = 1:length(intersect(bin{1,i},bin_mid))
            tmp = intersect(bin{1,i},bin_mid);
            signal_bin_mid{k,i}(q) = M(signalT(k),tmp(q));
        end
        for q = 1:length(intersect(bin{1,i},bin_P))
            tmp = intersect(bin{1,i},bin_P);
            signal_bin_P{k,i}(q) = M(signalT(k),tmp(q));
        end
    end
end
%fixing empty indices
if size(Area_bin,2) ~= length(bin)
    for p = 1:length(Area_bin)
        for k = size(Area_bin,2)+1:length(bin)
            Area_bin{p,k} = 0;
        end
    end
end
if size(signal_bin,2) ~= length(bin)
    for i = 1:length(signalT)
        for k = size(signal_bin,2)+1:length(bin)
            signal_bin{i,k} = 0;
        end
    end
end
if size(Area_bin_A,2) ~= length(bin)
    for k = size(Area_bin_A,2)+1:length(bin)
        Area_bin_A{1,k} = 0;
    end
end
if size(Area_bin_mid,2) ~= length(bin)
    for k = size(Area_bin_mid,2)+1:length(bin)
        Area_bin_mid{1,k} = 0;
    end
end
if size(Area_bin_P,2) ~= length(bin)
    for k = size(Area_bin_P,2)+1:length(bin)
        Area_bin_P{1,k} = 0;
    end
end
if size(signal_bin_A,2) ~= length(bin)
    for i = 1:length(signalT)
        for k = size(signal_bin_A,2)+1:length(bin)
            signal_bin_A{i,k} = 0;
        end
    end
end
if size(signal_bin_mid,2) ~= length(bin)
    for i = 1:length(signalT)
        for k = size(signal_bin_mid,2)+1:length(bin)
            signal_bin_mid{i,k} = 0;
        end
    end
end
if size(signal_bin_P,2) ~= length(bin)
    for i = 1:length(signalT)
        for k = size(signal_bin_P,2)+1:length(bin)
            signal_bin_P{i,k} = 0;
        end
    end
end

%averaging output/signal
for k = 1:length(Area_bin)
    for i = 1:length(bin)
         aveArea_bin(k,i) = mean(Area_bin{k,i});
         stdArea_bin(k,i) = std(Area_bin{k,i});%/sqrt(length(Area_bin{k,i}));
    end
end
for k = 1:length(signalT)
    for i = 1:length(bin)
%         aveArea_bin(i) = mean(Area_bin{1,i});
        aveArea_bin_A(i) = mean(Area_bin_A{1,i});
        aveArea_bin_mid(i) = mean(Area_bin_mid{1,i});
        aveArea_bin_P(i) = mean(Area_bin_P{1,i});
%         stdArea_bin(i) = std(Area_bin{1,i})/sqrt(length(Area_bin{1,i}));
        stdArea_bin_A(i) = std(Area_bin_A{1,i})/sqrt(length(Area_bin_A{1,i}));
        stdArea_bin_mid(i) = std(Area_bin_mid{1,i})/sqrt(length(Area_bin_mid{1,i}));
        stdArea_bin_P(i) = std(Area_bin_P{1,i})/sqrt(length(Area_bin_P{1,i}));
        tmp1(i) = mean(signal_bin{k,i});
        tmp2(i) = mean(signal_bin_A{k,i});
        tmp3(i) = mean(signal_bin_mid{k,i});
        tmp4(i) = mean(signal_bin_P{k,i});
        tmp5(i) = std(signal_bin{k,i})/sqrt(length(signal_bin{k,i}));
        tmp6(i) = std(signal_bin_A{k,i}/sqrt(length(signal_bin_A{k,i})));
        tmp7(i) = std(signal_bin_mid{k,i}/sqrt(length(signal_bin_mid{k,i})));
        tmp8(i) = std(signal_bin_P{k,i}/sqrt(length(signal_bin_P{k,i})));
        
    end
    aveSignal_bin(k,:) = tmp1(:);
    aveSignal_bin_A(k,:) = tmp2(:);
    aveSignal_bin_mid(k,:) = tmp3(:);
    aveSignal_bin_P(k,:) = tmp4(:);
    stdSignal_bin(k,:) = tmp5(:);
    stdSignal_bin_A(k,:) = tmp6(:);
    stdSignal_bin_mid(k,:) = tmp7(:);
    stdSignal_bin_P(k,:) = tmp8(:);
end

aveArea_bin(isnan(aveArea_bin)) = 0;
stdArea_bin(isnan(stdArea_bin)) = 0;
aveArea_bin_A(isnan(aveArea_bin_A)) = 0;
aveArea_bin_mid(isnan(aveArea_bin_mid)) = 0;
aveArea_bin_P(isnan(aveArea_bin_P)) = 0;
aveSignal_bin(isnan(aveSignal_bin)) = 0;
aveSignal_bin_A(isnan(aveSignal_bin_A)) = 0;
aveSignal_bin_mid(isnan(aveSignal_bin_mid)) = 0;
aveSignal_bin_P(isnan(aveSignal_bin_P)) = 0;
stdSignal_bin(isnan(stdSignal_bin)) = 0;
stdSignal_bin_A(isnan(stdSignal_bin_A)) = 0;
stdSignal_bin_mid(isnan(stdSignal_bin_mid)) = 0;
stdSignal_bin_P(isnan(stdSignal_bin_P)) = 0;


figure, plot([1:bins]/bins,aveArea_bin(end,:))
xlabel('Location')
ylabel('Nuclei Output')

figure, errorbar([1:bins]/bins,aveSignal_bin(end,:),stdSignal_bin(end,:))
title('Signal Output Across Embryo')
xlabel('Location')
ylabel('Nuclei Signal')

figure, plot([1:bins]/bins,aveSignal_bin_A(2,:))
title('Signal Output Across Embryo in Anterior Region')
xlabel('Location')
ylabel('Nuclei Signal')

figure, plot([1:bins]/bins,aveSignal_bin_mid(2,:))
title('Signal Output Across Embryo in Middle Region')
xlabel('Location')
ylabel('Nuclei Signal')

figure, errorbar([1:bins]/bins,aveSignal_bin_P(52,:),stdSignal_bin_P(52,:))
title('Signal Output Across Embryo in Posterior Region')
xlabel('Location')
ylabel('Nuclei Signal')

imgname = [imagepath 't48Rdldeln3.tif'];
long = length(imfinfo(imgname));
img = imread(imgname,145*2);
figure, imshow(img*10); hold on;
% plot(lineage_cx(74,(activeM([222 342 412]))),lineage_cy(74,(activeM([222 342 412]))),'wo');
% plot(lineage_cx(time,(442)),lineage_cy(time,(442)),'bo');

plot([0:512],liney3,'w',[0:512],liney4,'w')
plot([0:512],liney1,'w',[0:512],liney2,'w')%,[0:512],liney3,'w',[0:512],liney4,'r'); 

% close all