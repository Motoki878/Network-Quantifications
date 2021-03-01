
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program products figure 5(d).
% Figure5 shows basic properties of neuronal networks, and especially
% (d) shows the total number of identified excitatory and inhibitory neurons.
%
% Author:          Masanori Shimono     2018/2019
%        Modified  Motoki Kajiwara      2019/2020
%        Cleaned  Ritsuki Nomura        2020/2021
% contact address: nori417@gmail.com
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off
close all;
clear all;

addpath(genpath('./'))

kk = 1;

for data_num = [1,2,3,4,5,6,7]
    
    data_name1  = ['data',num2str(data_num)];
    data_name3  = data_name1;
    
    str = ['load  ./orig_data/',data_name1,'/Aft_Limited/folder_name4;'];     eval(str);
    str = ['load  ./orig_data/',data_name1,'/Aft_Limited/cell_categ;'];     eval(str);
    str = ['load  ./orig_data/',data_name1,'/Aft_Limited/pdf;'];     eval(str);
    str = ['load  ./orig_data/',data_name1,'/Aft_Limited/layer_lim;'];     eval(str);
    
    cortex = [layer.L1; layer.L23; layer.L4; layer.L5; layer.L6];
    cortex_L15 = [layer.L1; layer.L23; layer.L4; layer.L5];
    
    cortex_L15_lim = round(length(layer.cortex_lim)*length(cortex_L15)/length(cortex));
    
    length(cortex);
    
    %% dividing between excitatory and inhibitory neurons
    exc_inh = 1-cell_categ.inh;
    
    %%
    cell_num(kk,1) = sum(exc_inh==1);
    cell_num(kk,2) = sum(exc_inh==0);
    
    pdf_size = length(pdf);
    
    kk = kk + 1;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(183)
subplot(2,2,1)
bar( mean(cell_num(:,[1,2]),1) ); hold on;
errorbar([1,2], mean(cell_num(:,[1,2]),1), std(cell_num(:,[1,2]))/sqrt(kk-1),'b.');

xticks([1,2]);  xticklabels({'Exc.','Inh.'});
ylabel(['Cell numbers'],'fontsize', 10,'fontname','Arial');
xlabel(['Cell types'],'fontsize', 10,'fontname','Arial');

ratio = cell_num(:,[2])./cell_num(:,[1]);
pd = fitdist(ratio,'Normal');
ci = paramci(pd);

title(['CI: [', num2str(ci(1,1)) ,',', num2str(ci(2,1))  ,']'],'fontsize', 10,'fontname','Arial');

subplot(2,2,2)
pie(mean(cell_num(:,[2,1]),1) ); hold on;
ylabel(['Cell numbers'],'fontsize', 10,'fontname','Arial');
xlabel(['subjects #'],'fontsize', 10,'fontname','Arial');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(183)
subplot(2,2,3)
errorbar([1], mean(cell_num(:,[1]),1) , std(cell_num(:,[1]))/sqrt(kk-1) ,'r+','LineWidth',2); hold on;
errorbar([2], mean(cell_num(:,[2]),1) , std(cell_num(:,[2]))/sqrt(kk-1) ,'b+','LineWidth',2); hold on;

for kkm = 1:7
    plot([0.9]+kkm/(7*5) , cell_num(kkm,[1]) ,'r.'); hold on;
    plot([1.9]+kkm/(7*5) , cell_num(kkm,[2]) ,'b.'); hold on;
end

xticks([1,2]);  xticklabels({'Exc.','Inh.'});
ylabel(['Cell numbers'],'fontsize', 18,'fontname','Arial');
xlabel(['Cell types'],'fontsize', 18,'fontname','Arial');

set(gca,'XLim',[0.5 2.5],'fontsize',12,'fontname','Arial');

%%
T1 = table(cell_num(:,1));
T2 = table(cell_num(:,2));

C1 = table({'excitatory'})
C2 = table({'inhibitory'})

filename = 'cell_num.xlsx';
writetable(T1,filename,'Sheet',1,'Range','B2')
writetable(T2,filename,'Sheet',2,'Range','B2')

writetable(C1,filename,'Sheet',1,'Range','A2')
writetable(C2,filename,'Sheet',2,'Range','A2')
