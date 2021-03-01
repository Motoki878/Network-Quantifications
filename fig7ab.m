
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program products figure 7(a)(b).
% This figure shows E/I cell category and FVS (Feedback Vertex Set), especially
% (a) shows the difference of ratios of FVS for excitatory or inhibitory neurons.
% (b) shoes the difference of ratios of FVS for excitatory or inhibitory neurons within individual layers.
%
% Author:          Masanori Shimono     2018/2019
%        Modified  Motoki Kajiwara      2019/2020
%        Cleaned  Ritsuki Nomura        2020/2021
% contact address: shimonomlab@gmail.com
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

load ./orig_data/2019Dec27Result/data
load ./orig_data/2019Dec27Result/data_test

kkm = 1;

for data_num = [1,2,3,4,5,6,7]
    
    data_name1  = ['data',num2str(data_num)];
    data_name3 = data_name1;
    
    str = ['load  ./orig_data/', data_name1 ,'/Aft_Limited/folder_name4;'];     eval(str);
    str = ['./orig_data/newfvslist_200827_ver2/nodetype',num2str(data_name3),'.txt'];
    FVS_index_all{kkm} = load(str);
    
    str = ['load  ./orig_data/', data_name1 ,'/Aft_Limited/layer_lim;'];     eval(str);
    layer_all{kkm} = layer;
    str = ['load  ./orig_data/', data_name1 ,'/Aft_Limited/cell_categ;'];     eval(str);
    cell_categ_all{kkm}  = cell_categ.exc;
    
    %%
    str = ['./orig_data/ave_graph/KC_',data_name3,'.txt'];     KC_all{kkm}  = load(str)';
    
    str = ['load  ./orig_data/ave_graph/motif_hist;'];          eval(str);
    str = ['load  ./orig_data/ave_graph/motif_hist_all;'];      eval(str);
    str = ['load  ./orig_data/ave_graph/motif_hist_rand;'];     eval(str);
    str = ['load  ./orig_data/ave_graph/motif_hist_rand_all;']; eval(str);
    
    %% %% we need to eliminate here
    
    FVSmax = length(FVS_index_all{kkm});
    
    layer.L1_lim  = layer.L1_lim( layer.L1_lim  <= FVSmax);
    layer.L23_lim = layer.L23_lim( layer.L23_lim <= FVSmax);
    layer.L4_lim  = layer.L4_lim( layer.L4_lim  <= FVSmax);
    layer.L5_lim  = layer.L5_lim( layer.L5_lim  <= FVSmax);
    layer.L6_lim  = layer.L6_lim( layer.L6_lim  <= FVSmax);
    layer.cortex_lim  = layer.cortex_lim( 1:FVSmax);
    KC_all{kkm} = KC_all{kkm}( 1 : FVSmax );
    cell_categ_all{kkm} = cell_categ_all{kkm}( 1:FVSmax );
    
    %%
    layvec.L1  = zeros(length(layer.cortex_lim),1);    layvec.L1( layer.L1_lim ) = 1;
    layvec.L23 = zeros(length(layer.cortex_lim),1);    layvec.L23( layer.L23_lim) = 1;
    layvec.L4  = zeros(length(layer.cortex_lim),1);    layvec.L4(  layer.L4_lim ) = 1;
    layvec.L5  = zeros(length(layer.cortex_lim),1);    layvec.L5(  layer.L5_lim ) = 1;
    layvec.L6  = zeros(length(layer.cortex_lim),1);    layvec.L6(  layer.L6_lim ) = 1;
    
    %%
    P_FC_K(kkm) = sum(sum(FVS_index_all{kkm}==3 & KC_all{kkm}>=1)./length(FVS_index_all{kkm}));
    P_FC_I(kkm) = sum(sum(FVS_index_all{kkm}==3 & cell_categ_all{kkm} == 1)./length(FVS_index_all{kkm}));
    P_FCI_K(kkm) = sum(sum(FVS_index_all{kkm}>=2 & KC_all{kkm}>=1)./length(FVS_index_all{kkm}));
    P_FCI_I(kkm) = sum(sum(FVS_index_all{kkm}>=2 & cell_categ_all{kkm} == 1)./length(FVS_index_all{kkm}));
    P_K_I(kkm) = sum(sum(KC_all{kkm}>=1 & cell_categ_all{kkm} == 1)./length(FVS_index_all{kkm}));
    P_FC(kkm)  = sum(sum(FVS_index_all{kkm}==3)./length(FVS_index_all{kkm}));
    P_FCI(kkm) = sum(sum(FVS_index_all{kkm}>=2)./length(FVS_index_all{kkm}));
    P_K(kkm)   = sum(sum(KC_all{kkm}>=1)./length(FVS_index_all{kkm}));
    P_I(kkm)   = sum(sum(cell_categ_all{kkm} == 1)./length(FVS_index_all{kkm}));
    
    FVS.ei_L(1,1,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==0 & (layvec.L1==1|layvec.L23==1))==3);
    FVS.ei_L(1,2,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==1 & (layvec.L1==1|layvec.L23==1))==3);
    FVS.ei_L(2,1,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==0 & layvec.L4==1)==3);
    FVS.ei_L(2,2,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==1 & layvec.L4==1)==3);
    FVS.ei_L(3,1,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==0 & layvec.L5==1)==3);
    FVS.ei_L(3,2,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==1 & layvec.L5==1)==3);
    FVS.ei_L(4,1,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==0 & layvec.L6==1)==3);
    FVS.ei_L(4,2,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==1 & layvec.L6==1)==3);
    FVS_I.ei_L(1,1,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==0 & (layvec.L1==1|layvec.L23==1))==2);
    FVS_I.ei_L(1,2,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==1 & (layvec.L1==1|layvec.L23==1))==2);
    FVS_I.ei_L(2,1,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==0 & layvec.L4==1)==2);
    FVS_I.ei_L(2,2,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==1 & layvec.L4==1)==2);
    FVS_I.ei_L(3,1,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==0 & layvec.L5==1)==2);
    FVS_I.ei_L(3,2,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==1 & layvec.L5==1)==2);
    FVS_I.ei_L(4,1,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==0 & layvec.L6==1)==2);
    FVS_I.ei_L(4,2,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==1 & layvec.L6==1)==2);
    FVS_CI.ei_L(1,1,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==0 & (layvec.L1==1|layvec.L23==1))>=2);
    FVS_CI.ei_L(1,2,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==1 & (layvec.L1==1|layvec.L23==1))>=2);
    FVS_CI.ei_L(2,1,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==0 & layvec.L4==1)>=2);
    FVS_CI.ei_L(2,2,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==1 & layvec.L4==1)>=2);
    FVS_CI.ei_L(3,1,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==0 & layvec.L5==1)>=2);
    FVS_CI.ei_L(3,2,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==1 & layvec.L5==1)>=2);
    FVS_CI.ei_L(4,1,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==0 & layvec.L6==1)>=2);
    FVS_CI.ei_L(4,2,kkm)  = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==1 & layvec.L6==1)>=2);
    FVS.inh(kkm) = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==0)==3);
    FVS.exc(kkm) = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==1)==3);
    FVS_I.inh(kkm) = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==0)==2);
    FVS_I.exc(kkm) = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==1)==2);
    FVS_CI.inh(kkm) = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==0)>=2);
    FVS_CI.exc(kkm) = mean(FVS_index_all{kkm}(cell_categ_all{kkm}==1)>=2);
    
    KC.inh(kkm) = mean(KC_all{kkm}(cell_categ_all{kkm}==0));
    KC.exc(kkm) = mean(KC_all{kkm}(cell_categ_all{kkm}==1));
    
    P_FC_K_rate(kkm) = P_FC_K(kkm)/(P_FC(kkm)*P_K(kkm));
    P_FC_I_rate(kkm) = P_FC_I(kkm)/(P_FC(kkm)*P_I(kkm));
    P_K_I_rate(kkm) =  P_K_I(kkm)/(P_I(kkm)*P_K(kkm));
    P_FCI_K_rate(kkm) = P_FCI_K(kkm)/(P_FCI(kkm)*P_K(kkm));
    P_FCI_I_rate(kkm) = P_FCI_I(kkm)/(P_FCI(kkm)*P_I(kkm));
    P_K_I_rate(kkm) =  P_K_I(kkm)/(P_I(kkm)*P_K(kkm));
    
    kkm = kkm+1;
end

kkm = kkm - 1;

%%
FVS.ei_L(isnan(FVS.ei_L)) = 0;
FVS_I.ei_L(isnan(FVS_I.ei_L)) = 0;
FVS_CI.ei_L(isnan(FVS_CI.ei_L)) = 0;

%%
figure(111)

subplot(2,1,1);
bar([1],[mean(FVS_CI.exc)],0.8, 'r'); hold on;
bar([2],[mean(FVS_CI.inh)],0.8, 'b'); hold on;
errorbar([1,2],[mean(FVS_CI.exc),mean(FVS_CI.inh)],[std(FVS_CI.exc),std(FVS_CI.inh)]/sqrt(kkm),'k.','LineWidth',2); hold on;

for kkm = 1:length(FVS_CI.inh)
    plot([1,2],[FVS_CI.exc(kkm),FVS_CI.inh(kkm)],'go-','MarkerSize',7,'LineWidth',1);
end

set(gca,'XLim',[0.2 2.8],'fontsize',15,'fontname','Arial');
xticks([1,2]);  xticklabels({'Exc.','Inh.'});

for ii = 1:length(FVS_CI.exc)
    plot([1,2], [FVS_CI.exc(ii), FVS_CI.inh(ii)] ,'g-'); hold on;
    plot([1,2], [FVS_CI.exc(ii), FVS_CI.inh(ii)] ,'go'); hold on;
end

[p3, h3] = signrank( FVS_CI.exc(:), FVS_CI.inh(:),'tail','left');
xlabel(['data index'],'fontname','Arial','fontsize',18);
title(['FVS (all) (p=',num2str(p3),')'],'fontsize',18,'fontname','Arial');

%%
zscoree1 = (mean(FVS_CI.exc) - mean(FVS_CI.inh))/(var(FVS_CI.exc)/sqrt(7));
zscorei1 = (mean(FVS_CI.exc) - mean(FVS_CI.inh))/(var(FVS_CI.inh)/sqrt(7));
disp( ['  FVS (exc): ',num2str(mean(FVS_CI.exc)) ,' +/- ', num2str(var(FVS_CI.exc)) ] );
disp( ['  FVS (Inh): ',num2str(mean(FVS_CI.inh)) ,' +/- ', num2str(var(FVS_CI.inh)) ] );
disp( ['  z-scores: ',num2str( zscoree1),'/',num2str( zscorei1)] );
disp( ['  p-values: ',num2str(p3)] );

%%
subplot(2,1,2);
bar([1:4]-0.15,[squeeze(mean(FVS_CI.ei_L(:,2,:),3))],0.27, 'r'); hold on;
bar([1:4]+0.15,[squeeze(mean(FVS_CI.ei_L(:,1,:),3))],0.27 ,'b'); hold on;
set(gca,'XLim',[0.3 4.7],'fontsize',15,'fontname','Arial');
xticks([1:4]);  xticklabels({'L1-3','L4','L5','L6'});
disp(['Layer Group1: layer1-3, Group2: layer4, Group3: layer5, Group4: layer6']);

for nnl = 1:4
    FVS_exc1 = squeeze(FVS_CI.ei_L(nnl,2,:))';
    FVS_inh1 = squeeze(FVS_CI.ei_L(nnl,1,:))';
    FVS_ei = [FVS_exc1; FVS_inh1];
    plot(nnl-0.375+[1,2]/4, FVS_ei ,'g-'); hold on;
    plot(nnl-0.375+[1,2]/4, FVS_ei ,'go'); hold on;
    
    %%
    disp(['Layer Group',num2str(nnl)]);
    
    %%
    zscoree1 = (mean(FVS_exc1) - mean(FVS_inh1))/(var(FVS_exc1)/sqrt(7));
    zscorei1 = (mean(FVS_exc1) - mean(FVS_inh1))/(var(FVS_inh1)/sqrt(7));
    
    [p3, h3] = signrank( FVS_exc1, FVS_inh1);
    
    disp( ['   FVS (exc): ',num2str(mean(FVS_exc1)) ,' +/- ', num2str(var(FVS_exc1)) ] );
    disp( ['   FVS (Inh): ',num2str(mean(FVS_inh1)) ,' +/- ', num2str(var(FVS_inh1)) ] );
    disp( ['   z-scores: ',num2str( zscoree1),'/',num2str( zscorei1)] );
    disp( ['   p-values: ',num2str(p3)] );
    
end

for kkl = 1:4
    errorbar(kkl-0.15, squeeze(mean(FVS_CI.ei_L(kkl,2,:),3)), std(squeeze(FVS_CI.ei_L(kkl,2,:)))/sqrt(kkm) ,'b.');
    [p03(kkl), h03(kkl)] = signrank( squeeze(FVS_CI.ei_L(kkl,2,:)), squeeze(FVS_CI.ei_L(kkl,1,:)) ,'tail','left');
    errorbar(kkl+0.15, squeeze(mean(FVS_CI.ei_L(kkl,1,:),3)), std(squeeze(FVS_CI.ei_L(kkl,1,:)))/sqrt(kkm) ,'b.');
end

for kkl = 1:4
    [p03(kkl), h03(kkl)] = signrank( squeeze(FVS_CI.ei_L(kkl,2,:)), squeeze(FVS_CI.ei_L(kkl,1,:)) ,'tail','left');
    
    title(['FVS layers (all) (p=',num2str(p03(kkl)),')'],'fontsize',12,'fontname','Arial');
end

%%
for kk = 1:length(motif_hist_all)
    motif_hist_all{kk} ;
end

clear T1 T2 T3 C1 C2 C3

%%
filename = 'FVS_layers.xlsx';

for subj1 = 1:7  
    T1 = table(squeeze(FVS_CI.ei_L(:,2,subj1)));
    T2 = table(squeeze(FVS_CI.ei_L(:,1,subj1)));
    T1_all = table(squeeze(FVS_CI.exc(subj1)));
    T2_all = table(squeeze(FVS_CI.inh(subj1)));
    
    C1 = table({'FVS (exc)';'layer1-3';'layer4';'layer5';'layer6';'all'});
    C2 = table({'FVS (inh)';'layer1-3';'layer4';'layer5';'layer6';'all'});
    
    writetable(C1,filename,'Sheet',subj1,'Range','A1');
    writetable(C2,filename,'Sheet',subj1,'Range','D1');
    writetable(T1_all,filename,'Sheet',subj1,'Range','B6');
    writetable(T2_all,filename,'Sheet',subj1,'Range','E6');
    writetable(T1,filename,'Sheet',subj1,'Range','B2');
    writetable(T2,filename,'Sheet',subj1,'Range','E2');   
end

