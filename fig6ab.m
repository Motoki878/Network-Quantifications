
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program products figure 6(a)(b). This figure shows E/I cell categories and k-core centralities, especially
% (a) shows the difference of averaged k-core values for excitatory and inhibitory neurons, and
% (b) shows the difference of averaged k-core values for excitatory and inhibitory neurons within individual layers.
% 
% Author:          Masanori Shimono     2018/2019
%        Modified  Motoki Kajiwara      2019/2020
%        Cleaned  Ritsuki Nomura        2020/2021
% contact address: shimonomlab@gmail.com
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off;
close all;   clear all;

kk = 1;

for data_num = [1,2,3,4,5,6,7]
    
    clear layer
    
    data_name1  = ['data',num2str(data_num)];
    data_name3  = data_name1;
    
    str = ['load  ./orig_data/',data_name1,'/Aft_Limited/folder_name4;'];     eval(str);
    
    %% Loading data
    
    str = ['load  ./orig_data/ave_graph/network_measure_str;'];     eval(str);
    str = ['load  ./orig_data/',data_name1,'/Aft_Limited/cell_categ;'];     eval(str);
    str = ['load  ./orig_data/',data_name1,'/Aft_Limited/pdf;'];     eval(str);
    str = ['load  ./orig_data/',data_name1,'/Aft_Limited/synapse_label;'];     eval(str);
    str = ['load  ./orig_data/',data_name1,'/Aft_Limited/FR;'];     eval(str);
    str = ['load  ./orig_data/',data_name1,'/Aft_Limited/conn_prob;'];     eval(str);
    str = ['load  ./orig_data/',data_name1,'/Aft_Limited/layer_lim;'];     eval(str);
    
    etd = ['load ./orig_data/',data_name1,'/points_NeuN;']; eval(etd);

    %%
    layer_L1 = layer.L1; layer_L23 = layer.L23;  layer_L4 = layer.L4;  layer_L5 = layer.L5; layer_L6 = layer.L6;
    
    %% layers -- select neurons within distance < 2000
    kkm1   = 1;    kkm23   = 1;    kkm4   = 1;    kkm5   = 1;    kkm6   = 1;
    kkm1_e = 1;    kkm23_e = 1;    kkm4_e = 1;    kkm5_e = 1;    kkm6_e = 1;
    kkm1_i = 1;    kkm23_i = 1;    kkm4_i = 1;    kkm5_i = 1;    kkm6_i = 1;
    
    layer.L1_lim = [];    layer.L1e_lim = [];   layer.L1i_lim = [];
    layer.L23_lim = [];   layer.L23e_lim = [];  layer.L23i_lim = [];
    layer.L4_lim = [];    layer.L4e_lim = [];   layer.L4i_lim = [];
    layer.L5_lim = [];    layer.L5e_lim = [];   layer.L5i_lim = [];
    layer.L6_lim = [];    layer.L6e_lim = [];   layer.L6i_lim = [];
    
    for kk00 = 1:length(layer.cortex_lim)
        kk0 = layer.cortex_lim(kk00);
        %%
        if length(layer.L1) > 0;  for kk10 = 1:length(layer.L1);     kk1 = layer.L1(kk10);
                if kk0 == kk1;                        layer.L1_lim(kkm1)    = kk00;  kkm1 = kkm1 + 1;
                    if     cell_categ.exc(kk00) == 1;  layer.L1e_lim(kkm1_e) = kk00;  kkm1_e = kkm1_e + 1;
                    elseif cell_categ.exc(kk00) == 0;  layer.L1i_lim(kkm1_i) = kk00;  kkm1_i = kkm1_i + 1;  end
                end;  end;  end;
        %%
        if length(layer.L23) > 0;     for kk10 = 1:length(layer.L23);    kk1 = layer.L23(kk10);
                if kk0 == kk1;                       layer.L23_lim(kkm23)    = kk00;  kkm23 = kkm23 + 1;
                    if     cell_categ.exc(kk00) == 1;  layer.L23e_lim(kkm23_e) = kk00;  kkm23_e = kkm23_e + 1;
                    elseif cell_categ.exc(kk00) == 0;  layer.L23i_lim(kkm23_i) = kk00;  kkm23_i = kkm23_i + 1;  end
                end;  end;  end;
        %%
        if length(layer.L4) > 0;  for kk10 = 1:length(layer.L4)-1;     kk1 = layer.L4(kk10);
                if kk0 == kk1;                       layer.L4_lim(kkm4)    = kk00;  kkm4 = kkm4 + 1;
                    if     cell_categ.exc(kk00) == 1;  layer.L4e_lim(kkm4_e) = kk00;  kkm4_e = kkm4_e + 1;
                    elseif cell_categ.exc(kk00) == 0;  layer.L4i_lim(kkm4_i) = kk00;  kkm4_i = kkm4_i + 1;  end
                end;  end;  end;
        %%
        if length(layer.L5) > 0;  for kk10 = 1:length(layer.L5);  kk1 = layer.L5(kk10);
                if kk0 == kk1;                       layer.L5_lim(kkm5)    = kk00;  kkm5 = kkm5 + 1;
                    if     cell_categ.exc(kk00) == 1;  layer.L5e_lim(kkm5_e) = kk00;  kkm5_e = kkm5_e + 1;
                    elseif cell_categ.exc(kk00) == 0;  layer.L5i_lim(kkm5_i) = kk00;  kkm5_i = kkm5_i + 1;  end
                end;   end;  end;
        %%
        if length(layer.L6) > 0;  for kk10 = 1:length(layer.L6);    kk1 = layer.L6(kk10);
                if kk0 == kk1;                       layer.L6_lim(kkm6)    = kk00;  kkm6 = kkm6 + 1;
                    if     cell_categ.exc(kk00) == 1;  layer.L6e_lim(kkm6_e) = kk00;  kkm6_e = kkm6_e + 1;
                    elseif cell_categ.exc(kk00) == 0;  layer.L6i_lim(kkm6_i) = kk00;  kkm6_i = kkm6_i + 1;  end
                end;  end;  end;
    end
    
    %% recover the stored values (19/8/22)
    layer.L1 = layer_L1; layer.L23 = layer_L23;   layer.L4 = layer_L4;    layer.L5 = layer_L5;    layer.L6 = layer_L6;
    
    %% Layers summary
    cortex     = [layer.L1_lim, layer.L23_lim, layer.L4_lim, layer.L5_lim, layer.L6_lim]';
    cortex_L15 = [layer.L1_lim, layer.L23_lim, layer.L4_lim, layer.L5_lim]';
    cortex_e     = [layer.L1e_lim, layer.L23e_lim, layer.L4e_lim, layer.L5e_lim, layer.L6e_lim]';
    cortex_i     = [layer.L1i_lim, layer.L23i_lim, layer.L4i_lim, layer.L5i_lim, layer.L6i_lim]';
    
    %% representative networks measures on diff. layers
    
    for net_loop = 1:23
        %% all             {measure} (subject, layer, mean/std)
        network_measure_str2.allsub_Lcateg{net_loop}(kk,1,1) = mean(network_measure_str.allsub{net_loop}{kk}([layer.L1_lim,layer.L23_lim]));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,1,2) =  std(network_measure_str.allsub{net_loop}{kk}([layer.L1_lim,layer.L23_lim]));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,2,1) = mean(network_measure_str.allsub{net_loop}{kk}(layer.L4_lim));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,2,2) =  std(network_measure_str.allsub{net_loop}{kk}(layer.L4_lim));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,3,1) = mean(network_measure_str.allsub{net_loop}{kk}(layer.L5_lim));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,3,2) =  std(network_measure_str.allsub{net_loop}{kk}(layer.L5_lim));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,4,1) = mean(network_measure_str.allsub{net_loop}{kk}(layer.L6_lim));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,4,2) =  std(network_measure_str.allsub{net_loop}{kk}(layer.L6_lim));
        
        %% excitatory           {measure} (subject, layer, mean/std)
        network_measure_str2.allsub_Lcateg{net_loop}(kk,1,3) = mean(network_measure_str.allsub{net_loop}{kk}([layer.L1e_lim,layer.L23e_lim]));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,1,4) =  std(network_measure_str.allsub{net_loop}{kk}([layer.L1e_lim,layer.L23e_lim]));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,2,3) = mean(network_measure_str.allsub{net_loop}{kk}(layer.L4e_lim));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,2,4) =  std(network_measure_str.allsub{net_loop}{kk}(layer.L4e_lim));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,3,3) = mean(network_measure_str.allsub{net_loop}{kk}(layer.L5e_lim));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,4,4) =  std(network_measure_str.allsub{net_loop}{kk}(layer.L5e_lim));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,4,3) = mean(network_measure_str.allsub{net_loop}{kk}(layer.L6e_lim));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,4,4) =  std(network_measure_str.allsub{net_loop}{kk}(layer.L6e_lim));
        
        %% inhhibitory           {measure} (subject, layer, mean/std)
        network_measure_str2.allsub_Lcateg{net_loop}(kk,1,5) = mean(network_measure_str.allsub{net_loop}{kk}([layer.L1i_lim,layer.L23i_lim]));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,1,6) =  std(network_measure_str.allsub{net_loop}{kk}([layer.L1i_lim,layer.L23i_lim]));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,2,5) = mean(network_measure_str.allsub{net_loop}{kk}(layer.L4i_lim));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,2,6) =  std(network_measure_str.allsub{net_loop}{kk}(layer.L4i_lim));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,3,5) = mean(network_measure_str.allsub{net_loop}{kk}(layer.L5i_lim));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,3,6) =  std(network_measure_str.allsub{net_loop}{kk}(layer.L5i_lim));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,4,5) = mean(network_measure_str.allsub{net_loop}{kk}(layer.L6i_lim));
        network_measure_str2.allsub_Lcateg{net_loop}(kk,4,6) =  std(network_measure_str.allsub{net_loop}{kk}(layer.L6i_lim));
        network_measure_str2.allsub_Lcateg{net_loop}(isnan(network_measure_str2.allsub_Lcateg{net_loop})) = 0;
        network_measure_str2.allsub{net_loop}(kk,1) = mean(network_measure_str.allsub{net_loop}{kk}(cortex_e));
        network_measure_str2.allsub{net_loop}(kk,2) = mean(network_measure_str.allsub{net_loop}{kk}(cortex_i));
        
    end
    
    %%
    
    neunum_layer(kk,:,:) = [length(layer.L1i_lim)+length(layer.L23i_lim) , length(layer.L1e_lim)+length(layer.L23e_lim);...
    length(layer.L4i_lim) , length(layer.L4e_lim);...
    length(layer.L5i_lim) , length(layer.L5e_lim);...
    length(layer.L6i_lim) , length(layer.L6e_lim)];

neunum_layer(kk,:,:) = [length(layer.L1i_lim)+length(layer.L23i_lim) , length(layer.L1e_lim)+length(layer.L23e_lim);...
length(layer.L4i_lim) , length(layer.L4e_lim);...
length(layer.L5i_lim) , length(layer.L5e_lim);...
length(layer.L6i_lim) , length(layer.L6e_lim)];

kk = kk+1;

end


%% eliminating data in which layer 1 was not included
neunum_layer(isnan(neunum_layer))=0;
non_rec_layer = ((neunum_layer(:,1,1)+neunum_layer(:,1,2))>0);

%% calulating densities

BC_bin_IE_ratio = squeeze(network_measure_str2.allsub_Lcateg{6}(:,:,5)./...
    (network_measure_str2.allsub_Lcateg{6}(:,:,3)));
BC_exc_mean = squeeze(mean(network_measure_str2.allsub_Lcateg{6}(:,:,3),1));
BC_inh_mean = squeeze(mean(network_measure_str2.allsub_Lcateg{6}(:,:,5),1));
BC_exc_std = squeeze(std(network_measure_str2.allsub_Lcateg{6}(:,:,3)/sqrt(kk),1));
BC_inh_std = squeeze(std(network_measure_str2.allsub_Lcateg{6}(:,:,5)/sqrt(kk),1));
BC_bin_IE_ratio = squeeze((BC_inh_mean+0.0001)./(BC_exc_mean+BC_inh_mean+0.0001));
BC_exc_all =  squeeze(network_measure_str2.allsub_Lcateg{6}(:,:,3));
BC_inh_all =  squeeze(network_measure_str2.allsub_Lcateg{6}(:,:,5));

for kkm = 1:size(BC_exc_all,2)
    [p, res] = signrank( BC_exc_all(:,kkm) , BC_inh_all(:,kkm) ,'tail','left');
    p_bc_all(kkm) = p;
end

p_bc_all

%%
neunum_layer(isnan(neunum_layer))=0;

%%

KC_bin_IE_ratio = squeeze(network_measure_str2.allsub_Lcateg{4}(:,:,5)./...
    network_measure_str2.allsub_Lcateg{4}(:,:,3));


disp(['Layer Group1: layer1-3, Group2: layer4, Group3: layer5, Group4: layer6']);

for kkn = [1:size(network_measure_str2.allsub_Lcateg{4},2)]
    clear aae aai
    aae = squeeze(network_measure_str2.allsub_Lcateg{4}(:,kkn,3));  % aae = aae(aae~=0);
    aai = squeeze(network_measure_str2.allsub_Lcateg{4}(:,kkn,5));  % aai = aai(aai~=0);
    
    %%
    disp(['Layer Group',num2str(kkn)]);
    
    KC_exc_mean(kkn) = squeeze(mean(aae,1));
    KC_inh_mean(kkn) = squeeze(mean(aai,1));
    KC_exc_std(kkn) = squeeze(std(aae,1))/sqrt(kkm);
    KC_inh_std(kkn) = squeeze(std(aai,1))/sqrt(kkm);
    KC_exc_var(kkn) = squeeze(var(aae,1));
    KC_inh_var(kkn) = squeeze(var(aai,1));
    
    zscorei(kkn) = ( KC_exc_mean(kkn) - KC_inh_mean(kkn))/(KC_inh_var(kkn)/sqrt(7));
    zscoree(kkn) = ( KC_exc_mean(kkn) - KC_inh_mean(kkn))/(KC_exc_var(kkn)/sqrt(7));
    
    disp(['  KC(E) = ',num2str( KC_exc_mean(kkn)),' +/- ', num2str( KC_exc_var(kkn))]);
    disp(['  KC(I) = ',num2str( KC_inh_mean(kkn)),' +/- ', num2str( KC_inh_var(kkn))]);
    disp(['  Std of KC(E/I) = ',num2str( KC_exc_std(kkn)),'/', num2str( KC_inh_std(kkn))]);
    disp(['  Z-score(E/I) = ',num2str( zscoree(kkn)),'/',num2str( zscorei(kkn))]);
    
    
    
end

KC_bin_IE_ratio = squeeze((KC_inh_mean+0.0001)./(KC_exc_mean+0.0001));

%%
KC_exc_all =  squeeze(network_measure_str2.allsub_Lcateg{4}(:,:,3));
KC_inh_all =  squeeze(network_measure_str2.allsub_Lcateg{4}(:,:,5));

for kkm = 1:size(KC_exc_all,2)
    
    disp(['Layer Group',num2str(kkm)]);
    [p, res] = signrank( KC_exc_all(:,kkm) , KC_inh_all(:,kkm), 'tail','left');
    p_kc_all(kkm) = p;
    disp(['  p-val  = ',num2str(  p_kc_all(kkm) )]);
    
end

p_kc_all

%%
neunum_layer(isnan(neunum_layer))=0;

%%
figure(189)

%%
subplot(2,1,1)
data0 = mean(network_measure_str.ave_allsub{4},1);
bar( [1], data0(:,1) , 'r' ); hold on;
bar( [2], data0(:,2) , 'b' ); hold on;

[p1, h1] = signrank(network_measure_str.ave_allsub{4}(:,1), network_measure_str.ave_allsub{4}(:,2)); hold on;
title(['K-Core Centrality (p=',num2str(p1),')'],'fontsize',10,'fontname','Arial');
xticks([1,2]);  xticklabels({'Exc.','Inh.'});

for ii = 1:size(network_measure_str.ave_allsub{4},1)
    plot([1,2], network_measure_str.ave_allsub{4}(ii,:),'g-'); hold on;
    plot([1,2], network_measure_str.ave_allsub{4}(ii,:),'go'); hold on;
end
errorbar([1,2], mean(network_measure_str.ave_allsub{4},1), std(network_measure_str.ave_allsub{4})/sqrt(kk-1) ,'k.'); hold on;

ylabel(['K-Core Centrality'],'fontsize',10,'fontname','Arial');
set(gca,'XLim',[0.5 2.5],'fontsize',12,'fontname','Arial');


subplot(2,1,2)
data1 = [KC_exc_mean; KC_inh_mean]';
bar([1:4]-0.15, data1(:,1) ,0.27, 'r' ); hold on;
bar([1:4]+0.15, data1(:,2) ,0.27, 'b' ); hold on;
errorbar([1:4]+0.15, [KC_inh_mean]' , [KC_inh_std]','bo' ); hold on;
errorbar([1:4]-0.15, [KC_exc_mean]' , [KC_exc_std]','bo' ); hold on;
xticks([1:4]);  xticklabels({'L1-3','L4','L5','L6'});
title(num2str(p_kc_all) ,'fontsize',12,'fontname','Arial');
ylabel(['K-Core Centrality'],'fontsize',12,'fontname','Arial');
xlabel(['layer #'],'fontsize',12,'fontname','Arial');

for nnl = 1:4
    KC_exc1 = squeeze(network_measure_str2.allsub_Lcateg{4}(:,nnl,3))';
    KC_inh1 = squeeze(network_measure_str2.allsub_Lcateg{4}(:,nnl,5))';
    KC_ei = [KC_exc1; KC_inh1];
    plot(nnl-0.375+[1,2]/4, KC_ei ,'g-'); hold on;
    plot(nnl-0.375+[1,2]/4, KC_ei ,'go'); hold on;
    KC_ei_lay_all(:,:,nnl) = KC_ei;
end

xlabel(['layer #'],'fontsize',12,'fontname','Arial');

%%

clear T1 T2 T3 C1 C2 C3

filename = 'KC_layers.xlsx';

for subj1 = 1:7
    
    T1 = table(squeeze(KC_ei_lay_all(1,subj1,:)));
    T2 = table(squeeze(KC_ei_lay_all(2,subj1,:)));
    T1_all = table(squeeze( network_measure_str2.allsub{4}(subj1,1) ));
    T2_all = table(squeeze( network_measure_str2.allsub{4}(subj1,2) ));
    
    C1 = table({'KC (exc)';'layer1-3';'layer4';'layer5';'layer6';'all'});
    C2 = table({'KC (inh)';'layer1-3';'layer4';'layer5';'layer6';'all'});
    
    writetable(C1,filename,'Sheet',subj1,'Range','A1');
    writetable(C2,filename,'Sheet',subj1,'Range','D1');
    writetable(T1_all,filename,'Sheet',subj1,'Range','B6');
    writetable(T2_all,filename,'Sheet',subj1,'Range','E6');
    writetable(T1,filename,'Sheet',subj1,'Range','B2');
    writetable(T2,filename,'Sheet',subj1,'Range','E2');
    
end

%%
clear high_BC_index

topo_index0 = 1;
for topo_index = [4,6]
    for kkn = 1:kk-1
        mm_e = 1;
        high_KC_BC_index.exc{kkn} = [];
        for nn = 1:length( network_measure_str.exc_allsub{topo_index}{kkn})
            
            if network_measure_str.exc_allsub{topo_index}{kkn}(nn) == max( network_measure_str.exc_allsub{topo_index}{kkn}  )
                high_KC_BC_index.exc{kkn} = [high_KC_BC_index.exc{kkn},nn];
            end
            
            high_KC_BC_index1.exc{kkn} = 0;
            if network_measure_str.exc_allsub{topo_index}{kkn}(nn) >= 1
                high_KC_BC_index1.exc{kkn} = [high_KC_BC_index.exc{kkn},nn];
            end
            
            high_KC_BC_index2.exc{kkn} = 0;
            if network_measure_str.exc_allsub{topo_index}{kkn}(nn) >= 2
                high_KC_BC_index2.exc{kkn} = [high_KC_BC_index.exc{kkn},nn];
            end
            
            high_KC_BC_index3.exc{kkn} = 0;
            if network_measure_str.exc_allsub{topo_index}{kkn}(nn) >= 3
                high_KC_BC_index3.exc{kkn} = [high_KC_BC_index.exc{kkn},nn];
            end
            
        end
        
        disp(['Max k-core (exc)', num2str(max( network_measure_str.exc_allsub{topo_index}{kkn}))]);
        
        high_KC_BC_index.exc_num{topo_index0}(kkn) = length(high_KC_BC_index.exc{kkn});
        high_KC_BC_index.exc_ratio{topo_index0}(kkn) = length(high_KC_BC_index.exc{kkn})/length( network_measure_str.exc_allsub{topo_index}{kkn});
        
        high_KC_BC_index1.exc_ratio{topo_index0}(kkn) = length(high_KC_BC_index1.exc{kkn})/length( network_measure_str.exc_allsub{topo_index}{kkn});
        high_KC_BC_index2.exc_ratio{topo_index0}(kkn) = length(high_KC_BC_index2.exc{kkn})/length( network_measure_str.exc_allsub{topo_index}{kkn});
        high_KC_BC_index3.exc_ratio{topo_index0}(kkn) = length(high_KC_BC_index3.exc{kkn})/length( network_measure_str.exc_allsub{topo_index}{kkn});
        
        disp(['Max k-core (inh)', num2str(max( network_measure_str.inh_allsub{topo_index}{kkn}))]);
        mm_i = 1;
        high_KC_BC_index.inh{kkn} = [];
        for nn = 1:length( network_measure_str.inh_allsub{topo_index}{kkn})
            
            if network_measure_str.inh_allsub{topo_index}{kkn}(nn) == max( network_measure_str.inh_allsub{topo_index}{kkn}  )
                high_KC_BC_index.inh{kkn} = [high_KC_BC_index.inh{kkn},nn];
            end
            
            high_KC_BC_index1.inh{kkn} = 0;
            if network_measure_str.inh_allsub{topo_index}{kkn}(nn) >= 1
                high_KC_BC_index1.inh{kkn} = [high_KC_BC_index.inh{kkn},nn];
            end
            
            high_KC_BC_index2.inh{kkn} = 0;
            if network_measure_str.inh_allsub{topo_index}{kkn}(nn) >= 2
                high_KC_BC_index2.inh{kkn} = [high_KC_BC_index.inh{kkn},nn];
            end
            
            high_KC_BC_index3.inh{kkn} = 0;
            if network_measure_str.inh_allsub{topo_index}{kkn}(nn) >= 3
                high_KC_BC_index3.inh{kkn} = [high_KC_BC_index.inh{kkn},nn];
            end
        end
        
        high_KC_BC_index.inh_num{topo_index0}(kkn) = length(high_KC_BC_index.inh{kkn});
        high_KC_BC_index.inh_ratio{topo_index0}(kkn) = length(high_KC_BC_index.inh{kkn})/length( network_measure_str.inh_allsub{topo_index}{kkn});
        
        high_KC_BC_index1.inh_ratio{topo_index0}(kkn) = length(high_KC_BC_index1.inh{kkn})/length( network_measure_str.inh_allsub{topo_index}{kkn});
        high_KC_BC_index2.inh_ratio{topo_index0}(kkn) = length(high_KC_BC_index2.inh{kkn})/length( network_measure_str.inh_allsub{topo_index}{kkn});
        high_KC_BC_index3.inh_ratio{topo_index0}(kkn) = length(high_KC_BC_index3.inh{kkn})/length( network_measure_str.inh_allsub{topo_index}{kkn});
        
    end
    topo_index0 = topo_index0 + 1;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zscoree0 = (mean(network_measure_str.ave_allsub{4}(:,1)) - mean(network_measure_str.ave_allsub{4}(:,2)))/(std(network_measure_str.ave_allsub{4}(:,1))/sqrt(7));
zscorei0 = (mean(network_measure_str.ave_allsub{4}(:,1)) - mean(network_measure_str.ave_allsub{4}(:,2)))/(std(network_measure_str.ave_allsub{4}(:,2))/sqrt(7));

disp( ['k-core (Exc.): ',num2str(mean(network_measure_str.ave_allsub{4}(:,1))) ,' +/- ', num2str(std(network_measure_str.ave_allsub{4}(:,1))) ] );
disp( ['k-core (Inh.): ',num2str(mean(network_measure_str.ave_allsub{4}(:,2))) ,' +/- ', num2str(std(network_measure_str.ave_allsub{4}(:,2))) ] );
disp( ['  z-scores: ',num2str( zscoree0),'/',num2str( zscorei0)] );
disp( ['  p-values: ',num2str( p1)]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zscoree = (mean(high_KC_BC_index.exc_ratio{1}) - mean(high_KC_BC_index.inh_ratio{1}))/(std(high_KC_BC_index.exc_ratio{1})/sqrt(7));
zscorei = (mean(high_KC_BC_index.exc_ratio{1}) - mean(high_KC_BC_index.inh_ratio{1}))/(std(high_KC_BC_index.inh_ratio{1})/sqrt(7));

[pkce, hkce] = signrank( high_KC_BC_index.exc_ratio{1}(:), high_KC_BC_index.inh_ratio{1}(:));
[pkci, hkci] = signrank( high_KC_BC_index.inh_ratio{1}(:), high_KC_BC_index.exc_ratio{1}(:));

zscoree1 = (mean(high_KC_BC_index1.exc_ratio{1}) - mean(high_KC_BC_index1.inh_ratio{1}))/(std(high_KC_BC_index1.exc_ratio{1})/sqrt(7));
zscorei1 = (mean(high_KC_BC_index1.exc_ratio{1}) - mean(high_KC_BC_index1.inh_ratio{1}))/(std(high_KC_BC_index1.inh_ratio{1})/sqrt(7));

[pkce1, hkce1] = signrank( high_KC_BC_index1.exc_ratio{1}(:), high_KC_BC_index1.inh_ratio{1}(:));
[pkci1, hkci1] = signrank( high_KC_BC_index1.inh_ratio{1}(:), high_KC_BC_index1.exc_ratio{1}(:));

zscoree2 = (mean(high_KC_BC_index2.exc_ratio{1}) - mean(high_KC_BC_index2.inh_ratio{1}))/(std(high_KC_BC_index2.exc_ratio{1})/sqrt(7));
zscorei2 = (mean(high_KC_BC_index3.exc_ratio{1}) - mean(high_KC_BC_index3.inh_ratio{1}))/(std(high_KC_BC_index2.inh_ratio{1})/sqrt(7));

[pkce2, hkce2] = signrank( high_KC_BC_index2.exc_ratio{1}(:), high_KC_BC_index2.inh_ratio{1}(:));
[pkci2, hkci2] = signrank( high_KC_BC_index2.inh_ratio{1}(:), high_KC_BC_index2.exc_ratio{1}(:));

zscoree3 = (mean(high_KC_BC_index3.exc_ratio{1}) - mean(high_KC_BC_index3.inh_ratio{1}))/(std(high_KC_BC_index3.exc_ratio{1})/sqrt(7));
zscorei3 = (mean(high_KC_BC_index3.exc_ratio{1}) - mean(high_KC_BC_index3.inh_ratio{1}))/(std(high_KC_BC_index3.inh_ratio{1})/sqrt(7));

[pkce3, hkce3] = signrank( high_KC_BC_index3.exc_ratio{1}(:), high_KC_BC_index3.inh_ratio{1}(:));
[pkci3, hkci3] = signrank( high_KC_BC_index3.inh_ratio{1}(:), high_KC_BC_index3.exc_ratio{1}(:));

disp( ['High-KC ratio (Exc.): ',num2str(mean(high_KC_BC_index.exc_ratio{1})) ,' +/- ', num2str(var(high_KC_BC_index.exc_ratio{1})) ] );
disp( ['High-KC ratio (Inh.): ',num2str(mean(high_KC_BC_index.inh_ratio{1})) ,' +/- ', num2str(var(high_KC_BC_index.inh_ratio{1})) ] );
disp( ['z-scores: ',num2str( zscoree),'/',num2str( zscorei)] );
disp( ['p-values: ',num2str( pkce),'/',num2str( pkci)] );

disp( ['High-KC ratio (Exc., k-core>=1): ',num2str(mean(high_KC_BC_index1.exc_ratio{1})) ,' +/- ', num2str(var(high_KC_BC_index1.exc_ratio{1})) ] );
disp( ['High-KC ratio (Inh., k-core>=1): ',num2str(mean(high_KC_BC_index1.inh_ratio{1})) ,' +/- ', num2str(var(high_KC_BC_index1.inh_ratio{1})) ] );
disp( ['z-scores: ',num2str( zscoree1),'/',num2str( zscorei1)] );
disp( ['p-values: ',num2str( pkce1),'/',num2str( pkci1)] );

disp( ['High-KC ratio (Exc., k-core>=2): ',num2str(mean(high_KC_BC_index2.exc_ratio{1})) ,' +/- ', num2str(var(high_KC_BC_index2.exc_ratio{1})) ] );
disp( ['High-KC ratio (Inh., k-core>=2): ',num2str(mean(high_KC_BC_index2.inh_ratio{1})) ,' +/- ', num2str(var(high_KC_BC_index2.inh_ratio{1})) ] );
disp( ['z-scores: ',num2str( zscoree2),'/',num2str( zscorei2)] );
disp( ['p-values: ',num2str( pkce2),'/',num2str( pkci2)] );

disp( ['High-KC ratio (Exc., k-core>=3): ',num2str(mean(high_KC_BC_index3.exc_ratio{1})) ,' +/- ', num2str(var(high_KC_BC_index3.exc_ratio{1})) ] );
disp( ['High-KC ratio (Inh., k-core>=3): ',num2str(mean(high_KC_BC_index3.inh_ratio{1})) ,' +/- ', num2str(var(high_KC_BC_index3.inh_ratio{1})) ] );
disp( ['z-scores: ',num2str( zscoree3),'/',num2str( zscorei3)] );
disp( ['p-values: ',num2str( pkce3),'/',num2str( pkci3)] );

