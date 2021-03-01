
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program products figure 5(b).
% Figure 5 shows basic properties of neuronal networks, and
% especially 5-(b) shows degree histgrams, histgrams of number of connections.
%
%
% Author:          Masanori Shimono     2018/2019
%        Modified  Motoki Kajiwara      2019/2020
%        Cleaned   Ritsuki Nomura       2020/2021
% contact address: shimonomlab@gmail.com
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;  clear all;

kk = 1;
for data_num = [1,2,3,4,5,6,7];
    
    data_name1  = ['data',num2str(data_num)];
    data_name3 = data_name1;
    
    %%
    str = ['load  ./orig_data/', data_name1 ,'/Aft_Limited/folder_name4;'];     eval(str);
    str = ['load  ./orig_data/', data_name1 ,'/Aft_Limited/degree;'];        eval(str);
    
    degree_all.eiin_sum(kk) = sum(degree.eiin);
    degree_all.eein_sum(kk) = sum(degree.eein);
    degree_all.iein_sum(kk) = sum(degree.iein);
    degree_all.iiin_sum(kk) = sum(degree.iiin);
    
    degree_all.eiout_sum(kk) = sum(degree.eiout);
    degree_all.eeout_sum(kk) = sum(degree.eeout);
    degree_all.ieout_sum(kk) = sum(degree.ieout);
    degree_all.iiout_sum(kk) = sum(degree.iiout);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    range_log = log10([1:2:29]);
    
    str = ['load ./orig_data/', data_name1 ,'/Aft_Limited/deghist_in_EE;'];     eval(str);
    val = log10(deghist_in.Values);  val(isinf(val)) = 0.001*rand(1);  val(val==0) = 0.001*rand(1); deghist_all_EE.in(kk, :) = deghist_in.Values; clear deghist_in
    deghist_log_all_EE.in(kk, :) = val; clear val
    
    str = ['load ./orig_data/', data_name1 ,'/Aft_Limited/deghist_in_EI;'];     eval(str);
    val = log10(deghist_in.Values);  val(isinf(val)) = 0.001*rand(1);  val(val==0) = 0.001*rand(1); deghist_all_EI.in(kk, :) = deghist_in.Values; clear deghist_in
    deghist_log_all_EI.in(kk, :) = val; clear val
    
    str = ['load ./orig_data/', data_name1 ,'/Aft_Limited/deghist_in_IE;'];     eval(str);
    val = log10(deghist_in.Values);  val(isinf(val)) = 0.001*rand(1);  val(val==0) = 0.001*rand(1); deghist_all_IE.in(kk, :) = deghist_in.Values; clear deghist_in
    deghist_log_all_IE.in(kk, :) = val; clear val
    
    str = ['load ./orig_data/', data_name1 ,'/Aft_Limited/deghist_in_II;'];     eval(str);
    val = log10(deghist_in.Values);  val(isinf(val)) = 0.001*rand(1);  val(val==0) = 0.001*rand(1); deghist_all_II.in(kk, :) = deghist_in.Values; clear deghist_in
    deghist_log_all_II.in(kk, :) = val; clear val
    
    str = ['load  ./orig_data/', data_name1 ,'/Aft_Limited/deghist_in;'];     eval(str);
    val = log10(deghist_in.Values);  val(isinf(val)) = 0.001*rand(1);  val(val==0) = 0.001*rand(1); deghist_all.in(kk, :) = deghist_in.Values; clear deghist_in
    deghist_log_all.in(kk, :) = val; clear val
    
    %%
    str = ['load ./orig_data/', data_name1 ,'/Aft_Limited/deghist_out_EE;'];     eval(str);
    val = log10(deghist_out.Values);   val(isinf(val)) = 0.001*rand(1); val(val==0) = 0.001*rand(1);  deghist_all_EE.out(kk, :) = deghist_out.Values; clear deghist_out
    deghist_log_all_EE.out(kk, :) = val; clear val
    
    str = ['load ./orig_data/', data_name1 ,'/Aft_Limited/deghist_out_EI;'];     eval(str);
    val = log10(deghist_out.Values);   val(isinf(val)) = 0.001*rand(1); val(val==0) = 0.001*rand(1);  deghist_all_EI.out(kk, :) = deghist_out.Values; clear deghist_out
    deghist_log_all_EI.out(kk, :) = val; clear val
    
    str = ['load ./orig_data/', data_name1 ,'/Aft_Limited/deghist_out_IE;'];     eval(str);
    val = log10(deghist_out.Values);   val(isinf(val)) = 0.001*rand(1); val(val==0) = 0.001*rand(1);  deghist_all_IE.out(kk, :) = deghist_out.Values; clear deghist_out
    deghist_log_all_IE.out(kk, :) = val; clear val
    
    str = ['load ./orig_data/', data_name1 ,'/Aft_Limited/deghist_out_II;'];     eval(str);
    val = log10(deghist_out.Values);   val(isinf(val)) = 0.001*rand(1); val(val==0) = 0.001*rand(1);  deghist_all_II.out(kk, :) = deghist_out.Values; clear deghist_out
    deghist_log_all_II.out(kk, :) = val; clear val
    
    str = ['load ./orig_data/', data_name1 ,'/Aft_Limited/deghist_out;'];     eval(str);
    val = log10(deghist_out.Values);   val(isinf(val)) = 0.001*rand(1); val(val==0) = 0.001*rand(1);  deghist_all.out(kk, :) = deghist_out.Values; clear deghist_out
    deghist_log_all.out(kk, :) = val; clear val
    
    kk = kk + 1;
    
end

NN = 9;

%% From pdf (combined between pdf_exc + pdf_inh)

figure(181);
aa1 = mean(deghist_all.in(:,1:NN),1);   aa1n = find(aa1~=0);
aa2 = mean(deghist_all.out(:,1:NN),1);  aa2n = find(aa2~=0);

bb1 = std(deghist_all.in(:,1:NN),1)./sqrt(kk-1);   bb1n = find(aa1~=0);
bb2 = std(deghist_all.out(:,1:NN),1)./sqrt(kk-1);  bb2n = find(aa2~=0);

plot(range_log(aa1n), log10(aa1(aa1n)),'m-','Linewidth',1.5);  hold on;
plot(range_log(aa2n), log10(aa2(aa2n)),'g--','Linewidth',1.5);

errorbar(range_log(aa1n), log10(aa1(aa1n)), log10(bb1(aa1n))./sqrt(kk-1),'mo','MarkerSize',5);
errorbar(range_log(aa2n), log10(aa2(aa2n)), log10(bb2(aa2n))./sqrt(kk-1),'go','MarkerSize',5);

%%
y = deghist_all.in(:,1:NN);
for nnl = 1:size(y,1);  plot(range_log(aa1n).*ones(1,length(aa1n)) , log10(y(nnl,aa1n)) ,'m.'); hold on; end

y = deghist_all.out(:,1:NN);
for nnl = 1:size(y,1);  plot(range_log(aa1n).*ones(1,length(aa1n)) , log10(y(nnl,aa1n)) ,'g.'); hold on; end

%%
legend('In degree','Out Degree','In degree (error)','Out Degree (error)','fontsize', 10,'fontname','Arial');
legend(': In degree',': Out degree','fontsize', 15,'fontname','Arial');

title(['Degree histogram'],'fontsize', 18,'fontname','Arial');
xlabel(['Degree [log_{10}]'],'fontsize', 20,'fontname','Arial');
ylabel(['Event Counts [log_{10}]'],'fontsize', 20,'fontname','Arial');
set(gca,'YLim',[-1.2 2.5],'fontsize',18,'fontname','Arial');


clear T1 T2 T3 C1 C2 C3

T1 = table(range_log(aa1n));
T2 = table(deghist_all.in(:,aa1n));
T3 = table(deghist_all.out(:,aa2n));

C1 = table({'log10(degree)'})
C2 = table({'Events (inputs)'})
C3 = table({'Events (outputs)'})

filename = 'degree.xlsx';
writetable(T1,filename,'Sheet',1,'Range','B2')
writetable(T2,filename,'Sheet',2,'Range','B2')
writetable(T3,filename,'Sheet',3,'Range','B2')

writetable(C1,filename,'Sheet',1,'Range','A2')
writetable(C2,filename,'Sheet',2,'Range','A2')
writetable(C3,filename,'Sheet',3,'Range','A2')

%% for answering reviewer1

x1 = degree_all.ieout_sum;
x2 = degree_all.iiout_sum;

[p,zscore_x2] = showstat_vals(x1,x2,'degieout','degiiout');
