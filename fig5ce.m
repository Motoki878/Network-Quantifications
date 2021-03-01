
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program products figure 5(c)(e).
% Figure 5 shows basic properties of neuronal networks, especially
% (c) shows histgrams of firing rate for excitatory and inhibitory neurons.
% (e) shows histgrams of connectivity strenghs for excitatory and inhibitory neurons.
% 
% Author:          Masanori Shimono     2018/2019
%        Modified  Motoki Kajiwara      2019/2020
%        Cleaned   Ritsuki Nomura       2020/2021
% contact address: shimonomlab@gmail.com
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

kk = 1;


for data_num = [1,2,3,4,5,6,7];
    
    data_name1  = ['data',num2str(data_num)];
    data_name3  = data_name1;
    
    %%
    str = ['load  ./orig_data/', data_name1 ,'/Aft_Limited/folder_name4;'];     eval(str);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    range_FR=[-6:0.1:4];
    str = ['load  ./orig_data/', data_name1 ,'/Aft_Limited/hist_FR;'];     eval(str);
    str = ['load ./orig_data/', data_name1 ,'/Aft_Limited/hist_FR_raw;'];     eval(str);
    val = hist_FR.Values;  val(isinf(val)) = 0;
    hist_FR_all.all(kk, :) = val;
    hist_FR_nor_all.all(kk, :) = (val-mean(val))./std(val); clear hist_FR
    
    
    if kk == 1
        hist_FR_raw_all = hist_FR_raw;
    else
        hist_FR_raw_all = [ hist_FR_raw_all , hist_FR_raw  ];
    end
    
    %%
    str = ['load  ./orig_data/', data_name1 ,'/Aft_Limited/hist_FR_exc;'];     eval(str);
    str = ['load ./orig_data/', data_name1 ,'/Aft_Limited/hist_FR_exc_raw;'];     eval(str);
    val = hist_FR_exc.Values;  val(isinf(val)) = 0;
    hist_FR_all.exc(kk, :) = val;
    hist_FR_nor_all.exc(kk, :) = (val-mean(val))./std(val); clear hist_FR_exc
    
    if kk == 1
        hist_FR_exc_raw_all = hist_FR_exc_raw;
    else
        hist_FR_exc_raw_all = [ hist_FR_exc_raw_all , hist_FR_exc_raw  ];
    end
    
    %%
    str = ['load  ./orig_data/', data_name1 ,'/Aft_Limited/hist_FR_inh;'];     eval(str);
    str = ['load ./orig_data/', data_name1 ,'/Aft_Limited/hist_FR_inh_raw;'];     eval(str);
    val = hist_FR_inh.Values;  val(isinf(val)) = 0;
    hist_FR_all.inh(kk, :) = val;
    hist_FR_nor_all.inh(kk, :) = (val-mean(val))./std(val); clear hist_FR_inh
    
    if kk == 1
        hist_FR_inh_raw_all = hist_FR_inh_raw;
    else
        hist_FR_inh_raw_all = [ hist_FR_inh_raw_all , hist_FR_inh_raw  ];
    end
    
    %%
    range_wgt=[-6:0.1/2.5:-4];
    
    str = ['load ./orig_data/', data_name1 ,'/Aft_Limited/hist_wgt_conn;'];     eval(str);
    str = ['load ./orig_data/', data_name1 ,'/Aft_Limited/hist_wgt_conn_raw;'];     eval(str);
    val = hist_wgt_conn.Values;  val(isinf(val)) = 0;
    hist_wgt_conn_all.all(kk, :) = val;
    
    if kk == 1
        hist_wgt_conn_raw_all = hist_wgt_conn_raw;
    else
        hist_wgt_conn_raw_all = [hist_wgt_conn_raw_all; hist_wgt_conn_raw];
    end
    
    %%
    str = ['load ./orig_data/', data_name1 ,'/Aft_Limited/hist_wgt_conn_exc;'];     eval(str);
    str = ['load ./orig_data/', data_name1 ,'/Aft_Limited/hist_wgt_conn_exc_raw;'];     eval(str);
    val = hist_wgt_conn_exc.Values;  val(isinf(val)) = 0;
    hist_wgt_conn_all.exc(kk, :) = val;
    
    if kk == 1
        hist_wgt_conn_exc_raw_all = hist_wgt_conn_exc_raw;
    else
        hist_wgt_conn_exc_raw_all = [hist_wgt_conn_exc_raw_all; hist_wgt_conn_exc_raw];
    end
    
    %%
    str = ['load ./orig_data/', data_name1 ,'/Aft_Limited/hist_wgt_conn_inh;'];     eval(str);
    str = ['load ./orig_data/', data_name1 ,'/Aft_Limited/hist_wgt_conn_inh_raw;'];     eval(str);
    val = hist_wgt_conn_inh.Values;  val(isinf(val)) = 0;
    hist_wgt_conn_all.inh(kk, :) = val;
    
    if kk == 1
        hist_wgt_conn_inh_raw_all = hist_wgt_conn_inh_raw;
    else
        hist_wgt_conn_inh_raw_all = [hist_wgt_conn_inh_raw_all; hist_wgt_conn_inh_raw];
    end
    
    kk = kk + 1;
    
end

[h,p,adstat,cv] = adtest(hist_FR_exc_raw_all  ,'Alpha',0.01);  h_gfit.FR_exc = h;   p_gfit.FR_exc = p; adstat_gfit.FR_exc = adstat;
[h,p,adstat,cv] = adtest(hist_FR_inh_raw_all  ,'Alpha',0.01);   h_gfit.FR_inh = h;   p_gfit.FR_inh = p;  adstat_gfit.FR_inh = adstat;

[h,p,adstat,cv] = adtest(hist_wgt_conn_exc_raw_all  ,'Alpha',0.01); h_gfit.wgt_exc = h;   p_gfit.wgt_exc = p;  adstat_gfit.wgt_exc = adstat;
[h,p,adstat,cv] = adtest(hist_wgt_conn_inh_raw_all  ,'Alpha',0.01);  h_gfit.wgt_inh = h;  p_gfit.wgt_inh = p;  adstat_gfit.wgt_inh = adstat;

hist_FR_all.all = hist_FR_all.exc + hist_FR_all.inh;
hist_FR_all.inh = 10*hist_FR_all.inh ;

%%
clear a
yy  = mean(hist_FR_all.exc);
for iin = 1:size(hist_FR_all.exc,1)
    [a0]=find(hist_FR_all.exc(iin,:)==max(hist_FR_all.exc(iin,:)));
    a(iin) = a0(1);
end


stds = (std(hist_FR_all.exc))/sqrt(kk-1);
y   = mean(hist_FR_all.exc,1);

x = range_FR([1:end-1]);
f = fit(x.',y.','gauss1');

plot(f,'r-');   hold on;

errorbar(x, y, stds,'ro');

%%
clear a
yy  = mean(hist_FR_all.inh);

x = range_FR([1:end-1]);
hist_FR_all.inh(:,x<-0.4) = 0;

stds = (std(hist_FR_all.inh))/sqrt(kk-1);
y   = mean(hist_FR_all.inh,1);

f = fit(x.',y.','gauss1');
plot(f,'b-');  hold on;

errorbar(x, y, stds,'bo');
set(gca,'YLim',[0 100],'fontsize',10,'fontname','Arial');

%%
dx = x(2)-x(1);
for nnl = 1:length(x)
    plot(x(nnl)*ones(1,size(hist_FR_all.exc,1)) ,hist_FR_all.exc(:,nnl) ,'r.'); hold on;
end

dx = x(2)-x(1);
for nnl = 1:length(x)
    plot(x(nnl)*ones(1,size(hist_FR_all.inh,1)) ,hist_FR_all.inh(:,nnl) ,'b.'); hold on;
end

legend('Exc (Gauss fit)','Exc (Raw data)','Inh (Gauss fit)','Inh (Raw data*10)','fontsize',10,'fontname','Arial');

%%

xlabel(['Firing Rate'],'fontsize', 20,'fontname','Arial');
ylabel(['Event Counts '],'fontsize', 20,'fontname','Arial');
set(gca,'XLim',[-1.5 2.0],'fontsize',18,'fontname','Arial');

title(['p(exc., inh.)= [', num2str(mean(p_gfit.FR_exc)) ,',', num2str(mean(p_gfit.FR_inh)) ,']']);

%%
figure(186);
clear a
yy  = mean(hist_wgt_conn_all.exc);
for iin = 1:size(hist_wgt_conn_all.exc,1)
    [a0]=find(hist_wgt_conn_all.exc(iin,:)==max(hist_wgt_conn_all.exc(iin,:)));
    a(iin) = a0(1);
end


stds = (std(hist_wgt_conn_all.exc))/sqrt(kk-1);
y   = mean(hist_wgt_conn_all.exc,1);

x = range_wgt([1:end-1]);
f = fit(x.',y.','gauss1');
plot(f,'r-');  hold on;
errorbar(x, y, stds,'ro');

dx = x(2)-x(1);
for nnl = 1:length(x)
    plot(x(nnl)*ones(1,size(hist_wgt_conn_all.exc,1)) ,hist_wgt_conn_all.exc(:,nnl) ,'r.'); hold on;
end

%%
clear a
yy  = mean(hist_wgt_conn_all.inh);
for iin = 1:size(hist_wgt_conn_all.inh,1)
    [a0]=find(hist_wgt_conn_all.inh(iin,:)==max(hist_wgt_conn_all.inh(iin,:)));
    a(iin) = a0(1);
end


stds = (std(hist_wgt_conn_all.inh))/sqrt(kk-1);
y   = mean(hist_wgt_conn_all.inh,1);

x = range_wgt([1:end-1]); %
f = fit(x.',y.','gauss1');
plot(f,'b-');  hold on;    hold on;

errorbar(x, y, stds,'bo');

dx = x(2)-x(1);
for nnl = 1:length(x)
    plot(x(nnl)*ones(1,size(hist_wgt_conn_all.inh,1)) ,hist_wgt_conn_all.inh(:,nnl) ,'b.'); hold on;
end

xlabel(['Conncevity Strength'],'fontsize', 20,'fontname','Arial');
ylabel(['Event Counts'],'fontsize', 20,'fontname','Arial');

set(gca,'XLim',[-6.0 -4.5],'fontsize',18,'fontname','Arial');

set(gca,'YLim',[0 55],'fontsize',10,'fontname','Arial');

legend('Exc (Gauss fit)','Exc (Raw data)','Inh (Gauss fit)','Inh (Raw data)','fontsize',10,'fontname','Arial');

title(['p(exc., inh.)= [', num2str(mean(p_gfit.wgt_exc)) ,',', num2str(mean(p_gfit.wgt_inh)) ,']']);

%%

T1 = table(range_FR([1:end-1]));
T2 = table(hist_FR_all.exc);
T3 = table(hist_FR_all.inh);

C1 = table({'FR'})
C2 = table({'Events (excitatory)'})
C3 = table({'Events (inhibitory)'})

filename = 'FR.xlsx';
writetable(T1,filename,'Sheet',1,'Range','B2')
writetable(T2,filename,'Sheet',2,'Range','B2')
writetable(T3,filename,'Sheet',3,'Range','B2')

writetable(C1,filename,'Sheet',1,'Range','A2')
writetable(C2,filename,'Sheet',2,'Range','A2')
writetable(C3,filename,'Sheet',3,'Range','A2')

%%

T1 = table(range_wgt([1:end-1]));
T2 = table(hist_wgt_conn_all.exc);
T3 = table(hist_wgt_conn_all.inh);

C1 = table({'Weights'});
C2 = table({'Events (excitatory)'});
C3 = table({'Events (inhibitory)'});

filename = 'weights.xlsx';
writetable(T1,filename,'Sheet',1,'Range','B2');
writetable(T2,filename,'Sheet',2,'Range','B2');
writetable(T3,filename,'Sheet',3,'Range','B2');

writetable(C1,filename,'Sheet',1,'Range','A2');
writetable(C2,filename,'Sheet',2,'Range','A2');
writetable(C3,filename,'Sheet',3,'Range','A2');

