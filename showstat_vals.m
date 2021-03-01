function [p,zscore_x2] = showstat_vals(x1,x2,name1,name2)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function write down several key statistical quantities
% of Wilcoxon signed rank test.
%
% Author:          Masanori Shimono     2020
%        Modified  Motoki Kajiwara      2020/2021
%        Cleaned   Ritsuki Nomura       2021
% contact address: shimonomlab@gmail.com
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[p, h] = signrank( x1 , x2 );% ,'tail','right');

zscore_x1 = ( mean(x1) - mean(x2))/(var(x2)/sqrt(7));
zscore_x2 = ( mean(x1) - mean(x2))/(var(x1)/sqrt(7));
         
disp(['------------------------------------------------------']);
disp([name1,':  ',num2str(mean(x1)),' +-',num2str(std(x1))]);
disp([name2,':  ',num2str(mean(x2)),' +-',num2str(std(x2))]);
disp( ['   z-scores: ',num2str( zscore_x1),'/',num2str( zscore_x2)] );
disp( ['   p-values: ',num2str(p)] );
    
x3= x1./(x1+x2);
disp([name1,'/(',name1,'+',name2,'):  ',num2str(mean(x3)),' +-',num2str(std(x3))]);