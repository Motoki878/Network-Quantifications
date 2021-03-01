
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program products data set for plotting the Venn diagram
% to compare between high KC nodes and FV nodes.
%
%Author:          Masanori Shimono     2018/2019
%        Modified  Motoki Kajiwara      2019/2020
%        Cleaned  Ritsuki Nomura        2020/2021
% contact address: shimonomlab@gmail.com
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nodeall = [102,363,250,172,449,387,219];
FVSall  = [  9,  5,  8,  6, 7,  10, 10];

%% only KC max core
KCmax   = [ 14, 67, 55, 12, 8, 210, 51];
FVScomm = [  5,  5,  6,  6, 1,  10,  9];

fkm = 100*mean(FVScomm./KCmax);
ffm = 100*mean(FVScomm./FVSall);
fks = 100*std(FVScomm./KCmax);
ffs = 100*std(FVScomm./FVSall);

[p, h] = signrank( FVScomm./KCmax, FVScomm./FVSall);

zscorefvs = ( mean(FVScomm./KCmax)- mean(FVScomm./FVSall))/(var(FVScomm./FVSall)/sqrt(7));
zscorekc  = ( mean(FVScomm./KCmax)- mean(FVScomm./FVSall))/(var(FVScomm./KCmax)/sqrt(7));

disp(['------------------------------------------------------']);
disp(['FVScomm:  ',num2str(mean(FVScomm)),' +-',num2str(std(FVScomm))]);
disp(['KCmax:  ',num2str(mean(KCmax)),' +-',num2str(std(KCmax))]);
disp(['FVSall: ',num2str(mean(FVSall)),' +-',num2str(std(FVSall))]);
disp(['    FVScomm/KCmax:  ',num2str(mean(FVScomm./KCmax)),' +-',num2str(std(FVScomm./KCmax))]);
disp(['    FVScomm/FVSall: ',num2str(mean(FVScomm./FVSall)),' +-',num2str(std(FVScomm./FVSall))]);
disp( ['   z-scores: ',num2str( zscorefvs),'/',num2str( zscorekc)] );
disp( ['   p-values: ',num2str(p)] );
disp(['FVScomm per KCmax:  ',num2str(fkm),' +-',num2str(fks),'[%]']);
disp(['FVScomm per FVSall: ',num2str(ffm),' +-',num2str(ffs),'[%]']);
disp(['------------------------------------------------------']);

%% till second KC core
KCmax2   = [14+39, 67+90, 55+72, 12+31, 8+274, 210+54, 51+52];
FVScomm2 = [ 5+ 4,  5+ 0,  6+ 2,  6+ 0,  1+ 6,  10+ 0,  9+ 1];

fkm2 = 100*mean(FVScomm2./KCmax2);
ffm2 = 100*mean(FVScomm2./FVSall);
fks2 = 100*std(FVScomm2./KCmax2);
ffs2 = 100*std(FVScomm2./FVSall);

[p2, h2] = signrank( FVScomm2./KCmax2, FVScomm2./FVSall);

zscorefvs2 = ( mean(FVScomm2./KCmax2)- mean(FVScomm2./FVSall))/(var(FVScomm2./FVSall)/sqrt(7));
zscorekc2  = ( mean(FVScomm2./KCmax2)- mean(FVScomm2./FVSall))/(var(FVScomm2./KCmax2)/sqrt(7));

disp(['FVScomm2:  ',num2str(mean(FVScomm2)),' +-',num2str(std(FVScomm2))]);
disp(['KCmax2:  ',num2str(mean(KCmax2)),' +-',num2str(std(KCmax2))]);
disp(['    FVScomm2/KCmax:  ',num2str(mean(FVScomm2./KCmax2)),' +-',num2str(std(FVScomm2./KCmax2))]);
disp(['    FVScomm2/FVSall: ',num2str(mean(FVScomm2./FVSall)),' +-',num2str(std(FVScomm2./FVSall))]);
disp( ['   z-scores: ',num2str( zscorefvs2),'/',num2str( zscorekc2)] );
disp( ['   p-values: ',num2str(p2)] );
disp(['FVScomm2 per KCmax2: ',num2str(fkm2),' +-',num2str(fks2),'[%]']);
disp(['FVScomm2 per FVSall   :    ',num2str(ffm2),' +-',num2str(ffs2),'[%]']);
disp(['------------------------------------------------------']);

