fileobs = '/Users/florentbrient/Dropbox/CMIP5/CERES2/albcld/ERSST/albcld_ts_-30_30_0_360_ocean_slct_50000_ERSST_historical_step1_TEST.mat';
ff      = open(fileobs);
pathout = '/Users/florentbrient/Dropbox/GitHub/Cloud-variability-time-frequency/data/observations/';

% Albedo evolution
filein1 = [pathout,'albcld_ceres.txt'];
albcld  = ff.ev_a'*100.0;
fileID = fopen(filein1,'w');
fprintf(fileID,'%9s\n','Albcld (%)');
fprintf(fileID,'%9.4f\n',albcld);
fclose(fileID);

% SST evolution
filein2 = [pathout,'sst_ersst.txt'];
sst     = ff.ev_b';
fileID = fopen(filein2,'w');
fprintf(fileID,'%9s\n',char(['SST (C)']));
fprintf(fileID,'%9.4f\n',sst);
fclose(fileID);

% Outputs
% iimax    = 5; %(deseason, intra, season, inter, decadal)
fileout = [pathout,'results_obs.txt'];
timefreq= {'deseason';'intra';'season';'inter';'decadal'}';
corr    = ff.corr_FFT';
slope1  = ff.slope_FFT'*100.0;
slope2  = ff.b2_FFT'*100.0;
results = [corr,slope1,slope2];
fileID = fopen(fileout,'w');
fprintf(fileID,'%13s %13s %13s %13s\n','Frequency','Corr coef','OLS regress',' robust regress');
for ij=1:length(timefreq)
  fprintf(fileID,'%13s %13.4f %13.4f %13.4f\n',char(timefreq(ij)),results(ij,:));
end
fclose(fileID);



