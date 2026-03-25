clc;clear all;close all;
global COCHBA;
load aud24
paras = [8 4 -2 -1];
BP = 0;  % I like the reconstructions with BP=0 but the performance is  not generalizable because of the DC
fullt = 0;
fullx = 0;
para1=[paras fullt fullx];
para2 = [0.900 fullt fullx BP];
xcome = unitseq(auread('_come.au'));
y = wav2aud(xcome,paras);
rv = [2 4 6 8 10 16 18 20 ];
sv=[2 4 6 8] ;
yhat = 0;
for sv= [2:2:8]
 cr = aud2cor(y, para1,rv, sv ,'_come');
[yh, para1, rv, sv, HH] = cor2aud('_come', cr, para2);
yhat = yhat + aud_fix(yh);
end
figure(1)
imagesc(y');axis xy;
figure(2)
imagesc((yhat)');axis xy
