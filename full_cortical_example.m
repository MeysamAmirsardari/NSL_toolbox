clc;clear;close all;
global COCHBA;
load aud24
paras = [8 4 -2 -1];
xcome = unitseq(auread('_come.au'));
y = wav2aud(xcome,paras);
fcorname = '_come.cor';
rv = 2.^(1:5);
sv = 2.^(-2:3);

cr = aud2cor(y, paras, rv, sv, fcorname);
%% SPECTROGRAM RECONSTRUCTION
yh = cor2aud(fcorname);
yh = aud_fix(yh);
