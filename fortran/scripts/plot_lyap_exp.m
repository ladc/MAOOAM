%
% Short Plot Script for the Lyapunov exponents
% run ./maooam_lyap to get file lyapunov_exponents.dat
%
clear all
close all
addpath(pathdef)

%Read namelist
nml=read_nml('params.nml','int_params.nml','modeselection.nml');
blle=0;
clle=0;
% Compute total dimension of system
ndim=(sum(nml.AMS(:,1)~=1)*2+sum(nml.AMS(:,1)==1)*3)*2+2*nml.NBOC;

% Open Lyapunov Exponent File and determsizeine size of output
sb=dir('BLV_exp.dat');
sf=dir('FLV_exp.dat');
sc=dir('CLV_exp.dat');
% Length Of Output (number of writeouts)
n_lenb=floor(sb.bytes/8/ndim);
n_lenf=floor(sf.bytes/8/ndim);
n_lenc=floor(sc.bytes/8/ndim);
%Read Data, Reshape into 'common sense' form and close file
fid=fopen('BLV_exp.dat');
blle=fread(fid,ndim*n_lenb,'real*8');
fclose(fid);
blle=reshape(blle,ndim,n_lenb);
fid=fopen('FLV_exp.dat');
flle=fread(fid,ndim*n_lenf,'real*8');
fclose(fid);
flle=reshape(flle,ndim,n_lenf);
fid=fopen('CLV_exp.dat');
clle=fread(fid,ndim*n_lenc,'real*8');
fclose(fid);
clle=reshape(clle,ndim,n_lenc);

[BLE_SV,VAR_BLE_SV,FLE_SV,VAR_FLE_SV,CLE_SV,VAR_CLE_SV] = read_SV;

% unit of LEs in 1/day
facLE=nml.F0*24*3600;


A=1;
B=size(blle,2);
%plot of spectrum
figBLE=figure;
subplot(2,1,1)
plot(mean(blle(:,A:B),2)*facLE,'+');
hold on;
plot(BLE_SV,'.');
% plot(mean(lle_direct)*facLE,'o');
title('Lyapunov Spectrum')
xlabel('Lyapunov exponent number')
ylabel('[1/day]')
annotation(gcf,'textbox',...
    [0.175 0.6 0.2 0.178],...
    'String',{['DT=',num2str(nml.DT)],...
              ['T\_RUN=',num2str(nml.T_RUN)],...
              ['T\_TRANS=',num2str(nml.T_TRANS)],...
              ['RESCALING_TIME=',num2str(nml.RESCALING_TIME)]},...
    'FitBoxToText','off','edgecolor','none');
legend({'MAOOAM' 'Paper S&V'})
title('BLV LE')

subplot(2,1,2)
plot(var(blle(:,A:B),[],2)*facLE^2,'+');
hold on;
plot(VAR_BLE_SV,'.');
% plot(var(lle_direct,[],1)*facLE^2,'o');
title('VAR [1/day]')
xlabel('Lyapunov exponent number')
ylabel('[1/day]')
legend({'MAOOAM' 'Paper S&V'})
title('VAR BLV LE')
export_fig blv.png

% 
A=find(flle(1,:)~=0,1,'first');
B=size(flle,2);

figFLE=figure;

subplot(2,1,1)
plot(mean(flle(:,A:B),2)*facLE,'+');
hold on;
plot(FLE_SV,'.');
title('Lyapunov Spectrum')
xlabel('Lyapunov exponent number')
ylabel('[1/day]')
annotation(gcf,'textbox',...
    [0.175 0.6 0.2 0.178],...
    'String',{['DT=',num2str(nml.DT)],...
              ['T\_RUN=',num2str(nml.T_RUN)],...
              ['T\_TRANS=',num2str(nml.T_TRANS)],...
              ['RESCALING_TIME=',num2str(nml.RESCALING_TIME)]},...
    'FitBoxToText','off','edgecolor','none');
legend({'MAOOAM' 'Paper S&V' 'DIV'})
title('FLV LE')

subplot(2,1,2)
plot(var(flle(:,A:B),[],2)*facLE^2,'+');
hold on;
plot(VAR_FLE_SV,'.');
title('VAR [1/day]')
xlabel('Lyapunov exponent number')
ylabel('[1/day]')
legend({'MAOOAM' 'Paper S&V'})
title('VAR FLV LE')
export_fig flv.png

A=find(clle(1,:)~=0,1,'first');
B=size(clle,2);

figCLE=figure;
subplot(2,1,1)
plot(mean(clle(:,A:B),2)*facLE,'+');
hold on;
plot(CLE_SV,'.');
title('Lyapunov Spectrum')
xlabel('Lyapunov exponent number')
ylabel('[1/day]')
annotation(gcf,'textbox',...
    [0.175 0.6 0.2 0.178],...
    'String',{['DT=',num2str(nml.DT)],...
              ['T\_RUN=',num2str(nml.T_RUN)],...
              ['T\_TRANS=',num2str(nml.T_TRANS)],...
              ['RESCALING_TIME=',num2str(nml.RESCALING_TIME)]},...
    'FitBoxToText','off','edgecolor','none');
legend({'MAOOAM' 'Paper S&V' 'DIV'})
title('CLV LE')
subplot(2,1,2)
plot(var(clle(:,A:B),[],2)*facLE^2,'+');
hold on;
plot(VAR_CLE_SV,'.');
title('VAR [1/day]')
xlabel('Lyapunov exponent number')
ylabel('[1/day]')
legend({'MAOOAM' 'Paper S&V'})
title('VAR CLV LE')
export_fig clv.png

% export_fig (figLE,'result_LE.png');
% 
% figts=figure;
% plot(blle(1,:)*facLE);
% hold on;
% plot(lle_direct*facLE,'--')
% title('DIV vs QR LE 1')
% export_fig(figts,'result_ts.png');
