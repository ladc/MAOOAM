%
% Short Plot Script for the Lyapunov exponents
% run ./maooam_lyap to get file lyapunov_exponents.dat
%
clear all
close all
addpath(pathdef)

%Read namelist
nml=read_nml('params.nml','int_params.nml','modeselection.nml');

% Compute total dimension of system
ndim=(sum(nml.AMS(:,1)~=1)*2+sum(nml.AMS(:,1)==1)*3)*2+2*nml.NBOC;

% Open Lyapunov Exponent File and determsizeine size of output
s=dir('lyapunov_exponents.dat');
% Length Of Output (number of writeouts)
n_len=floor(s.bytes/8/ndim);
%Read Data, Reshape into 'common sense' form and close file
fid=fopen('lyapunov_exponents.dat');
lle=fread(fid,ndim*n_len,'real*8');
fclose(fid);
lle=reshape(lle,ndim,n_len);
lle_direct=importdata('lyapunov_exponents_div.dat');
%
[BLE_SV,VAR_BLE_SV,FLE_SV,VAR_FLE_SV,CLE_SV,VAR_CLE_SV] = read_SV;

% unit of LEs in 1/day
facLE=nml.F0*24*3600;

% plot of spectrum
figLE=figure;
subplot(2,1,1)
plot(mean(lle,2)*facLE,'+');
hold on;
plot(BLE_SV,'.');
plot(mean(lle_direct)*facLE,'o');
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
subplot(2,1,2)
plot(var(lle,[],2)*facLE^2,'+');
hold on;
plot(VAR_BLE_SV,'.');
plot(var(lle_direct,[],1)*facLE^2,'o');
title('VAR [1/day]')
xlabel('Lyapunov exponent number')
ylabel('[1/day]')
legend({'MAOOAM' 'Paper S&V' 'DIV'})
export_fig (figLE,'result_LE.png');

figts=figure;
plot(lle(1,:)*facLE);
hold on;
plot(lle_direct*facLE,'--')
title('DIV vs QR LE 1')
export_fig(figts,'result_ts.png');
