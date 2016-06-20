%
% Short Plot Script for the Lyapunov exponents
% run ./maooam_lyap to get file lyapunov_exponents.dat
%

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


% unit of LEs in 1/day
facLE=nml.F0*24*3600;

% plot of spectrum
figure;
plot(mean(lle,2)*facLE);
title('Lyapunov Spectrum')
xlabel('Lyapunov exponent number')
ylabel('[1/day]')
annotation(gcf,'textbox',...
    [0.6 0.6 0.2 0.3],...
    'String',{['DT=',num2str(nml.DT)],...
              ['T\_RUN=',num2str(nml.T_RUN)],...
              ['T\_TRANS=',num2str(nml.T_TRANS)],...
              ['TW=',num2str(nml.TW)]},...
    'FitBoxToText','off','edgecolor','none');
