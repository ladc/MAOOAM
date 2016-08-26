%
% Short Plot Script for the Lyapunov exponents
% run ./maooam_lyap to get file lyapunov_exponents.dat
%
clear all
close all
%addpath(pathdef)

BLV=true;
FLV=false;
CLV=true;
%Read namelist
nml=read_nml('params.nml','int_params.nml','modeselection.nml');

% Compute total dimension of system
ndim=(sum(nml.AMS(:,1)~=1)*2+sum(nml.AMS(:,1)==1)*3)*2+2*nml.NBOC;


% unit of LEs in 1/day
facLE=nml.F0*24*3600;

sb=dir('compare_correlation_BLV.dat');
fid=fopen('compare_correlation_BLV.dat');
corrBLV=reshape(fread(fid,sb.bytes/8,'real*8'),ndim,sb.bytes/8/ndim);
fclose(fid);

sb=dir('compare_correlation_FLV.dat');
fid=fopen('compare_correlation_FLV.dat');
corrFLV=reshape(fread(fid,sb.bytes/8,'real*8'),ndim,sb.bytes/8/ndim);
fclose(fid);

sb=dir('compare_correlation_CLV.dat');
fid=fopen('compare_correlation_CLV.dat');
corrCLV=reshape(fread(fid,sb.bytes/8,'real*8'),ndim,sb.bytes/8/ndim);
fclose(fid);


sb=dir('compare_diff_BLE.dat');
fid=fopen('compare_diff_BLE.dat');
diffBLE=reshape(fread(fid,sb.bytes/8,'real*8'),ndim,sb.bytes/8/ndim);
fclose(fid);

sb=dir('compare_diff_FLE.dat');
fid=fopen('compare_diff_FLE.dat');
diffFLE=reshape(fread(fid,sb.bytes/8,'real*8'),ndim,sb.bytes/8/ndim);
fclose(fid);

sb=dir('compare_diff_CLE.dat');
fid=fopen('compare_diff_CLE.dat');
diffCLE=reshape(fread(fid,sb.bytes/8,'real*8'),ndim,sb.bytes/8/ndim);
fclose(fid);