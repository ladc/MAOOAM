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
offsets=[0 1];
%Read namelist
nml=read_nml('params.nml','int_params.nml','modeselection.nml');

% Compute total dimension of system
ndim=(sum(nml.AMS(:,1)~=1)*2+sum(nml.AMS(:,1)==1)*3)*2+2*nml.NBOC;

% Open Lyapunov Exponent File and determsizeine size of output
if BLV; sb=dir('BLV_exp_0.dat'); end;
if FLV; sf=dir('FLV_exp_0.dat'); end;
if CLV; sc=dir('CLV_exp_0.dat'); end;
% Length Of Output (number of writeouts)
if BLV; n_lenb=floor(sb.bytes/8/ndim); end;
if FLV; n_lenf=floor(sf.bytes/8/ndim); end;
if CLV; n_lenc=floor(sc.bytes/8/ndim); end;
%Read Data, Reshape into 'common sense' form and close file
o=0;
for offset=offsets
    o=o+1;
    if BLV
        fid=fopen(['BLV_exp_',num2str(offset),'.dat']);
        blle{o}=fread(fid,ndim*n_lenb,'real*8');
        fclose(fid);
        blle{o}=reshape(blle{o},ndim,n_lenb);
    end
    if FLV
        fid=fopen(['FLV_exp_',num2str(offset),'.dat']);
        flle{o}=fread(fid,ndim*n_lenf,'real*8');
        fclose(fid);
        flle{o}=reshape(flle{o},ndim,n_lenf);
    end
    if CLV
        fid=fopen(['CLV_exp_',num2str(offset),'.dat']);
        clle{o}=fread(fid,ndim*n_lenc,'real*8');
        fclose(fid);
        clle{o}=reshape(clle{o},ndim,n_lenc);
    end
end

% unit of LEs in 1/day
facLE=nml.F0*24*3600;
if BLV
    figure;
    plot(blle{1}(1,:))
    hold on;
    plot(blle{2}(1,:),'--')
    title('Comparison of time series of first BLE')
    %
    % Comaprison of convergence times
    %
    LE=mean([blle{1} blle{2}],2);
    d=abs(diff(LE));
    dummy=abs(blle{1}(:,:)-blle{2}(:,:))./abs(blle{2}(:,:));
    distance=[abs(LE(1)-LE(2)) arrayfun(@(x,y) min(x,y),d(1:end-1),d(2:end))' abs(LE(end)-LE(end-1))];
    f=@(x,p) find(dummy(x,:)<p/100,1);
    
    
    i=0;
    for p=10.^[-2 -3 -4 -5 -6 -7 ]
        i=i+1
        for l=1:ndim
            dummy=f(l,p)*nml.DT;
            if isempty(dummy); dummy=-n_lenb;end;
            dropoff(l,i)=dummy;
        end
    end
    ps=10.^[-2 -3 -4 -5 -6 -7 ];
    figC=figure;
    for p=10.^[-2 -3 -4 -5 -6 -7 ]
        subplot(3,2,-log10(p)-1)
        plot(1:36,dropoff(:,p==ps))
        hold on;
        plot(1:36,1./distance)
        xlabel('CLVs','fontsize',3)
        ylabel('Convergence Time','fontsize',3)
        legend({sprintf( '%s\n%s', 'Time until difference ',  ['is below ',num2str(p,'%1.0e')] ) ...
            sprintf( '%s\n%s', 'Inverse of smallest', 'distance to neighbor LEs' )},'fontsize',3)
        set(gca,'fontsize',3)
    end
    % export_fig('convergenceBLE','-png','-pdf','-m6')
end

if CLV
    figure;
    plot(clle{1}(1,:))
    hold on;
    plot(clle{2}(1,:),'--')
    title('Comparison of time series of first CLE')
    %
    % Comaprison of convergence times
    %
    LE=mean([clle{1} clle{2}],2);
    d=abs(diff(LE));
    dummy=abs(clle{1}(:,:)-clle{2}(:,:))./abs(clle{2}(:,:));
    distance=[abs(LE(1)-LE(2)) arrayfun(@(x,y) min(x,y),d(1:end-1),d(2:end))' abs(LE(end)-LE(end-1))];
    f=@(x,p) find(dummy(x,:)<p/100,1);
    f2=@(x,p) find(dummy(x,:)<p/100,1,'last');
    
    i=0
    for p=10.^[-2 -3 -4 -5 -6 -7 ]
        i=i+1
        for l=1:ndim
            dummy=f(l,p)*nml.DT;
            if isempty(dummy); dummy=-n_lenc;end;
            dropoff(l,i)=dummy;
            dummy=f2(l,p)*nml.DT;
            if isempty(dummy); dummy=-n_lenc;end;
            dropon(l,i)=dummy;
        end
    end
    
    ps=10.^[-2 -3 -4 -5 -6 -7 ];
    figC=figure;
    for p=10.^[-2 -3 -4 -5 -6 -7 ]
        subplot(3,2,-log10(p)-1)
        plot(1:36,dropoff(:,p==ps))
        hold on;
        plot(1:36,dropon(:,p==ps))
        plot(1:36,1./distance)
        xlabel('CLVs','fontsize',3)
        ylabel('Convergence Time','fontsize',3)
        legend({sprintf( '%s\n%s', 'Forward Time until difference ',  ['is below ',num2str(p,'%1.0e')] ) ...
            sprintf( '%s\n%s', 'Backward Time until difference ',  ['is below ',num2str(p,'%1.0e')] ) ...
            sprintf( '%s\n%s', 'Inverse of smallest', 'distance to neighbor LEs' )},'fontsize',3)
        set(gca,'fontsize',3)
    end
    export_fig('convergenceCLE','-png','-pdf','-m6')
end
