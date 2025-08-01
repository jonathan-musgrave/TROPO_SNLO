% plots pump from sweep;
if options.SweepDir == 1
SweepDir_str = 'normal';
else
SweepDir_str = 'reverse';
end
LW = 2;
tunits = 1e-12; tunits_str = 'ps';
wunits = 1e12*2*pi; wunits_str = 'THz';
CL = [-60,0];
    
Ip = abs(LP).^2;
Is = abs(LS).^2;
Ii = abs(LI).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Figure 1 Intracavity Power %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(fig_num+1);clf;
hold on;
plot(detune_ar./alpha,sum(Ip,2),'Color',[0,0,0.7],linewidth = LW)
plot(detune_ar./alpha,sum(Is,2),'Color',[0,0.7,0],linewidth = LW)
plot(detune_ar./alpha,sum(Ii,2),'Color',[0.7,0,0],linewidth = LW)
xlabel(strcat('\Delta_',options.SweepResonance,' (\alpha_p)'))
set(gca,'XDir',SweepDir_str)
ylabel('Power (W)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Figure 2 Pump Signal And Idler Normalized temporal Evolution %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:Nrt
    IpN = Ip(i,:)./max(Ip(i,:));
    IsN = Is(i,:)./max(Is(i,:));
    IiN = Ii(i,:)./max(Ii(i,:));
end
figure(fig_num+2);clf;
s1 = subplot(1,3,1);
hold on;
imagesc(t./tunits,detune_ar./alpha,IpN);colorbar;
ylabel(strcat('\Delta_',options.SweepResonance,' (\alpha_p)'))
set(gca,'YDir',SweepDir_str)
xlabel(['time (',tunits_str,')'])
title('Pump Evolution')
s2 = subplot(1,3,2);
hold on;
imagesc(t./tunits,detune_ar./alpha,IsN);colorbar;
ylabel(strcat('\Delta_',options.SweepResonance,' (\alpha_p)'))
set(gca,'YDir',SweepDir_str)
xlabel(['time (',tunits_str,')'])
title('Signal Evolution')
s3 = subplot(1,3,3);
hold on;
imagesc(t./tunits,detune_ar./alpha,IiN);colorbar;
ylabel(strcat('\Delta_',options.SweepResonance,' (\alpha_p)'))
set(gca,'YDir',SweepDir_str)
xlabel(['time (',tunits_str,')'])
title('Idler Evolution')
linkaxes([s1,s2,s3],'xy')
xlim([min(t),max(t)]./tunits)
ylim([min(detune_ar),max(detune_ar)]./alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 3 Pump Signal And Idler Normalized frequency Evolution %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cp = max(max(Ip)); Cs = max(max(Is)); Ci = max(max(Ii)); 
Iwp = fftshift(fft(fftshift(LP,2)));
Iws = fftshift(fft(fftshift(LS,2)));
Iwi = fftshift(fft(fftshift(LI,2)));
for i = 1:Nrt
    IpN = 10*log10(abs(Iwp(i,:)./max(Iwp(i,:))).^2);
    IsN = 10*log10(abs(Iws(i,:)./max(Iws(i,:))).^2);
    IiN = 10*log10(abs(Iwi(i,:)./max(Iwi(i,:))).^2);
end
figure(fig_num+3);clf;
s1 = subplot(1,3,1);
hold on;
imagesc(w./wunits,detune_ar./alpha,IpN);c = colorbar;c.Title.String = 'db';
ylabel(strcat('\Delta_',options.SweepResonance,' (\alpha_p)'))
set(gca,'YDir',SweepDir_str)
xlabel(['frequnecy (',wunits_str,')'])
title('Pump Evolution')
set(gca,'CLIM',CL)
s2 = subplot(1,3,2);
hold on;
imagesc(w./wunits,detune_ar./alpha,IsN);c = colorbar;c.Title.String = 'db';
ylabel(strcat('\Delta_',options.SweepResonance,' (\alpha_p)'))
set(gca,'YDir',SweepDir_str)
xlabel(['frequnecy (',wunits_str,')'])
title('Signal Evolution')
set(gca,'CLIM',CL)
s3 = subplot(1,3,3);
hold on;
imagesc(w./wunits,detune_ar./alpha,IiN);c = colorbar;c.Title.String = 'db';
ylabel(strcat('\Delta_',options.SweepResonance,' (\alpha_p)'))
set(gca,'YDir',SweepDir_str)
set(gca,'CLIM',CL)
xlabel(['frequnecy (',wunits_str,')'])
title('Idler Evolution')
linkaxes([s1,s2,s3],'xy')
xlim([min(w),max(w)]./wunits)
ylim([min(detune_ar),max(detune_ar)]./alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 4 Pump Signal And Idler NON Normalized temporal Evolution %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cp = max(max(Ip)); Cs = max(max(Is)); Ci = max(max(Ii)); 
for i = 1:Nrt
    IpN = Ip(i,:)./Cp;
    IsN = Is(i,:)./Cs;
    IiN = Ii(i,:)./Ci;
end
figure(fig_num+4);clf;
s1 = subplot(1,3,1);
hold on;
imagesc(t./tunits,detune_ar./alpha,IpN);colorbar;
ylabel(strcat('\Delta_',options.SweepResonance,' (\alpha_p)'))
set(gca,'YDir',SweepDir_str)
xlabel(['time (',tunits_str,')'])
title('Pump Evolution')
s2 = subplot(1,3,2);
hold on;
imagesc(t./tunits,detune_ar./alpha,IsN);colorbar;
ylabel(strcat('\Delta_',options.SweepResonance,' (\alpha_p)'))
set(gca,'YDir',SweepDir_str)
xlabel(['time (',tunits_str,')'])
title('Signal Evolution')
s3 = subplot(1,3,3);
hold on;
imagesc(t./tunits,detune_ar./alpha,IiN);colorbar;
ylabel(strcat('\Delta_',options.SweepResonance,' (\alpha_p)'))
set(gca,'YDir',SweepDir_str)
xlabel(['time (',tunits_str,')'])
title('Idler Evolution')
linkaxes([s1,s2,s3],'xy')
xlim([min(t),max(t)]./tunits)
ylim([min(detune_ar),max(detune_ar)]./alpha)
