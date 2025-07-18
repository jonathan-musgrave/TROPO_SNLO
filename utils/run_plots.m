%% PLOTTING AND PROCCESSING
frep = 400E6; % Freespace frep
n = 2.210;
tact = L./(c/n);
CC = 1;
CCI =  1;
CCP = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IP = abs(LP).^2*CCP;
IS = abs(LS).^2*CC;
II = abs(LI).^2*CCI;
PZp = abs(AZp).^2*CCP;
PZs = abs(AZs).^2*CC;
PZi = abs(AZi).^2*CCI;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%------stability of peak power-----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Intensity_P = IP;
Intensity_S = IS;
Intensity_I = II;
AA=zeros(Nrt,1);
for ind = 1:Nrt
    PP_P(ind)= max(Intensity_P(ind,:));
    PP_S(ind)= max(Intensity_S(ind,:));
    PP_I(ind)= max(Intensity_I(ind,:));
end
roundtrip=1:1:Nrt;
FS=10;
FS2=10;
LW=2;
figure(fig_num+1);clf;hold on;
plot(roundtrip,PP_P,'-b','linewidth',LW)
plot(roundtrip,PP_S,'-g','linewidth',LW)
plot(roundtrip,PP_I,'-r','linewidth',LW)
xlabel('round-trip number','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
ylabel('peak power (W)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
legend('Pump','Signal','Idler')
%set(gca,'xtick',5000:5000:20000);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind = Nrt;
tp = 144E-15;
sechSS = max(IS(ind,:))*abs(sech(t/tp)).^2;
XL = [-6,6]./10;
phiYL = [-5,5];
inst_wp = -gradient(unwrap(angle(LP(ind,:))),dt)/2/pi/1E12;
inst_ws = -gradient(unwrap(angle(LS(ind,:))),dt)/2/pi/1E12;
inst_wi = -gradient(unwrap(angle(LI(ind,:))),dt)/2/pi/1E12;
% figure(fig_num+2);clf;
% subplot(2,2,1);
% yyaxis('left')
% % plot(t*1e12,IP(ind,:),'-','Color',[0,0,0.7],'linewidth',LW);hold on;
% plot(t*1e12,IS(ind,:),'-','Color',[0,0.7,0],'linewidth',LW);hold on;
% plot(t*1e12,II(ind,:),'-','Color',[0.7,0,0],'linewidth',LW);hold on;
% xlim([-15,15]./15)
% set(gca,'YCOLOR','k')
% ylim([-.05,max(IS(ind,:)).*1.2])
% ylabel('Power (W)')
% yyaxis('right')
% grid on;
% % plot(t*1e12,inst_wp,'--','Color',[0,0,0.7],'linewidth',LW);hold on;
% plot(t*1e12,inst_ws,'--','Color',[0,0.7,0],'linewidth',LW);hold on;
% plot(t*1e12,inst_wi,'--','Color',[0.7,0,0],'linewidth',LW);hold on;
% xlim([-15,15]./10)
% ylim([-0.5,0.5]*3/2)
% set(gca,'YCOLOR','k')
% ylabel('Chirp (THz)')
% xlabel('Time (ps)')
% 
% 
%%%%%%%
SPP = abs(fftshift(ifft(ifftshift(LP(ind,:))))).^2;
SPP = SPP./max(SPP);
SPS = abs(fftshift(ifft(ifftshift(LS(ind,:))))).^2;
SPS = SPS./max(SPS);
SPI = abs(fftshift(ifft(ifftshift(LI(ind,:))))).^2;
SPI = SPI./max(SPI);

% subplot(2,2,2);
% BV = -50;
% LW = 1;
% % plot(w/2/pi/1e12,10*log10(SPP),'Color',[0,0,0.7],'linewidth',LW)
% % hold on
% stem(w/2/pi/1e12,10*log10(SPS)+20,'Color',[0,0.7,0],'linewidth',LW,'basevalue',BV,'marker','none')
% hold on
% stem(w/2/pi/1e12,10*log10(SPI),'Color',[0.7,0,0],'linewidth',LW,'basevalue',BV,'marker','none')
% ylabel('Power (dB)')
% ylim([BV,30])
% xlim([-1,1]*9)
% xlabel('Frequency (THz)')
% grid on;
% 
% 
% 
% subplot(2,2,3);
% yyaxis('left')
% plot(t*1e12,IP(ind,:),'-','Color',[0,0,0.7],'linewidth',LW);hold on;
% ylim([0,0.2])
% set(gca,'YCOLOR','k')
% ylabel('Power (W)')
% yyaxis('right')
% plot(t*1e12,inst_wp,'--','Color',[0,0,0.7],'linewidth',LW);hold on;
% xlim([-15,15]./15)
% ylim([-1,1])
% set(gca,'YCOLOR','k')
% ylabel('Chirp (THz)')
% xlabel('Time (ps)')
% grid on;
% 
% subplot(2,2,4);
% BV = -50;
% LW = 1;
% % plot(w/2/pi/1e12,10*log10(SPP),'Color',[0,0,0.7],'linewidth',LW)
% % hold on
% stem(w/2/pi/1e12,10*log10(SPP)+25,'Color',[0,0,0.7],'linewidth',LW,'basevalue',BV,'marker','none')
% hold on
% % stem(w/2/pi/1e12,10*log10(SPI),'Color',[0.7,0,0],'linewidth',LW,'basevalue',BV,'marker','none')
% ylabel('Power (dB)')
% ylim([BV,30])
% xlim([-1,1]*9)
% xlabel('Frequency (THz)')
% grid on;





try
    XL = 5*fwhm(t*1e12,IP(ind,:));
catch
   XL = max(t*1e12); 
end
figure(fig_num+3);clf;
subplot(1,1,1);
yyaxis left
plot(t*1e12, IP(ind,:),'Color',[0,0,0.7],'linewidth',LW)
ylabel('peak power (W)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
hold on
set(gca,'ylim',[0 max(IP(ind,:))]*1.2);
subplot(1,1,1);
yyaxis('right')
plot(t*1e12, II(ind,:),'r-','linewidth',LW)
hold on
set(gca,'ylim',[0 max(IS(ind,:))]*1.2);
ylabel('peak power (W)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'Ycolor','k')
subplot(1,1,1);
yyaxis('right');
if lam0.s == lam0.i
plot(t*1e12, IS(ind,:),'--','Color',[0,0.7,0],'linewidth',LW)
else
plot(t*1e12, IS(ind,:),'-','Color',[0,0.7,0],'linewidth',LW)
end
hold on
% plot(t*1e12, sechSS,'g--','linewidth',LW)
hold off
ylabel('peak power (W) p','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'Ycolor',[0,0,0.7])
set(gca,'ylim',[0 max([II(ind,:),IS(ind,:)])*1.2]);
[~,x0inx] = max(IP(ind,:));
xlim([-1,1].*max(t)*1e12)
ylabel('peak power (W) [s,i] ','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'YColor','k')
%
inst_wp = -gradient(unwrap(angle(LP(ind,:))),dt)/2/pi/1E12;
inst_ws = -gradient(unwrap(angle(LS(ind,:))),dt)/2/pi/1E12;
inst_wi = -gradient(unwrap(angle(LI(ind,:))),dt)/2/pi/1E12;
%}
% %{
% inst_wp = unwrap(angle(LP(ind,:)))/pi;
% inst_ws = unwrap(angle(LS(ind,:)))/pi;
% inst_wi = unwrap(angle(LI(ind,:)))/pi;
% %}
% subplot(3,1,1);
% yyaxis right
% plot(t*1e12, inst_wp,'k--','linewidth',LW)
% hold on
% set(gca,'xlim',XL);
% set(gca,'ylim',phiYL);
% subplot(3,1,2);
% yyaxis right
% plot(t*1e12, inst_ws,'r--','linewidth',LW)
% hold on
% set(gca,'xlim',XL);
% set(gca,'ylim',phiYL);
% subplot(3,1,3);
% yyaxis right
% plot(t*1e12, inst_wi,'b--','linewidth',LW)
% hold off
% set(gca,'Ycolor','k')
% ylabel('chirp (THz)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
% xlabel('fast time (ps)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
% set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
% set(gca,'xlim',XL);
% set(gca,'ylim',phiYL);
% hco=legend('pump','signal','idler','sech^{2} fitting','location','northeast');legend boxoff;
% set(hco,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %
% [fval,fpos]=find(IS(Nrt,:)>=max(IS(Nrt,:))/2);
% pulse_duration1 = length(fpos)*dt*1e15
% [fval,fpos]=find(II(Nrt,:)>=max(II(Nrt,:))/2);
% pulse_duration2 = length(fpos)*dt*1e15
% bandwidth_s = 0.315/pulse_duration1*1E3
% bandwidth_i = 0.315/pulse_duration2*1E3
% max_BW_s = detunes/2/pi/(beta_offset3*L)/1E12
% max_BW_i = detunei/2/pi/(beta_offset3*L)/1E12
% %}

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inst_w = -gradient(unwrap(angle(LS(ind,:))),dt);
% figure(fig_num+4);clf;
% yyaxis left
% plot(t*1e12,IS(ind,:),'linewidth',LW);
% xlabel('time (ps)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
% ylabel('peak power (W)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
% set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
% %set(gca,'xlim',time_range)
% %set(gca,'ycolor','r');
% yyaxis right
% plot(t*1e12,inst_w/2/pi/1e12,'linewidth',LW);
% ylabel('instaneous frequency (THz)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
% set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----spectrum in linear coordinates---%%%%%%%%%%%%%%%%%%%%%%%%%%
% ind=Nrt;
% 
% SPP = abs(fftshift(ifft(ifftshift(LP(ind,:))))).^2;
% SPP = SPP./max(SPP);
% SPS = abs(fftshift(ifft(ifftshift(LS(ind,:))))).^2;
% SPS = SPS./max(SPS);
% SPI = abs(fftshift(ifft(ifftshift(LI(ind,:))))).^2;
% SPI = SPI./max(SPI);

%{
figure(100)
TLP = abs(fftshift(fft(ifftshift(sqrt(SPS))))).^2./max(abs(fftshift(fft(ifftshift(sqrt(SPS))))).^2);
plot(t*1E12,TLP,'r')
hold on
plot(t*1E12,IS(Nrt,:)./max(IS(Nrt,:)),'g')
hold off
[fval,fpos]=find(TLP>=max(TLP)/2);
pulse_duration3 = length(fpos)*dt*1e15

[fval,fpos]=find(SPS>=-10);
BW1 = length(fpos)*dw/2/pi/1E12
[fval,fpos]=find(SPI>=-10);
BW2 = length(fpos)*dw/2/pi/1E12
%}
%
% [fval,fpos]=max(SPS);
% fres1 = ((Nw+1)/2-fpos)*dw/2/pi/1E12
% [fval,fpos]=max(SPI);
% fres2 = ((Nw+1)/2-fpos)*dw/2/pi/1E12
% 
% beta2_eff = (beta2_s + beta2_i);
% abs(beta_offset3/beta2_eff)/2/pi/1E12;
% 
% 
% tp = 128E-15;
% sechSP = abs(fftshift(ifft(ifftshift(sech(t/tp))))).^2;
% figure(fig_num+5);clf;
% plot(w/2/pi/1e12,10*log10(SPP),'k','linewidth',LW)
% hold on
% plot(w/2/pi/1e12,10*log10(SPS),'r','linewidth',LW)
% hold on
% plot(w/2/pi/1e12,10*log10(SPI),'b','linewidth',LW)
% hold on
% plot(w/2/pi/1e12,10*log10(sechSP./max(sechSP)),'g--','linewidth',LW)
% hold off
% ylabel('power (dB)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
% xlabel('frequency (THz)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
% set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
% set(gca,'xlim',[-5 5]);
% set(gca,'ylim',[-50 0]);
% hco=legend('pump','signal','idler','sech^{2} fitting','location','south');legend boxoff;
% set(hco,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

tp = 128E-15;
sechSP = abs(fftshift(ifft(ifftshift(sech(t/tp))))).^2;
figure(fig_num+6);clf;
BV = -50;
%stem(w/2/pi/1e12,10*log10(SPP),'k','Marker', 'none','basevalue',BV)
%hold on
subplot(1,2,1);
stem(w/2/pi/1e12,20+10*log10(SPI),'Color',[0.7,0,0],'Marker', 'none','basevalue',BV)
hold on
stem(w/2/pi/1e12,10*log10(SPS),'Color',[0,0.7,0],'Marker', 'none','basevalue',BV)
hold off
ylabel('power (dB)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
xlabel('frequency - {\it\omega_{q=p,s,i}}/(2\pi) (THz)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
set(gca,'xlim',[min(w),max(w)]./2/pi/1e12);
set(gca,'ylim',[-40 30]);
%hco=legend('idler','signal','sech^{2}','location','northeast');legend boxoff;
hco=legend('idler','signal','location','northeast');legend boxoff;
set(hco,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
subplot(1,2,2)
stem(w/2/pi/1e12,20+10*log10(SPP),'Color',[0,0,0.7],'Marker', 'none','basevalue',BV)
%hold on
%plot(w/2/pi/1e12,10*log10(sechSP./max(sechSP)),'g--','linewidth',LW)
hold off
ylabel('power (dB)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
xlabel('frequency - {\it\omega_{q=p,s,i}}/(2\pi) (THz)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
set(gca,'xlim',[min(w),max(w)]./2/pi/1e12);
set(gca,'ylim',[-50 20]);
%hco=legend('idler','signal','sech^{2}','location','northeast');legend boxoff;
hco=legend('pump');legend boxoff;
set(hco,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%{
x = (dk - beta_offset1.*w - beta2_p/2.*w.^2)*L;
%ko = -alpha2/L/2 + 1j*(GVM*w+GVD2/2*w.^2);
%x = (dk + 1j*ko)*L;
I = (1-exp(-1j*x)-1j*x)./x.^2;
Ir = fftshift(fft(ifftshift(real(I))));
Ii = -fftshift(fft(ifftshift(imag(I))));
Irr = LS(ind,:).*LI(ind,:);
Iii = LS(ind,:).*LI(ind,:);
I1 = fftshift(ifft(ifftshift(LS(ind,:))));
I2 = fftshift(ifft(ifftshift(LI(ind,:))));
I3 = fftshift(ifft(ifftshift(LS(ind,:).*LI(ind,:))));
I33 = fftshift(ifft(ifftshift(LS(ind,:).*LI(ind,:)))).*I;
I4 = abs(Irr);
I5 = abs(Iii);


figure(666)
plot(t*1e12,Irr,'r-','linewidth',LW)
hold on
plot(t*1e12,Iii,'g-','linewidth',LW)
hold off


figure(777)
yyaxis left
plot(w/2/pi/1e12,abs(I1).^2,'r-','linewidth',LW)
hold on
plot(w/2/pi/1e12,abs(I2).^2,'g-','linewidth',LW)
hold on
plot(w/2/pi/1e12,abs(I3).^2,'b-','linewidth',LW)
hold off
ylabel('power','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
xlabel('frequency shift (THz)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
set(gca,'xlim',[-10 10]);
amp = [max(abs(I1)) max(abs(I2)) max(abs(I3))]*1E3

yyaxis right
plot(w/2/pi/1e12,real(I),'m','linewidth',LW)
hold on
plot(w/2/pi/1e12,-imag(I),'k','linewidth',LW)
hold off
hco=legend('FFT(S)','FFT(I)','FFT(S*I)','real(I(\Omega))','imag(I(\Omega))','location','southwest');legend boxoff;
set(hco,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)


figure(888)
plot(t*1E12,I4,'b','linewidth',LW)
hold on
plot(t*1E12,abs(LS(Nrt,:)),'r','linewidth',LW)
hold on
plot(t*1E12,abs(LI(Nrt,:)),'g','linewidth',LW)
hold off
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)


figure(999)
plot(t*1E12,(real(LS(ind,:))),'r-','linewidth',LW)
hold on
plot(t*1E12,(imag(LS(ind,:))),'r--','linewidth',LW)
hold on
plot(t*1E12,(real(LI(ind,:))),'g-','linewidth',LW)
hold on
plot(t*1E12,(imag(LI(ind,:))),'g--','linewidth',LW)
hold off

plot(real(LS(ind,:)),imag(LS(ind,:)),'r','linewidth',LW)
hold on
plot(real(LI(ind,:)),imag(LI(ind,:)),'g','linewidth',LW)
hold on
plot(real(LS(ind,:).*LI(ind,:)),imag(LS(ind,:).*LI(ind,:)),'b','linewidth',LW)
hold off
set(gca,'xlim',[-0.5 0.5]);
set(gca,'ylim',[-0.5 0.5]);
axis square


figure(111)
plot(real(I1),imag(I1),'r','linewidth',LW)
hold on
plot(real(I2),imag(I2),'g','linewidth',LW)
hold on
plot(real(I3),imag(I3),'b','linewidth',LW)
hold off
axis square
%}

%{
AAA = fftshift(ifft(ifftshift(LS(ind,:)))).*fftshift(ifft(ifftshift(LI(ind,:))));
convA = conv(AAA,Ir./max(abs(Ir)),'same');
figure(777)
plot(t*1E12,abs(convA)./max(abs(convA)))
%}









trange = tww/2-dt;
%%%%%%%%%%%%%%%%%%%%%%%%%-----pulse evolutionfor each roundtrip-----%%%%%%%%%%%%%%%%%%%%%%%%%

roundtrip = 1:1:Nrt;
tickrange = floor((Nw+1)/2)-floor(trange/dt):1:floor((Nw+1)/2)+floor(trange/dt);
tt = t(tickrange)*1E12;
[TM,RTM] = meshgrid(tt,roundtrip/1000);
PS = zeros(Nrt,Nw);
PI = zeros(Nrt,Nw);
PP = zeros(Nrt,Nw);
for i = 1:1:Nrt
    PS(i,:) = IS(i,:)./max(IS(i,:));
    PI(i,:) = II(i,:)./max(II(i,:));
    PP(i,:) = IP(i,:)./max(IP(i,:));
end
figure(fig_num+7);clf;
ax1 = subplot(1,3,1);
hh=pcolor(TM,RTM,PP(roundtrip,tickrange));
set(hh,'edgecolor','k','Marker','*')
shading interp
set(gca,'LineWidth',LW,'FontSize',FS)
% load redblue;
colormap(ax1,'hot');
colorbar
ylabel('roundtrip','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
xlabel('time (ps)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
%set(gca,'ytick',100:100:500);
hco=colorbar;
set(hco,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold');
title('Pump Evolution')
set(gca,'CLIM',[0,1])
%set(get(hco,'Title'),'string','signal W');
ax2 = subplot(1,3,2);
hh=pcolor(TM,RTM,PS(roundtrip,tickrange));
set(hh,'edgecolor','k','Marker','*')
shading interp
set(gca,'LineWidth',LW,'FontSize',FS)
% load redblue;
colormap(ax2,'parula');
colorbar
ylabel('roundtrip','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
xlabel('time (ps)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
%set(gca,'ytick',100:100:500);
hco=colorbar;
set(hco,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold');
title('Signal Evolution')
set(gca,'CLIM',[0,1])
%set(get(hco,'Title'),'string','signal W');

ax3 = subplot(1,3,3);
hh=pcolor(TM,RTM,PI(roundtrip,tickrange));
set(hh,'edgecolor','k','Marker','*')
shading interp
set(gca,'LineWidth',LW,'FontSize',FS)
% load redblue;
colormap(ax3,'redblue');
colorbar
ylabel('roundtrip','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
xlabel('time (ps)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
%set(gca,'ytick',100:100:500);
hco=colorbar;
set(hco,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold');
title('Idler Evolution')
%set(get(hco,'Title'),'string','signal W');


% 
% %%%%%%%%%%%%%%%%%%%%%%%%%-----pulse evolutionfor each roundtrip-----%%%%%%%%%%%%%%%%%%%%%%%%%
% roundtrip = 1:1:Nrt;
% tickrange = floor((Nw+1)/2)-floor(trange/dt):1:floor((Nw+1)/2)+floor(trange/dt);
% tt = t(tickrange)*1E12;
% [TM,RTM] = meshgrid(tt,roundtrip/1000);
% figure(66);clf;
% hh=pcolor(TM,RTM,II(roundtrip,tickrange));
% set(hh,'edgecolor','k','Marker','*')
% shading interp
% set(gca,'LineWidth',LW,'FontSize',FS)
% % load seismic;
% colormap(hot);
% colorbar
% ylabel('roundtrip','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
% xlabel('time (ps)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
% set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
% %set(gca,'ytick',100:100:500);
% hco=colorbar;
% set(hco,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold');
% title('Idler Evolution')
% %set(get(hco,'Title'),'string','signal W');



% 
% roundtrip = 1:1:Nrt;
% frange = 10E12;
% wtickrange = floor((Nw+1)/2)-floor(frange/(dw/2/pi)):1:floor((Nw+1)/2)+floor(frange/(dw/2/pi));
% ww = w(wtickrange)/2/pi/1E12;
% [WM,RTM] = meshgrid(ww,roundtrip);
% SPEC = zeros(Nrt,Nw);
% for ind = 1:Nrt
%     SPEC(ind,:) = abs(fftshift(ifft(ifftshift(LS(ind,:))))).^2./max(abs(fftshift(ifft(ifftshift(LS(ind,:))))).^2);
% end
% figure(7)
% hh=pcolor(WM,RTM,10*log10(SPEC(roundtrip,wtickrange)));
% set(hh,'edgecolor','k','Marker','*')
% shading interp
% set(gca,'LineWidth',LW,'FontSize',FS)
% ylabel('roundtrip','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
% xlabel('rrequency (THz)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
% set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
% %set(gca,'ytick',100:100:500);
% hco=colorbar;
% set(hco,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold');
% %set(get(hco,'Title'),'string','Power (dBm)');
% 
% 
% 
% trange = 0.4E-12;
% %%%%%%%%%%%%%%%%%%%%%%%%%-----pulse evolutionfor each roundtrip-----%%%%%%%%%%%%%%%%%%%%%%%%%
% tickrange = floor((Nw+1)/2)-floor(trange/dt)-1:1:floor((Nw+1)/2)+floor(trange/dt)+1;
% tt = t(tickrange)*1E12;
% zz = z*1E3;
% [TM,RZ] = meshgrid(tt,zz);
% figure(8)
% hh=pcolor(TM,RZ,PZs(:,tickrange));
% set(hh,'edgecolor','k','Marker','*')
% shading interp
% set(gca,'LineWidth',LW,'FontSize',FS)
% % load seismic;
% colormap(jet);
% colorbar
% ylabel('distance (mm)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
% xlabel('time (ps)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
% set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
% %set(gca,'ytick',1000:1000:5000);
% hco=colorbar;
% set(hco,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold');
% %set(get(hco,'Title'),'string','signal W');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%-----pulse evolutionfor each roundtrip-----%%%%%%%%%%%%%%%%%%%%%%%%%
% tickrange = floor((Nw+1)/2)-floor(trange/dt)-1:1:floor((Nw+1)/2)+floor(trange/dt)+1;
% tt = t(tickrange)*1E12;
% [TM,RZ] = meshgrid(tt,zz);
% figure(9)
% hh=pcolor(TM,RZ,PZi(:,tickrange));
% set(hh,'edgecolor','k','Marker','*')
% shading interp
% set(gca,'LineWidth',LW,'FontSize',FS)
% % load seismic;
% colormap(hsv);
% colorbar
% ylabel('distance (mm)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
% xlabel('time (ps)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
% set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
% %set(gca,'ytick',1000:1000:5000);
% hco=colorbar;
% set(hco,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold');
% %set(get(hco,'Title'),'string','signal W');