clear all;
nk = 2^10;
dk = linspace(-10,10,nk)*pi;
GDDi = -1096e-27;
nw = 2.^10;
w = (-nw/2:nw/2-1)*4*2*pi./nw*1e12;
khat = -pi/160/2*0+ 1j*GDDi.*w.^2./2;
II = zeros(nk,nw);
L = 1;
for i = 1:nk
   DK = dk(i);
   II(i,:) = (exp((1j*DK+khat).*L) - 1 - (khat+1j*DK).*L)./(khat+1j*DK).^2;
   % II(i,:) = (1-exp(-(1j*DK+khat).*L) - (khat+1j*DK).*L)./(khat+1j*DK).^2;
   P(i,:) = real(II(i,:));
   Q(i,:) = imag(II(i,:));
   PbyQ(i) = P(i,nw/2+1)./Q(i,nw/2+1);

   P(i,:) = P(i,:)./max(P(i,:));
   Q(i,:) = Q(i,:)./max(abs(Q(i,:)));
end

figure(1);clf;
subplot(2,1,1);
imagesc(w,dk,P);colorbar;
subplot(2,1,2);
imagesc(w,dk,Q);colorbar;

figure(2);clf;
plot(dk./pi,(P(:,nw/2+1)));hold on;

plot(dk./pi,(Q(:,nw/2+1)))
grid on;

figure(3);clf;plot(dk./pi,1./abs(PbyQ));set(gca,'YSCALE','log'); ylim([1e-2,1e4]);grid on;