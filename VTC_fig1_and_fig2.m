%% 参数
clear 
clc
close
x=linspace(-6,15,26);
t=10.^(x./20);
%% Qlbkw自己的上下界
Qlbkw1=exp((-2./pi).*t.^2)./4;						%自己的一项
Qlbkw2=exp((-2.*(3.^0.5)./pi).*t.^2)./6+exp((-3.^0.5./pi).*t.^2)./6;		%自己的两项
Qlbkw3=exp(-6.*(3.^0.5-1).*t.^2/pi)/12+exp(-2.*(3-3.^0.5).*t.^2./pi)/12+...
    exp(-3^0.5.*t.^2./pi)./6;       									%自己的三项

%% chniani的上下界
Qlbcds3=exp(-2.*t.^2./3)./6+exp(-t.^2)/12+exp(-2.*t.^2)./12;			%chi的
Qub=exp(-2.*t.^2)./6+exp(-t.^2)./12+exp(-t.^2./2)./4;
%%  abreu的上下界
f1=(t.^2+3-((t.^2-1).^2+8).^0.5).^0.5;
f2=(4-f1.^2).^0.5;
Qlba3=3.*f1./(4.*pi+pi.*f2).*exp(-2.*t.^2./(2-f2))+...
    (4+(pi-2).*f2-2.*f1)./(16+4.*(pi-2).*f2).*(exp(-(t.^2)./2)+exp(-2.*t.^2./f1.^2));			%ab的
%% 实际的q函数
exactQ=qfunc(t);

%% 草稿
% f1=(t.^2+3-((t.^2-1).^2+8).^0.5).^0.5;
% f2=(4-f1.^2).^0.5;
% al1=3.*f1./(4+f2);
% al2=pi.*f1./(4+(pi-2).*f2);
% Qlba4=3.*f1./(4.*pi+pi.*f2).*exp(-1.*t.^2./(2.*sin(al1./2).^2))+...
%     (4+(pi-2).*f2-2.*f1)./(16+4.*(pi-2).*f2).*(exp(-(t.^2./(2.*sin(al1).^2)))+exp(-(t.^2./(2.*sin(pi./2).^2))));


%% 画画
% imshow(uint8(null),'border','tight','initialmagnification','fit')
semilogy(x,Qlbkw1,'-k')
hold on
semilogy(x,Qlbkw2,'.-k')
semilogy(x,Qlbkw3,'k-+')
semilogy(x,Qlbcds3,'--b')
semilogy(x,Qlba3,'m-.')
semilogy(x,Qub,'o:b')
semilogy(x,exactQ,'r')
% imshow(uint8(data),'border','tight','initialmagnification','fit')
xlabel('x(dB)=20*log(10)x')
ylabel('Q(x)')
legend('Q_{LB-KW-1}_(x)','Q_{LB-KW-2}_(x)','Q_{LB-KW-3}_(x)',...
    'Q_{LB-CDS-3}_(x)','Q_{LB-A-3}_(x)','Q_{UB-CDS-3}_(x)',...
    'exactQ','Location','southwest','FontSize',16)
axis([-6 5 0.03 0.3])

figure
semilogy(x,Qlbkw1,'-k')
hold on
semilogy(x,Qlbkw2,'.-k')
semilogy(x,Qlbkw3,'k-+')
semilogy(x,Qlbcds3,'--b')
semilogy(x,Qlba3,'m-.')
semilogy(x,Qub,'o:b')
semilogy(x,exactQ,'r')
% imshow(uint8(data),'border','tight','initialmagnification','fit')
xlabel('x(dB)=20*log(10)x')
ylabel('Q(x)')
legend('Q_{LB-KW-1}_(x)','Q_{LB-KW-2}_(x)','Q_{LB-KW-3}_(x)',...
    'Q_{LB-CDS-3}_(x)','Q_{LB-A-3}_(x)','Q_{UB-CDS-3}_(x)',...
    'exactQ','Location','southwest','FontSize',16)
axis([8 15 1*10^(-10) 5*10^(-3)])



