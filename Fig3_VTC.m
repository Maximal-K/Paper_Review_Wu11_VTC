clear 
clc
close

syms y m a
x=linspace(10,40,20);
t=10.^(x./10);
%% 2−term Chiani lower bound
pl1=(1+2.*t./1).^(-1)./6+(1+2.*t./3).^(-1)./6;
pl2=(1+2.*t./2.5).^(-2.5)./6+(1+2.*t./7.5).^(-2.5)./6;
%pu2=(1+t./5).^(-2.5)./6+(1+2.*t./7.5).^(-2.5)./6;
%% exact by numerical integration
F1=1./sqrt(2.*pi).*exp(-a.^2./2);
f=int(F1,sqrt(m),100);
f1=@(m)f./(t.*gamma(1)).*exp(-m./t);
f2=@(m)f.*2.5.^2.5.*m.^1.5./(t.^2.5.*gamma(2.5)).*exp(-2.5.*m./t);     
j1=int(f1(m),0,10);
j2=int(f2(m),0,10);


%% our 2−term lower bound
four=(exp((-2.*(3.^0.5)./pi).*y)./6+...
    exp((-3.^0.5./pi).*y)./6)./t.*exp(-y./t);
four2=(exp((-2.*(3.^0.5)./pi).*y)./6+...
    exp((-3.^0.5./pi).*y)./6).*2.5.^2.5.*(y.^1.5)./(t.^2.5.*gamma(2.5)).*exp(-2.5.*y./t);
Iour=int(four,y,0,10000);    
Iour2=int(four2,y,0,10000);
%% 2−term Chiani upper bound
pu1=(1+t).^(-1)./4+(1+t./2).^(-1)./4;
pu2=(1+t./2.5).^(-2.5)./4+(1+t./5).^(-2.5)./4;
%% 画图
figure
p1=semilogy(x,pu1,'b-.');
hold on
p2=semilogy(x,j1,'r');
p3=semilogy(x,Iour,'k.-','MarkerSize',13);
p4=semilogy(x,pl1,'--b');
%semilogy(x,I2,'b')
semilogy(x,j2,'r');
semilogy(x,Iour2,'k.-','MarkerSize',13);
semilogy(x,pl2,'--b');
semilogy(x,pu2,'b-.')
legend([p3 p4 p1 p2],{'our 2-term lower bound','2-term Chiani lower bound',...
    '2-term Chiani upper bound','exact by numerical integration'},...
    'Location','southwest','FontSize',16)
xlabel('SNR(dB)')
ylabel('l')
set(gca,'FontSize',16)
%text(18,2*10^(-4),'o','color','k','FontSize')
axis([10 40 10^(-6) 0.04]);
