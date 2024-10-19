%% 参数

clear 
clc
close
syms s
x=linspace(2,19,18);
t=10.^(x./10);
%% exact by numerical integration
exact1=1./pi.*exp(-t.*sin(pi./8).^2./sin(s).^2);
exact2=1./pi.*exp(-t.*sin(pi./16).^2./sin(s).^2);
intexact1=int(exact1,s,0,7.*pi./8);
intexact2=int(exact2,s,0,15.*pi./16);
%% our 2−term lower bound our 3−term lower bound
Qlbkw21=(7./24).*exp(-sin(pi./8).^2.*t.*(cot(7.*pi./24)-cot(14.*pi./24))./(7.*pi./24))...
    +(7./24).*exp(-sin(pi./8).^2.*t.*(cot(14.*pi./24)-cot(21.*pi./24))./(7.*pi./24));
Qlbkw22=(15./48).*exp(-sin(pi./16).^2.*t.*(cot(15.*pi./48)-cot(30.*pi./48))./(15.*pi./48))...
    +(15./48).*exp(-sin(pi./16).^2.*t.*(cot(30.*pi./48)-cot(45.*pi./48))./(15.*pi./48));
Qlbkw31=(7./32).*exp(-sin(pi./8).^2.*t.*(cot(7.*pi./32)-cot(14.*pi./32))./(7.*pi./32))...
    +(7./32).*exp(-sin(pi./8).^2.*t.*(cot(14.*pi./32)-cot(21.*pi./32))./(7.*pi./32))...
    +(7./32).*exp(-sin(pi./8).^2.*t.*(cot(21.*pi./32)-cot(28.*pi./32))./(7.*pi./32));
Qlbkw32=(7./32).*exp(-sin(pi./16).^2.*t.*(cot(7.*pi./32)-cot(14.*pi./32))./(7.*pi./32))...
    +(7./32).*exp(-sin(pi./16).^2.*t.*(cot(14.*pi./32)-cot(21.*pi./32))./(7.*pi./32))...
    +(7./32).*exp(-sin(pi./16).^2.*t.*(cot(21.*pi./32)-cot(28.*pi./32))./(7.*pi./32));

%% Chiani 3−term lower bound
Qchiani31=(7./32).*exp(-sin(pi./8).^2.*t./(sin(7.*pi./32)).^2)...
        +(7./32).*exp(-sin(pi./8).^2.*t/(sin(21.*pi./32)).^2)...
        +(5./32).*exp(-sin(pi./8).^2.*t./(sin(28*pi./32)).^2);

Qchiani32=(7./32).*exp(-sin(pi./16).^2.*t./(sin(7.*pi./32)).^2)...
        +(7./32).*exp(-sin(pi./16).^2.*t/(sin(21.*pi./32)).^2)...
        +(5./32).*exp(-sin(pi./16).^2.*t./(sin(28*pi./32)).^2);
% Qchiani31=(7./32).*exp(-sin(pi./8).^2.*t./(sin(7.*pi./32)).^2)...
%         +(7./32).*exp(-sin(pi./8).^2.*t/(sin(21.*pi./32)).^2)...
%         +(5./32).*exp(-sin(pi./8).^2.*t./(sin(28*pi./32)).^2);

% Qchiani32=(7./32).*exp(-sin(pi./16).^2.*t./(sin(7.*pi./32)).^2)...
%         +(7./32).*exp(-sin(pi./16).^2.*t/(sin(21.*pi./32)).^2)...
%         +(5./32).*exp(-sin(pi./16).^2.*t./(sin(28*pi./32)).^2);


%% Chiani 3−term upper bound

Qchiani_8=(14./32).*exp(-sin(pi./8).^2.*t./(sin(14.*pi./32)).^2)...
    +(7./32).*exp(-sin(pi./8).^2.*t/(sin(14.*pi./32)).^2)...
    +(7./32).*exp(-sin(pi./8).^2.*t./(sin(21*pi./32)).^2);

Qchiani_16=(30./64).*exp(-sin(pi./16).^2.*t./(sin(30.*pi./64)).^2)...
    +(30./64).*exp(-sin(pi./16).^2.*t/(sin(30.*pi./64)).^2)...
    +(15./64).*exp(-sin(pi./16).^2.*t./(sin(60*pi./64)).^2);

%% Abreu 3 or 4−term lower bound

b=2.*t.*(sin(1.*pi./8)).^2;
f1=(3+b-((b-1).^2+8).^0.5).^0.5;


f2=(4-f1.^2).^0.5;
al1=3.*f1./(4+f2);
al2=pi.*f1./(4+(pi-2).*f2);
% ab1=(al2)./pi.*exp(-sin(pi./8).^2.*t.^1./((sin(al2./2)).^2))...
%    +(al2-al1)./(2.*pi).*(exp(-sin(pi./8).^2.*t.^1./((sin(al2)).^2))+exp(-sin(pi./8).^2.*t.^1./((sin(14.*pi./32)).^2)))...  
%    +(5.*pi./32)./(2.*pi).*(exp(-sin(pi./8).^2.*t.^1./((sin(14.*pi./32)).^2))+exp(-sin(pi./8).^2.*t.^1./((sin(21.*pi./32)).^2)));
ab1=al1.*exp(-t*(sin(pi/8))^2./(sin(al1./2)).^2)./pi+...
    (1/2-al2./pi).*(exp(-t*(sin(pi/8))^2./(sin(pi-al1)).^2)+...
    exp(-t*(sin(pi/8))^2./(sin(al1)).^2))+...
    (al1./pi-1/8).*exp(-t*(sin(pi/8))^2./(sin(15/16*pi-al1./2)).^2);

%  %ab2=(al1)./pi.*exp(-sin(pi./16).^2.*t.^1./((sin(al1./2)).^2))...
%   + (al2-al1)./(2.*pi).*(exp(-sin(pi./16).^2.*t.^1./((sin(al2)).^2))+exp(-sin(pi./16).^2.*t.^1./((sin(al1./2)).^2)))...    
%   +(13.*pi./64)./(2.*pi).*(exp(-sin(pi./16).^2.*t.^1./((sin(45.*pi./64)).^2))+exp(-(sin(pi./16)).^2.*t.^1./((sin(16.*pi./32)).^2)));
b=2.*t.*(sin(1.*pi./16)).^2;
f1=(3+b-((b-1).^2+8).^0.5).^0.5;
f2=(4-f1.^2).^0.5;
al1=3.*f1./(4+f2);
al2=pi.*f1./(4+(pi-2).*f2);
ab2= al1.*exp(-t*(sin(pi/16))^2./(sin(al1./2)).^2)./pi+...
    (1/2-al2./pi).*(exp(-t*(sin(pi/16))^2./(sin(pi-al1)).^2)+...
    exp(-t*(sin(pi/16))^2./(sin(al1)).^2))+...
    (al1./pi-1/16).*exp(-t*(sin(pi/16))^2./(sin(31/32*pi-al1./2)).^2);


%% 画画
p1=semilogy(x,intexact1,'r');
    hold on
semilogy(x,intexact2,'r');
p2=semilogy(x,Qlbkw21,'k.-','MarkerSize',10);
semilogy(x,Qlbkw22,'k.-','MarkerSize',10);
p3=semilogy(x,Qlbkw31,'k');
semilogy(x,Qlbkw32,'k');
p4=semilogy(x,ab1,'m-.','linewidth',1.0);
semilogy(x,ab2,'m-.','linewidth',1.0);
p5=semilogy(x,Qchiani31,'b--','linewidth',1.0);
semilogy(x,Qchiani32,'b--','linewidth',1.0);
p6=semilogy(x,Qchiani_8,'b:','linewidth',1.6);
semilogy(x,Qchiani_16,'b:','linewidth',1.6);
xlabel('SNR(dB)')
ylabel('SEP')
axis([2 19 10^(-5) 1]);
legend([p2 p3 p5 p4 p6 p1],{'our 2-term lower bound','our 3-term lower bound',...
    'Chiani 3-term lower bound','Abreu 3 or 4−term lower bound','Chiani 3-term upper bound',...
    'exact by numerical integration'},'Location',"southwest",'fontsize',13)
grid on
