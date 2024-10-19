%% 参数
clear 
clc
close
syms s
x=linspace(6,20,14);
t=10.^(x./10);
%% 实际-exact by numerical integrations
 exact1=1./pi.*exp(-t.*sin(pi./8).^2./(1+cos(pi./8).*cos(s)));
exact2=1./pi.*exp(-t.*sin(pi./16).^2./(1+cos(pi./16).*cos(s)));
intexact1=int(exact1,s,0,7.*pi./8);
intexact2=int(exact2,s,0,15.*pi./16);

%% Chiani 2−term lower bound
Qchiani21=(7./24).*exp(-t.*sin(pi./8).^2./(1+cos(pi./8).*cos(14.*pi./24)))...
        +(7./24).*exp(-t.*sin(pi./8).^2./(1+cos(pi./8).*cos(7.*pi./24)));
  
Qchiani22=5/16.*exp(-t.*sin(pi./16).^2./(1+cos(pi./16).*cos(5.*pi./16)))...
        +5/16.*exp(-t.*sin(pi./16).^2./(1+cos(pi./16).*cos(10.*pi./16)));

%% Chiani 2−term upper bound
Qchianiup21=(12./24).*exp(-t.*sin(pi./8).^2./(1+cos(pi./8).*cos(0.*pi./24)))...
           +(9./24).*exp(-t.*sin(pi./8).^2./(1+cos(pi./8).*cos(12.*pi./24)));
Qchianiup22=(24./48).*exp(-t.*sin(pi./16).^2./(1+cos(pi./16).*cos(0.*pi./24)))...
           +(21./48).*exp(-t.*sin(pi./16).^2./(1+cos(pi./16).*cos(12.*pi./24)));
       
 %% our 2−term lower bound
 ak=2.*sqrt(2.*(2+sqrt(2)));
 t8=tan(pi./16);
%  bk1=int(ak,s,7.*pi./24,12.*pi./24);
%  bk2=int(ak,s,12.*pi./24,21.*pi./24);
%    Qlbkw21=(16./32).*exp(t.*ak.*sin(pi./8).^2.*(atan(t8.*tan(8.*pi./64))-atan(t8.*tan(16.*pi./64)))./(16.*pi./32))...
%           +(12./32).*exp(t.*ak.*sin(pi./8).^2.*(atan(t8.*tan(16.*pi./64))-atan(t8.*tan(28.*pi./64)))./(12.*pi./32));
%       
      Qlbkw21=(14./32).*exp(t.*ak.*sin(pi./8).^2.*(atan(t8.*tan(14.*pi./64))-atan(t8.*tan(21.*pi./64)))./(14.*pi./32))...
        +(12./32).*exp(t.*ak.*sin(pi./8).^2.*(atan(t8.*tan(21.*pi./64))-atan(t8.*tan(28.*pi./64)))./(12.*pi./32));
 %%  草稿
%  Qlbkw22=(7./24).*exp(t.*ak.*sin(pi./8).^2.*(atan(t8.*tan(7.*pi./48))-atan(t8.*tan(14.*pi./48)))./(12.*pi./24))...
%           +(7./24).*exp(t.*ak.*sin(pi./8).^2.*(atan(t8.*tan(14.*pi./48))-atan(t8.*tan(21.*pi./48)))./(7.*pi./24));


%% 论文里面的系数

k=2.*sqrt(2.*(2+sqrt(2)));
m1=tan(pi./32);
s1=sin(pi./16);
Qlbkw22=(36./64).*exp(t.*k.*s1.^2.*(atan(m1.*tan(30.*pi./128))-atan(m1.*tan(45.*pi./128)))./(16.*pi./64))...
    +(18./64).*exp(t.*k.*s1.^2.*(atan(m1.*tan(45.*pi./128))-atan(m1.*tan(60.*pi./128)))./(18.*pi./64));
%     Qlbkw22=(7./24).*exp(t.*sin(pi./8).*2.*(atan(tan(7.*pi./48).*btk)-atan(tan(14.*pi./48).*btk))./(14.*pi./24))...
%           +(14./24).*exp(t.*sin(pi./8).*2.*(atan(tan(14.*pi./48).*btk)-atan(tan(21.*pi./48).*btk))./(14.*pi./24));
% Qlbkw21=(30./64).*exp(t.*sin(pi./16).*2.*(atan(tan(30.*pi./128).*btk)-atan(tan(45.*pi./128).*btk))./(30.*pi./64))...
%           +(28./64).*exp(t.*sin(pi./16).*2.*(atan(tan(45.*pi./128).*btk)-atan(tan(60.*pi./128).*btk))./(28.*pi./64));


     %% 画画
p5=semilogy(x,intexact1,'r');
    hold on
    semilogy(x,intexact2,'r');
p2=semilogy(x,Qchiani21,'b--');
    semilogy(x,Qchiani22,'b--');
p4=semilogy(x,Qchianiup21,'b:','linewidth',1.6);
    semilogy(x,Qchianiup22,'b:','linewidth',1.6);
p1=semilogy(x,Qlbkw21,'k.-','MarkerSize',10);
  semilogy(x,Qlbkw22,'k.-','MarkerSize',10);
   legend([p1 p2 p4 p5],{'our 2-term lower bound','Chiani 2-term lower bound',...
    'Chiani 2-term upper bound','exact by numerical integration'},...
    'Location',"southwest",'fontsize',13);
xlabel('SNR(dB)')
ylabel('SEP')
axis([6 20 10^(-3) 1]);