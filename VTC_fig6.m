% % 参数设置  
clc
clear
syms al            % al设置为角度
x=linspace(15,40,25);
t=10.^(x./10);
%% exact by numerical integration  通式
s=@(m)(-sin(pi./m).^2./sin(al).^2);         %关于m的表达式 
MGF=@(k,m)(1+k)./((1+k+(-s(m)).*t).*pi).*...
    exp(-k.*(-s(m)).*t./(1+k+(-s(m)).*t));  %MGF关于m，k的表达式
Pes=@(k,m)int(MGF(k,m),al,0,pi-pi./m);      %最终的通项表达式关于k，m
%% our 2−term lower bound 
%       k=2,m=8的情况
bk21_28=sin(pi./8).^2.*(cot(7.*pi./32)-cot(14.*pi./32))./(7.*pi./32);
bk22_28=sin(pi./8).^2.*(cot(14.*pi./32)-cot(21.*pi./32))./(7.*pi./32);
MGF_OUR2_28=@(k)7./32.*(1+k)./(1+k+bk21_28.*t).*...
    exp(-k.*(bk21_28).*t./(1+k+bk21_28.*t))+7./32.*(1+k)./(1+k+bk22_28.*t).*...
    exp(-k.*(bk22_28).*t./(1+k+bk22_28.*t));
%       k=10,m=8的情况
bk21_108=sin(pi./8).^2.*(cot(7.*pi./24)-cot(14.*pi./24))./(7.*pi./24);
bk22_108=sin(pi./8).^2.*(cot(14.*pi./24)-cot(21.*pi./24))./(7.*pi./24);
MGF_OUR2_108=@(k) 5./24.*(1+k)./(1+k+bk21_108.*t).*...
    exp(-k.*(bk21_108).*t./(1+k+bk21_108.*t))+10./24.*(1+k)./(1+k+bk22_108.*t).*...
    exp(-k.*(bk22_108).*t./(1+k+bk22_108.*t));
%       k=2,m=16的情况
bk21_16=sin(pi./16).^2.*(cot(15.*pi./48)-cot(30.*pi./48))./(15.*pi./48);
bk22_16=sin(pi./16).^2.*(cot(30.*pi./48)-cot(45.*pi./48))./(15.*pi./48);
MGF_OUR2_216=@(k) 12./48.*(1+k)./(1+k+bk21_16.*t).*...
    exp(-k.*(bk21_16).*t./(1+k+bk21_16.*t))+15./48.*(1+k)./(1+k+bk22_16.*t).*...
    exp(-k.*(bk22_16).*t./(1+k+bk22_16.*t));
%       k=10,m=16的情况
MGF_OUR2_1016=@(k) 12./48.*(1+k)./(1+k+bk21_16.*t).*...
    exp(-k.*(bk21_16).*t./(1+k+bk21_16.*t))+15./48.*(1+k)./(1+k+bk22_16.*t).*...
    exp(-k.*(bk22_16).*t./(1+k+bk22_16.*t));


%% our 3−term lower bound 通式
alk3=@(m) (pi-pi./m)./4;
ak3=@(m) alk3(m)./pi;
bk31=@(m)sin(pi./m).^2.*(cot(alk3(m))-cot(alk3(m).*2))./((pi-pi./m)./4);
bk32=@(m)sin(pi./m).^2.*(cot(alk3(m).*2)-cot(alk3(m).*3))./((pi-pi./m)./4);
bk33=@(m)sin(pi./m).^2.*(cot(alk3(m).*3)-cot(alk3(m).*4))./((pi-pi./m)./4);
MGF_OUR3=@(k,m) ak3(m).*(1+k)./(1+k+bk31(m).*t).*...
    exp(-k.*(bk31(m)).*t./(1+k+bk31(m).*t))+ak3(m).*(1+k)./(1+k+bk32(m).*t).*...
    exp(-k.*(bk32(m)).*t./(1+k+bk32(m).*t))+ak3(m).*(1+k)./(1+k+bk33(m).*t).*...
    exp(-k.*(bk33(m)).*t./(1+k+bk33(m).*t));

%% Chiani 3−term upper bound 通式
alk3_CU=@(m) (pi-pi./m)./4;
ak3_CU=@(m) alk3_CU(m)./pi;
bk31_CU=@(m)sin(pi./m).^2./(sin(2.*alk3_CU(m))).^2;
bk32_CU=@(m)sin(pi./m).^2./(sin(2.*alk3_CU(m))).^2;
bk33_CU=@(m)sin(pi./m).^2./(sin(3.*alk3_CU(m))).^2;
MGF_CHIANI_UB3=@(k,m) 2.*(ak3_CU(m)).*(1+k)./(1+k+bk31_CU(m).*t).*...
    exp(-k.*(bk31_CU(m)).*t./(1+k+bk31_CU(m).*t))+(ak3_CU(m)).*(1+k)./(1+k+bk32_CU(m).*t).*...
    exp(-k.*(bk32_CU(m)).*t./(1+k+bk32_CU(m).*t))+(ak3_CU(m)).*(1+k)./(1+k+bk33_CU(m).*t).*...
    exp(-k.*(bk33_CU(m)).*t./(1+k+bk33_CU(m).*t)); 

%% Chiani 3−term lower bound 有点小问题，基本无误，但是有两条线不是很齐
%m=8的系数
bk31_CL_8=(sin(pi./8)).^2./(sin(1.*7.*pi./32)).^2;
bk32_CL_8=(sin(pi./8)).^2./(sin(3.*7.*pi./32)).^2;
bk33_CL_8=(sin(pi./8)).^2./(sin(4.*7.*pi./32)).^2;
%    m=8，k=2
MGF_CHIANI_LB28=@(k) (7./32).*(1+k)./(1+k+bk31_CL_8.*t).*...
    exp(-k.*(bk31_CL_8).*t./(1+k+bk31_CL_8.*t))+(6./32).*(1+k)./(1+k+bk32_CL_8.*t).*...
    exp(-k.*(bk32_CL_8).*t./(1+k+bk32_CL_8.*t))+(7./32).*(1+k)./(1+k+bk33_CL_8.*t).*...
    exp(-k.*(bk33_CL_8).*t./(1+k+bk33_CL_8.*t));
%     m=8，k=10
MGF_CHIANI_LB108=@(k) (2./32).*(1+k)./(1+k+bk31_CL_8.*t).*...
    exp(-k.*(bk31_CL_8).*t./(1+k+bk31_CL_8.*t))+(6./32).*(1+k)./(1+k+bk32_CL_8.*t).*...
    exp(-k.*(bk32_CL_8).*t./(1+k+bk32_CL_8.*t))+(21./32).*(1+k)./(1+k+bk33_CL_8.*t).*...
    exp(-k.*(bk33_CL_8).*t./(1+k+bk33_CL_8.*t));
%m=16的系数
bk31_CL_16=(sin(pi./16)).^2./(sin(1.*15.*pi./64)).^2;
bk32_CL_16=(sin(pi./16)).^2./(sin(3.*15.*pi./32)).^2;
bk33_CL_16=(sin(pi./16)).^2./(sin(4.*15.*pi./32)).^2;
%    m=16，k=2
MGF_CHIANI_LB216=@(k) (10./64).*(1+k)./(1+k+bk31_CL_16.*t).*...
    exp(-k.*(bk31_CL_16).*t./(1+k+bk31_CL_16.*t))+(12./64).*(1+k)./(1+k+bk32_CL_16.*t).*...
    exp(-k.*(bk32_CL_16).*t./(1+k+bk32_CL_16.*t))+(12./64).*(1+k)./(1+k+bk33_CL_16.*t).*...
    exp(-k.*(bk33_CL_16).*t./(1+k+bk33_CL_16.*t));
%    m=16，k=10 有点问题
MGF_CHIANI_LB1016=@(k) (10./64).*(1+k)./(1+k+bk31_CL_16.*t).*...
    exp(-k.*(bk31_CL_16).*t./(1+k+bk31_CL_16.*t))+(10./64).*(1+k)./(1+k+bk32_CL_16.*t).*...
    exp(-k.*(bk32_CL_16).*t./(1+k+bk32_CL_16.*t))+(13./64).*(1+k)./(1+k+bk33_CL_16.*t).*...
    exp(-k.*(bk33_CL_16).*t./(1+k+bk33_CL_16.*t));

%% 画图
%exact by numerical integration 
    p5=semilogy(x,Pes(2,8),'r');
hold
    semilogy(x,Pes(2,16),'r');
    semilogy(x,Pes(10,8),'r');
    semilogy(x,Pes(10,16),'r');
%our 2−term lower bound
    p1=semilogy(x,MGF_OUR2_28(2),'k.-','MarkerSize',10);
    semilogy(x,MGF_OUR2_108(10),'k.-','MarkerSize',10);
    semilogy(x,MGF_OUR2_216(2),'k.-','MarkerSize',10);
    semilogy(x,MGF_OUR2_1016(10),'k.-','MarkerSize',10);
    %   our 3−term lower bound
    p2=semilogy(x,MGF_OUR3(2,8),'k');
    semilogy(x,MGF_OUR3(2,16),'k');
    semilogy(x,MGF_OUR3(10,8),'k');
    semilogy(x,MGF_OUR3(10,16),'k');
    %    Chiani 3−term lower bound
    p3=semilogy(x,MGF_CHIANI_LB28(2),'b--','linewidth',1.0);
    semilogy(x,MGF_CHIANI_LB108(10),'b--','linewidth',1.0);
    semilogy(x,MGF_CHIANI_LB216(2),'b--','linewidth',1.0);
    semilogy(x,MGF_CHIANI_LB1016(10),'b--','linewidth',1.0);
    %   Chiani 3−term upper bound
    p4=semilogy(x,MGF_CHIANI_UB3(2,8),'b:','linewidth',1.2);
    semilogy(x,MGF_CHIANI_UB3(2,16),'b:','linewidth',1.2);
    semilogy(x,MGF_CHIANI_UB3(10,8),'b:','linewidth',1.2);
    semilogy(x,MGF_CHIANI_UB3(10,16),'b:','linewidth',1.2);
legend([p1 p2 p3 p4 p5],{'our 2-term lower bound','our 3−term lower bound',...
    'Chiani 3-term lower bound','Chiani 3−term upper bound',...
    'exact by numerical integration'},'Location',"southwest",'fontsize',14);

xlabel('SNR(dB)γ')
ylabel('SEP')
axis([15 40 10^(-8) 0.3]);