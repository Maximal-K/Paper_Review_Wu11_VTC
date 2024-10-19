%% 参数
clear 
clc
close
syms al;
bt=@(m) sqrt((1-cos(pi./m))./(1+cos(pi./m)));
x=linspace(15,40,25);
t=10.^(x./10);
linewid=1;

%% exact by numerical integration 通式无误
s=@(m)(-sin(pi./m).^2./(1+cos(pi./m).*cos(al)));
MGF=@(k,m)(1+k)./((1+k+(-s(m)).*t).*pi).*...
    exp(-k.*(-s(m)).*t./(1+k+(-s(m)).*t));  
Pes=@(k,m)int(MGF(k,m),al,0,pi-pi./m); 


%% Chiani 2−term lower bound 通式无误
alk2_CL=@(m) (pi-pi./m)./3;
ak2_CL=@(m) alk2_CL(m)./pi;
bk21_CL=@(m)sin(pi./m).^2./(1+cos(pi./m).*cos(alk2_CL(m)));
bk22_CL=@(m)sin(pi./m).^2./(1+cos(pi./m).*cos(2.*alk2_CL(m)));

MGF_CHIANI_LB2=@(k,m) (ak2_CL(m)).*(1+k)./(1+k+bk21_CL(m).*t).*...
    exp(-k.*(bk21_CL(m)).*t./(1+k+bk21_CL(m).*t))+(ak2_CL(m)).*(1+k)./(1+k+bk22_CL(m).*t).*...
    exp(-k.*(bk22_CL(m)).*t./(1+k+bk22_CL(m).*t)); 

%% Chiani 2−term upper bound 通式无误
alk2_CU=@(m) (pi-pi./m)./3;
ak2_CU_21=1./2;
ak2_CU_22=@(m) 3.*alk2_CU(m)./pi-1./2;
bk21_CU=@(m)sin(pi./m).^2./(1+cos(pi./m).*cos(0));
bk22_CU=@(m)sin(pi./m).^2./(1+cos(pi./m).*cos(1.*pi./2));
MGF_CHIANI_UB2=@(k,m) (ak2_CU_21).*(1+k)./(1+k+bk21_CU(m).*t).*...
    exp(-k.*(bk21_CU(m)).*t./(1+k+bk21_CU(m).*t))+(ak2_CU_22(m)).*(1+k)./(1+k+bk22_CU(m).*t).*...
    exp(-k.*(bk22_CU(m)).*t./(1+k+bk22_CU(m).*t)); 

%% our 2−term lower bound

alk2=@(m) (pi-pi./m)./3;
ak=@(m) alk2(m)./pi;
bk21=@(m)2.*sin(pi./m).*(atan(tan(1.*alk2(m)./2).*bt(m))-...
    atan(tan(0.*alk2(m)./2).*bt(m)))./(alk2(m));
bk22=@(m)2.*sin(pi./m).*(atan(tan(2.*alk2(m)./2).*bt(m))-...
    atan(tan(1.*alk2(m)./2).*bt(m)))./(alk2(m));
MGF_OUR_2=@(k,m) (ak(m)).*(1+k)./(1+k+bk21(m).*t).*...
    exp(-k.*(bk21(m)).*t./(1+k+bk21(m).*t))+(ak(m)).*(1+k)./(1+k+bk22(m).*t).*...
    exp(-k.*(bk22(m)).*t./(1+k+bk22(m).*t)); 
MGF_OUR_108=@(k,m) (6./24).*(1+k)./(1+k+bk21(m).*t).*...
    exp(-k.*(bk21(m)).*t./(1+k+bk21(m).*t))+(9./24).*(1+k)./(1+k+bk22(m).*t).*...
    exp(-k.*(bk22(m)).*t./(1+k+bk22(m).*t)); 

MGF_OUR_28=@(k,m) (8./24).*(1+k)./(1+k+bk21(m).*t).*...
    exp(-k.*(bk21(m)).*t./(1+k+bk21(m).*t))+(7./24).*(1+k)./(1+k+bk22(m).*t).*...
    exp(-k.*(bk22(m)).*t./(1+k+bk22(m).*t)); 
MGF_OUR_216=@(k,m) (8./24).*(1+k)./(1+k+bk21(m).*t).*...
    exp(-k.*(bk21(m)).*t./(1+k+bk21(m).*t))+(7./24).*(1+k)./(1+k+bk22(m).*t).*...
    exp(-k.*(bk22(m)).*t./(1+k+bk22(m).*t)); 


%% 画画
%exact by numerical integration 
p4=semilogy(x,Pes(2,8),'r');
hold on
   semilogy(x,Pes(2,16),'r');
   semilogy(x,Pes(10,8),'r');
semilogy(x,Pes(10,16),'r');
% Chiani 2−term upper bound 通式无误
p3=semilogy(x,MGF_CHIANI_UB2(2,8),'b:','linewidth',1);
 semilogy(x,MGF_CHIANI_UB2(2,16),'b:','linewidth',1);
 semilogy(x,MGF_CHIANI_UB2(10,8),'b:','linewidth',1);
semilogy(x,MGF_CHIANI_UB2(10,16),'b:','linewidth',1);
% Chiani 2−term lower bound 通式无误
p2=semilogy(x,MGF_CHIANI_LB2(2,8),'b--','linewidth',1);
 semilogy(x,MGF_CHIANI_LB2(2,16),'b--','linewidth',1);
 semilogy(x,MGF_CHIANI_LB2(10,8),'b--','linewidth',1);
 semilogy(x,MGF_CHIANI_LB2(10,16),'b--','linewidth',1);
% our 2−term lower bound
p1=semilogy(x,MGF_OUR_28(2,8),'k.-','MarkerSize',10);
semilogy(x,MGF_OUR_216(2,16),'k.-','MarkerSize',10);
semilogy(x,MGF_OUR_108(10,8),'k.-','MarkerSize',10);
semilogy(x,MGF_OUR_2(10,16),'k.-','MarkerSize',10);
legend([p1 p2 p3 p4],{'our 2-term lower bound',...
    'Chiani 2-term lower bound','Chiani 2−term upper bound',...
    'exact by numerical integration'},'Location',"southwest",'fontsize',13);

set(gca,'fontsize',17,'fontname','Times');
xlabel('SNR(dB)γ');
ylabel('SEP');
axis([15 40 10^(-7) 0.3]);