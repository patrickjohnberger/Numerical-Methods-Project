clear
clc
[A,B,C]=xlsread('PRFData.xlsx');
dH1=A(1);
k01=A(2);
Ea1R=A(3);
dH2=A(4);
k02=A(5);
Ea2R=A(6);
CPa=A(7);
CPb=A(8);
CPc=A(9);
Fa0=A(10);
Fb0=A(11);
Fc0=A(12);
Ca0=A(13);
Cb0=A(14);
Cc0=A(15);
T0=A(16);
Ua=A(17);
Ta=A(18);

V0=0;
Vf=1;
dV=.05;
V=[V0:dV:Vf];

CT0=Ca0+Cb0+Cc0;

Ca=@(Fa,Fb,Fc,T) CT0*(Fa/(Fa+Fb+Fc))*(T0/T);
Cb=@(Fa,Fb,Fc,T) CT0*(Fb/(Fa+Fb+Fc))*(T0/T);
Cc=@(Fa,Fb,Fc,T) CT0*(Fc/(Fa+Fb+Fc))*(T0/T);
k1=@(T) k01*exp(-Ea1R/(T));
k2=@(T) k02*exp(-Ea2R/(T));
r1=@(Fa,Fb,Fc,T)k1(T)*Ca(Fa,Fb,Fc,T);
r2=@(Fa,Fb,Fc,T)k2(T)*(Ca(Fa,Fb,Fc,T))^2;

dFadV=@(V,Fa,Fb,Fc,T) -r1(Fa,Fb,Fc,T)-2*r2(Fa,Fb,Fc,T);
dFbdV=@(V,Fa,Fb,Fc,T) r1(Fa,Fb,Fc,T);
dFcdV=@(V,Fa,Fb,Fc,T) r2(Fa,Fb,Fc,T);
dTdV=@(V,Fa,Fb,Fc,T) ((Ua*(Ta-T)+...
    r1(Fa,Fb,Fc,T)*(-dH1)+...
    r2(Fa,Fb,Fc,T)*(-dH2)))/(Fa*CPa+Fb*CPb+Fc*CPc);

Fa=zeros(1,length(V));
Fb=Fa;
Fc=Fa;
T=Fa;
Fa(1)=Fa0;
Fb(1)=Fb0;
Fc(1)=Fc0;
T(1)=T0;
for i=2:length(V)
   j=i-1;
    Fanew=Fa(j)+dFadV(V(j),Fa(j),Fb(j),Fc(j),T(j))*dV;
    Fbnew=Fb(j)+dFbdV(V(j),Fa(j),Fb(j),Fc(j),T(j))*dV;
    Fcnew=Fc(j)+dFcdV(V(j),Fa(j),Fb(j),Fc(j),T(j))*dV;
    Tnew=T(j)+dTdV(V(j),Fa(j),Fb(j),Fc(j),T(j))*dV;
    dFa=(dFadV(V(j),Fa(j),Fb(j),Fc(j),T(j))...
        +dFadV(V(i),Fanew,Fbnew,Fcnew,Tnew))/2;
    dFb=(dFbdV(V(j),Fa(j),Fb(j),Fc(j),T(j))...
        +dFbdV(V(i),Fanew,Fbnew,Fcnew,Tnew))/2;
    dFc=(dFcdV(V(j),Fa(j),Fb(j),Fc(j),T(j))...
        +dFcdV(V(i),Fanew,Fbnew,Fcnew,Tnew))/2;
    dT=(dTdV(V(j),Fa(j),Fb(j),Fc(j),T(j))...
        +dTdV(V(i),Fanew,Fbnew,Fcnew,Tnew))/2;
    Fa(i)=Fa(j)+dFa*dV;
    Fb(i)=Fb(j)+dFb*dV;
    Fc(i)=Fc(j)+dFc*dV;
    T(i)=T(j)+dT*dV;
end


figure(1)
plot(V,T)
figure(2)
plot(V,Fa,'r',V,Fb,'b',V,Fc,'g')
