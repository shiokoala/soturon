T = readtable("rao.csv");
T(1,:)
% A = table2array(T(:,1));
Fn = T{:,1};
omega = T{:,4};
Fx = T{:,10};
Fy = T{:,13};
Fz = T{:,16};
Mx = T{:,19};
My = T{:,22};
Mz = T{:,25};


P=cat(2,Fn,omega);

fnq= 0.0:0.025:1.0;
oq = 0.0:0.25:15;
[Xq,Yq] = meshgrid(fnq,oq);

FFx=scatteredInterpolant(P,Fx);
FFy=scatteredInterpolant(P,Fy);
FFz=scatteredInterpolant(P,Fz);
FMx=scatteredInterpolant(P,Mx);
FMy=scatteredInterpolant(P,My);
FMz=scatteredInterpolant(P,Mz);


VFxq=FFx(Xq,Yq);
VFyq=FFy(Xq,Yq);
VFzq=FFz(Xq,Yq);
VMxq=FMx(Xq,Yq);
VMyq=FMy(Xq,Yq);
VMzq=FMz(Xq,Yq);

mesh(Xq,Yq,VMzq);

