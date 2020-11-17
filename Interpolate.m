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

Fx_deg = T{:,11};
Fy_deg = T{:,14};
Fz_deg = T{:,17};
Mx_deg = T{:,20};
My_deg = T{:,23};
Mz_deg = T{:,26};


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

FFxd=scatteredInterpolant(P,Fx_deg);
FFyd=scatteredInterpolant(P,Fy_deg);
FFzd=scatteredInterpolant(P,Fz_deg);
FMxd=scatteredInterpolant(P,Mx_deg);
FMyd=scatteredInterpolant(P,My_deg);
FMzd=scatteredInterpolant(P,Mz_deg);

VFxq=FFx(Xq,Yq);
VFyq=FFy(Xq,Yq);
VFzq=FFz(Xq,Yq);
VMxq=FMx(Xq,Yq);
VMyq=FMy(Xq,Yq);
VMzq=FMz(Xq,Yq);

VFxd=FFxd(Xq,Yq);
VFyd=FFyd(Xq,Yq);
VFzd=FFzd(Xq,Yq);
VMxd=FMxd(Xq,Yq);
VMyd=FMyd(Xq,Yq);
VMzd=FMzd(Xq,Yq);

mesh(Xq,Yq,VMyd);

% writematrix(VFxd,'Fxd.csv') 
% writematrix(VFzd,'Fzd.csv') 
% writematrix(VMyd,'Myd.csv') 

