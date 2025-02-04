clc; clear; close all;


%% Synthesis
%% i = 1
% indeterminate 
pvar x y v theta
z = [x; y; v; theta];
prog1 = sosprogram(z);
mx = [monomials(z,0:4)];

zdc = [46;34;.2;0];


% dynamics
costheta = 1-theta^2/2 ;
sintheta = theta - theta^3/factorial(3);


fz = [v*costheta; v*sintheta; 0; 0];
g = [0 0; 0 0; 0 1; 1 0];
[prog1, u1] = sospolyvar(prog1,mx(1:end));
[prog1, u2] = sospolyvar(prog1,mx(1:end));
uz = [u1; u2];
zdot = fz + g*uz;
C = [1 0; -1 0; 0 1; 0 -1];
c =  [1; 1; 2; 2];

% HOCBF 1
r1 = 7;
xo1 = 35;
yo1 = 25;
psi1_0 = (x-xo1)^2 + (y-yo1)^2 - r1^2;

% HOCBF 2
r2 = 3;
xo2 = 41;
yo2 = 10;
psi2_0 = (x-xo2)^2 + (y-yo2)^2 - r2^2;

% HOCBF 3
r3 = 4;
xo3 = 10;
yo3 = 40;
psi3_0 = (x-xo3)^2 + (y-yo3)^2 - r3^2;

% HOCBF 4
xo4 = 0;
psi4_0 = x-xo4;

% HOCBF 5
xo5 = 50;
psi5_0 = xo5-x;

% HOCBF 6
yo6 = 0;
psi6_0 = y-yo6;


% HOCBF 7
yo7 = 50;
psi7_0 = yo7-y;
dotpsi7_0 = [0,-1,0,0]*fz;
ddotpsi7_0 = [0, 0, -sintheta, -v*costheta]*zdot;

Vxz = 2*(x-zdc(1))*v*costheta + 1*(x-zdc(1))^2;
Vxdotz = [2*v*costheta+2*(x-zdc(1)), 0, 2*(x-zdc(1))*costheta, -2*(x-zdc(1))*v*sintheta]*zdot(1:4,1);
gammax = 1;

Vyz = 2*(y-zdc(2))*v*sintheta + 1*(y-zdc(2))^2;
Vydotz = [0, 2*v*sintheta+2*(y-zdc(2)), 2*(y-zdc(2))*sintheta, 2*(y-zdc(2))*v*costheta]*zdot(1:4,1);
gammay = 1;

Vvz = (v-zdc(3))^2;
Vvdotz = [0, 0, 2*(v-zdc(3)), 0]*zdot(1:4,1);
gammav = 1;

% SOS
[prog1, s1_0] = sossosvar(prog1,mx);
[prog1, s1_1] = sossosvar(prog1,mx);
[prog1, s1_2] = sossosvar(prog1,mx);
[prog1, s1_3] = sossosvar(prog1,mx);
[prog1, s1_4] = sossosvar(prog1,mx);
[prog1, s1_5] = sossosvar(prog1,mx);
[prog1, s1_6] = sossosvar(prog1,mx);
[prog1, s1_7] = sossosvar(prog1,mx);
[prog1, s1_8] = sossosvar(prog1,mx);
[prog1, s1_9] = sossosvar(prog1,mx);
[prog1, s1_10] = sossosvar(prog1,mx);
[prog1, s1_11] = sossosvar(prog1,mx);

[prog1, su_0] = sossosvar(prog1,mx);
[prog1, su_1] = sossosvar(prog1,mx);
[prog1, su_2] = sossosvar(prog1,mx);
[prog1, su_3] = sossosvar(prog1,mx);
[prog1, su_4] = sossosvar(prog1,mx);
% [prog1, su_5] = sossosvar(prog1,mx);
% [prog1, su_6] = sossosvar(prog1,mx);
% [prog1, su_7] = sossosvar(prog1,mx);
% [prog1, su_8] = sossosvar(prog1,mx);

[prog1, su_9] = sossosvar(prog1,mx);
[prog1, su_10] = sossosvar(prog1,mx);
[prog1, su_11] = sossosvar(prog1,mx);
[prog1, su_12] = sossosvar(prog1,mx);
% [prog1, su_13] = sossosvar(prog1,mx);
% [prog1, su_14] = sossosvar(prog1,mx);
% [prog1, su_15] = sossosvar(prog1,mx);
% [prog1, su_16] = sossosvar(prog1,mx);

[prog1, su_17] = sossosvar(prog1,mx);
[prog1, su_18] = sossosvar(prog1,mx);
[prog1, su_19] = sossosvar(prog1,mx);
[prog1, su_20] = sossosvar(prog1,mx);
% [prog1, su_21] = sossosvar(prog1,mx);
% [prog1, su_22] = sossosvar(prog1,mx);
% [prog1, su_23] = sossosvar(prog1,mx);
% [prog1, su_24] = sossosvar(prog1,mx);

[prog1, su_25] = sossosvar(prog1,mx);
[prog1, su_26] = sossosvar(prog1,mx);
[prog1, su_27] = sossosvar(prog1,mx);
[prog1, su_28] = sossosvar(prog1,mx);
% [prog1, su_29] = sossosvar(prog1,mx);
% [prog1, su_30] = sossosvar(prog1,mx);
% [prog1, su_31] = sossosvar(prog1,mx);
% [prog1, su_32] = sossosvar(prog1,mx);
[prog1, su_33] = sossosvar(prog1,mx);
[prog1, su_34] = sossosvar(prog1,mx);
[prog1, su_35] = sossosvar(prog1,mx);

[prog1, su_36] = sossosvar(prog1,mx);
[prog1, su_37] = sossosvar(prog1,mx);
[prog1, su_38] = sossosvar(prog1,mx);
[prog1, su_39] = sossosvar(prog1,mx);
[prog1, su_40] = sossosvar(prog1,mx);
[prog1, su_41] = sossosvar(prog1,mx);
[prog1, su_42] = sossosvar(prog1,mx);
[prog1, su_43] = sossosvar(prog1,mx);
[prog1, su_44] = sossosvar(prog1,mx);
[prog1, su_45] = sossosvar(prog1,mx);
[prog1, su_46] = sossosvar(prog1,mx);
[prog1, su_47] = sossosvar(prog1,mx);
[prog1, su_48] = sossosvar(prog1,mx);
[prog1, su_49] = sossosvar(prog1,mx);
[prog1, su_50] = sossosvar(prog1,mx);
[prog1, su_51] = sossosvar(prog1,mx);
[prog1, su_52] = sossosvar(prog1,mx);
[prog1, su_53] = sossosvar(prog1,mx);
[prog1, su_54] = sossosvar(prog1,mx);
[prog1, su_55] = sossosvar(prog1,mx);
[prog1, su_56] = sossosvar(prog1,mx);
[prog1, su_57] = sossosvar(prog1,mx);
[prog1, su_58] = sossosvar(prog1,mx);
[prog1, su_59] = sossosvar(prog1,mx);
[prog1, su_60] = sossosvar(prog1,mx);
[prog1, su_61] = sossosvar(prog1,mx);
[prog1, su_62] = sossosvar(prog1,mx);
[prog1, su_63] = sossosvar(prog1,mx);


[prog1, sV_0] = sossosvar(prog1,mx);
[prog1, sV_1] = sossosvar(prog1,mx);
[prog1, sV_2] = sossosvar(prog1,mx);
[prog1, sV_3] = sossosvar(prog1,mx);
[prog1, sV_4] = sossosvar(prog1,mx);
[prog1, sV_5] = sossosvar(prog1,mx);
% [prog1, sV_6] = sossosvar(prog1,mx);
% [prog1, sV_7] = sossosvar(prog1,mx);
% [prog1, sV_8] = sossosvar(prog1,mx);
% [prog1, sV_9] = sossosvar(prog1,mx);
[prog1, sV_10] = sossosvar(prog1,mx);
[prog1, sV_11] = sossosvar(prog1,mx);
[prog1, sV_12] = sossosvar(prog1,mx);
[prog1, sV_13] = sossosvar(prog1,mx);
[prog1, sV_14] = sossosvar(prog1,mx);
[prog1, sV_15] = sossosvar(prog1,mx);
% [prog1, sV_16] = sossosvar(prog1,mx);
% [prog1, sV_17] = sossosvar(prog1,mx);
% [prog1, sV_18] = sossosvar(prog1,mx);
% [prog1, sV_19] = sossosvar(prog1,mx);
[prog1, sV_20] = sossosvar(prog1,mx);
[prog1, sV_21] = sossosvar(prog1,mx);
[prog1, sV_22] = sossosvar(prog1,mx);
[prog1, sV_23] = sossosvar(prog1,mx);
[prog1, sV_24] = sossosvar(prog1,mx);
[prog1, sV_25] = sossosvar(prog1,mx);
[prog1, sV_26] = sossosvar(prog1,mx);
[prog1, sV_27] = sossosvar(prog1,mx);
[prog1, sV_28] = sossosvar(prog1,mx);
[prog1, sV_29] = sossosvar(prog1,mx);
[prog1, sV_30] = sossosvar(prog1,mx);
[prog1, sV_31] = sossosvar(prog1,mx);
[prog1, sV_32] = sossosvar(prog1,mx);
[prog1, sV_33] = sossosvar(prog1,mx);
[prog1, sV_34] = sossosvar(prog1,mx);
[prog1, sV_35] = sossosvar(prog1,mx);
[prog1, sV_36] = sossosvar(prog1,mx);
[prog1, sV_37] = sossosvar(prog1,mx);
[prog1, sV_38] = sossosvar(prog1,mx);
[prog1, sV_39] = sossosvar(prog1,mx);
[prog1, sV_40] = sossosvar(prog1,mx);
[prog1, sV_41] = sossosvar(prog1,mx);
[prog1, sV_42] = sossosvar(prog1,mx);
[prog1, sV_43] = sossosvar(prog1,mx);



Zmin = [0, 0, 0, -3.14];
Zmax = [50, 50, 2, 3.14];

dpvar b7_1 
beta7_1prime = b7_1;
prog1 = sosdecvar(prog1,b7_1);
prog1 = sosineq(prog1,beta7_1prime);

dpvar  rhox
prog1 = sosdecvar(prog1,rhox);
prog1 = sosineq(prog1,rhox);
CLlFcertx = -Vxdotz + rhox - sV_0 - sV_1*psi7_0 - sV_2*(v-Zmin(3)) - sV_3*(theta-Zmin(4)) - sV_4*(Zmax(3)-v) - sV_5*(Zmax(4)-theta) - sV_26*psi1_0 - sV_27*psi2_0 - sV_28*psi3_0 - sV_29*psi4_0 - sV_30*psi5_0 - sV_31*psi6_0;
prog1 = soseq(prog1, CLlFcertx);

dpvar  rhoy
prog1 = sosdecvar(prog1,rhoy);
prog1 = sosineq(prog1,rhoy);
CLlFcerty = -Vydotz + rhoy - sV_20 - sV_21*psi7_0 - sV_22*(v-Zmin(3)) - sV_23*(theta-Zmin(4)) - sV_24*(Zmax(3)-v) - sV_25*(Zmax(4)-theta) - sV_32*psi1_0 - sV_33*psi2_0 - sV_34*psi3_0 - sV_35*psi4_0 - sV_36*psi5_0 - sV_37*psi6_0;
prog1 = soseq(prog1, CLlFcerty);

dpvar  rhov
prog1 = sosdecvar(prog1,rhov);
prog1 = sosineq(prog1,rhov);
CLlFcertv = -Vvdotz + rhov - sV_10 - sV_11*psi7_0 - sV_12*(v-Zmin(3)) - sV_13*(theta-Zmin(4)) - sV_14*(Zmax(3)-v) - sV_15*(Zmax(4)-theta) - sV_38*psi1_0 - sV_39*psi2_0 - sV_40*psi3_0 - sV_41*psi4_0 - sV_42*psi5_0 - sV_43*psi6_0;
prog1 = soseq(prog1, CLlFcertv);

CBFcert1 = ddotpsi7_0 - beta7_1prime - s1_0 - s1_1*psi7_0 - s1_2*(v-Zmin(3)) - s1_3*(theta-Zmin(4)) - s1_4*(Zmax(3)-v) - s1_5*(Zmax(4)-theta) - s1_6*psi1_0 - s1_7*psi2_0 - s1_8*psi3_0 - s1_9*psi4_0 - s1_10*psi5_0 - s1_11*psi6_0;
prog1 = soseq(prog1, CBFcert1);

inputreq = -C*uz+c;
inputcert1 = inputreq(1) - su_0 - su_1*(v-Zmin(3)) - su_2*(theta-Zmin(4)) - su_3*(Zmax(3)-v) - su_4*(Zmax(4)-theta) - su_36*psi1_0 - su_37*psi2_0 - su_38*psi3_0 - su_39*psi4_0 - su_40*psi5_0 - su_41*psi6_0 - su_42*psi7_0;
prog1 = soseq(prog1, inputcert1);
inputcert2 = inputreq(2) - su_33 - su_9*(v-Zmin(3)) - su_10*(theta-Zmin(4)) - su_11*(Zmax(3)-v) - su_12*(Zmax(4)-theta) - su_43*psi1_0 - su_44*psi2_0 - su_45*psi3_0 - su_46*psi4_0 - su_47*psi5_0 - su_48*psi6_0 - su_49*psi7_0;
prog1 = soseq(prog1, inputcert2);
inputcert3 = inputreq(3) - su_34 -  su_17*(v-Zmin(3)) - su_18*(theta-Zmin(4)) - su_19*(Zmax(3)-v) - su_20*(Zmax(4)-theta) - su_50*psi1_0 - su_51*psi2_0 - su_52*psi3_0 - su_53*psi4_0 - su_54*psi5_0 - su_55*psi6_0 - su_56*psi7_0;
prog1 = soseq(prog1, inputcert3);
inputcert4 = inputreq(4) - su_35  - su_25*(v-Zmin(3)) - su_26*(theta-Zmin(4)) - su_27*(Zmax(3)-v) - su_28*(Zmax(4)-theta) - su_57*psi1_0 - su_58*psi2_0 - su_59*psi3_0 - su_60*psi4_0 - su_61*psi5_0 - su_62*psi6_0 - su_63*psi7_0;
prog1 = soseq(prog1, inputcert4);


% Solve

prog1 = sossetobj(prog1,-b7_1)
prog1 = sossolve(prog1)

fprintf('Feasibility Ratio: %g\n',prog1.solinfo.info.feasratio);
fprintf('Numerical issues: %g\n',prog1.solinfo.info.numerr);
fprintf('Primal infeasible? %g\n',prog1.solinfo.info.pinf);

b7_1sol = double(sosgetsol(prog1,b7_1))
rhoxsol = double(sosgetsol(prog1,rhox))
rhoysol = double(sosgetsol(prog1,rhoy))
rhovsol = double(sosgetsol(prog1,rhov))
uzsol = sosgetsol(prog1,uz)



%% i=2
% clear prog1 z

% zeta1 max
fun = @(z) -1*(yo7-z(2));
function [c,ceq] = psi0(z)
ceq = [];
c = -1*(50-z(2));
end
nonlcon = @psi0
[zmax,zeta7max] = fmincon(fun, [0;0],[],[],[],[],[Zmin(1);Zmin(2)],[Zmax(1),Zmax(2)],nonlcon);
zeta7bar = -1*zeta7max;
beta7bar = b7_1sol*zeta7bar;


% indeterminate 
pvar x y v theta
z = [x; y; v; theta];
prog2 = sosprogram(z);
mx = [monomials(z,0:4)];


% input
[prog2, u1] = sospolyvar(prog2,mx(1:end));
[prog2, u2] = sospolyvar(prog2,mx(1:end));
uz = [u1; u2];


% dynamics
costheta = 1-theta^2/2 ;
sintheta = theta - theta^3/factorial(3);
cosprime = - theta ;
sinprime = 1 - theta^2/2;
fz = [v*costheta; v*sintheta; 0; 0];
g = [0 0; 0 0; 0 1; 1 0];
zdot = fz + g*uz;


% SOS
[prog2, s2_0] = sossosvar(prog2,mx);
[prog2, s2_1] = sossosvar(prog2,mx);
[prog2, s2_2] = sossosvar(prog2,mx);
[prog2, s2_3] = sossosvar(prog2,mx);
[prog2, s2_4] = sossosvar(prog2,mx);
[prog2, s2_5] = sossosvar(prog2,mx);
[prog2, s2_6] = sossosvar(prog2,mx);
[prog2, s2_7] = sossosvar(prog2,mx);
[prog2, s2_8] = sossosvar(prog2,mx);
[prog2, s2_9] = sossosvar(prog2,mx);
[prog2, s2_10] = sossosvar(prog2,mx);
[prog2, s2_11] = sossosvar(prog2,mx);


[prog2, sV_0] = sossosvar(prog2,mx);
[prog2, sV_1] = sossosvar(prog2,mx);
[prog2, sV_2] = sossosvar(prog2,mx);
[prog2, sV_3] = sossosvar(prog2,mx);
[prog2, sV_4] = sossosvar(prog2,mx);
[prog2, sV_5] = sossosvar(prog2,mx);
% [prog2, sV_6] = sossosvar(prog2,mx);
% [prog2, sV_7] = sossosvar(prog2,mx);
% [prog2, sV_8] = sossosvar(prog2,mx);
% [prog2, sV_9] = sossosvar(prog2,mx);
% [prog2, sV_10] = sossosvar(prog2,mx);
[prog2, sV_11] = sossosvar(prog2,mx);
[prog2, sV_12] = sossosvar(prog2,mx);
[prog2, sV_13] = sossosvar(prog2,mx);
[prog2, sV_14] = sossosvar(prog2,mx);
[prog2, sV_15] = sossosvar(prog2,mx);
% [prog2, sV_16] = sossosvar(prog2,mx);
% [prog2, sV_17] = sossosvar(prog2,mx);
% [prog2, sV_18] = sossosvar(prog2,mx);
% [prog2, sV_19] = sossosvar(prog2,mx);
% [prog2, sV_20] = sossosvar(prog2,mx);
[prog2, sV_21] = sossosvar(prog2,mx);
[prog2, sV_22] = sossosvar(prog2,mx);
[prog2, sV_23] = sossosvar(prog2,mx);
[prog2, sV_24] = sossosvar(prog2,mx);
[prog2, sV_25] = sossosvar(prog2,mx);
[prog2, sV_26] = sossosvar(prog2,mx);
[prog2, sV_27] = sossosvar(prog2,mx);
[prog2, sV_28] = sossosvar(prog2,mx);
[prog2, sV_29] = sossosvar(prog2,mx);
[prog2, sV_30] = sossosvar(prog2,mx);
[prog2, sV_31] = sossosvar(prog2,mx);
[prog2, sV_32] = sossosvar(prog2,mx);
[prog2, sV_33] = sossosvar(prog2,mx);
[prog2, sV_34] = sossosvar(prog2,mx);
[prog2, sV_35] = sossosvar(prog2,mx);
[prog2, sV_36] = sossosvar(prog2,mx);
[prog2, sV_37] = sossosvar(prog2,mx);
[prog2, sV_38] = sossosvar(prog2,mx);
[prog2, sV_39] = sossosvar(prog2,mx);
[prog2, sV_40] = sossosvar(prog2,mx);
[prog2, sV_41] = sossosvar(prog2,mx);
[prog2, sV_42] = sossosvar(prog2,mx);
[prog2, sV_43] = sossosvar(prog2,mx);



[prog2, su_0] = sossosvar(prog2,mx);
[prog2, su_1] = sossosvar(prog2,mx);
[prog2, su_2] = sossosvar(prog2,mx);
[prog2, su_3] = sossosvar(prog2,mx);
[prog2, su_4] = sossosvar(prog2,mx);
% [prog2, su_5] = sossosvar(prog2,mx);
% [prog2, su_6] = sossosvar(prog2,mx);
% [prog2, su_7] = sossosvar(prog2,mx);
% [prog2, su_8] = sossosvar(prog2,mx);

[prog2, su_9] = sossosvar(prog2,mx);
[prog2, su_10] = sossosvar(prog2,mx);
[prog2, su_11] = sossosvar(prog2,mx);
[prog2, su_12] = sossosvar(prog2,mx);
% [prog2, su_13] = sossosvar(prog2,mx);
% [prog2, su_14] = sossosvar(prog2,mx);
% [prog2, su_15] = sossosvar(prog2,mx);
% [prog2, su_16] = sossosvar(prog2,mx);

[prog2, su_17] = sossosvar(prog2,mx);
[prog2, su_18] = sossosvar(prog2,mx);
[prog2, su_19] = sossosvar(prog2,mx);
[prog2, su_20] = sossosvar(prog2,mx);
% [prog2, su_21] = sossosvar(prog2,mx);
% [prog2, su_22] = sossosvar(prog2,mx);
% [prog2, su_23] = sossosvar(prog2,mx);
% [prog2, su_24] = sossosvar(prog2,mx);

[prog2, su_25] = sossosvar(prog2,mx);
[prog2, su_26] = sossosvar(prog2,mx);
[prog2, su_27] = sossosvar(prog2,mx);
[prog2, su_28] = sossosvar(prog2,mx);
% [prog2, su_29] = sossosvar(prog2,mx);
% [prog2, su_30] = sossosvar(prog2,mx);
% [prog2, su_31] = sossosvar(prog2,mx);
% [prog2, su_32] = sossosvar(prog2,mx);
[prog2, su_33] = sossosvar(prog2,mx);
[prog2, su_34] = sossosvar(prog2,mx);
[prog2, su_35] = sossosvar(prog2,mx);

[prog2, su_36] = sossosvar(prog2,mx);
[prog2, su_37] = sossosvar(prog2,mx);
[prog2, su_38] = sossosvar(prog2,mx);
[prog2, su_39] = sossosvar(prog2,mx);
[prog2, su_40] = sossosvar(prog2,mx);
[prog2, su_41] = sossosvar(prog2,mx);
[prog2, su_42] = sossosvar(prog2,mx);
[prog2, su_43] = sossosvar(prog2,mx);
[prog2, su_44] = sossosvar(prog2,mx);
[prog2, su_45] = sossosvar(prog2,mx);
[prog2, su_46] = sossosvar(prog2,mx);
[prog2, su_47] = sossosvar(prog2,mx);
[prog2, su_48] = sossosvar(prog2,mx);
[prog2, su_49] = sossosvar(prog2,mx);
[prog2, su_50] = sossosvar(prog2,mx);
[prog2, su_51] = sossosvar(prog2,mx);
[prog2, su_52] = sossosvar(prog2,mx);
[prog2, su_53] = sossosvar(prog2,mx);
[prog2, su_54] = sossosvar(prog2,mx);
[prog2, su_55] = sossosvar(prog2,mx);
[prog2, su_56] = sossosvar(prog2,mx);
[prog2, su_57] = sossosvar(prog2,mx);
[prog2, su_58] = sossosvar(prog2,mx);
[prog2, su_59] = sossosvar(prog2,mx);
[prog2, su_60] = sossosvar(prog2,mx);
[prog2, su_61] = sossosvar(prog2,mx);
[prog2, su_62] = sossosvar(prog2,mx);
[prog2, su_63] = sossosvar(prog2,mx);



% HOCBF 1
b1_1sol = 0.1953;
a1_1sol = 1.0889;
A1 = 0.0147;
psi1_0 = (x-xo1)^2 + (y-yo1)^2 - r1^2;
alpha1 = A1*psi1_0;
psi1_1 = 2*v*(x-xo1)*costheta + 2*v*(y-yo1)*sintheta + alpha1;
dotpsi1_1 = [2*v*costheta + 2*A1*(x-xo1), 2*v*sintheta + 2*A1*(y-yo1), 2*(x-xo1)*costheta + 2*(y-yo1)*sintheta, -2*v*(x-xo1)*sintheta + 2*v*(y-yo1)*costheta]*zdot;

% HOCBF 2
b2_1sol = 34.3503;
a2_1sol = 1.1963e+04;
A2 = 0.1969;
psi2_0 = (x-xo2)^2 + (y-yo2)^2 - r2^2;
alpha2 = A2*psi2_0;
psi2_1 = 2*v*(x-xo2)*costheta + 2*v*(y-yo2)*sintheta + alpha2;
dotpsi2_1 = [2*v*costheta + 2*A2*(x-xo2), 2*v*sintheta + 2*A2*(y-yo2), 2*(x-xo2)*costheta + 2*(y-yo2)*sintheta, -2*v*(x-xo2)*sintheta + 2*v*(y-yo2)*costheta]*zdot;

% HOCBF 3
b3_1sol = 0.0152;
a3_1sol = 1;
A3 = 0.0042;
psi3_0 = (x-xo3)^2 + (y-yo3)^2 - r3^2;
alpha3 = A3*psi3_0;
psi3_1 = 2*v*(x-xo3)*costheta + 2*v*(y-yo3)*sintheta + alpha3;
dotpsi3_1 = [2*v*costheta + 2*A3*(x-xo3), 2*v*sintheta + 2*A3*(y-yo3), 2*(x-xo3)*costheta + 2*(y-yo3)*sintheta, -2*v*(x-xo3)*sintheta + 2*v*(y-yo3)*costheta]*zdot;

% HOCBF 4
b4_1sol = 1.4277;
a4_1sol = 0.5963;
A4 = 0.2390;
psi4_0 = x-xo4;
alpha4 = A4*psi4_0;
psi4_1 = v*costheta + alpha4;
dotpsi4_1 = [A4, 0, costheta, -v*sintheta]*zdot;

% HOCBF 5
b5_1sol = 1.4830;
a5_1sol = 1;
A5 = 0.2436;
psi5_0 = xo5-x;
alpha5 = A5*psi5_0;
psi5_1 = -v*costheta + alpha5;
dotpsi5_1 = [-A5, 0, -costheta, v*sintheta]*zdot;

% HOCBF 6
b6_1sol = 0.0052;
a6_1sol = 7.3323;
A6 = 0.0144;
psi6_0 = y-yo6;
alpha6 = A6*psi6_0;
psi6_1 = v*sintheta + alpha6;
dotpsi6_1 = [0, A6, sintheta, v*costheta]*zdot;

% HOCBF 7
psi7_0 = yo7-y;
A7 = (1/zeta7bar)*sqrt(2*beta7bar);
alpha7 = A7*psi7_0;
psi7_1 = -v*sintheta + alpha7;
dotpsi7_1 = [0, -A7, -sintheta, -v*costheta]*zdot;

Vxz = 2*(x-zdc(1))*v*costheta + 1*(x-zdc(1))^2;
Vxdotz = [2*v*costheta+2*(x-zdc(1)), 0, 2*(x-zdc(1))*costheta, -2*(x-zdc(1))*v*sintheta]*zdot(1:4,1);
gammax = 1;

Vyz = 2*(y-zdc(2))*v*sintheta + 1*(y-zdc(2))^2;
Vydotz = [0, 2*v*sintheta+2*(y-zdc(2)), 2*(y-zdc(2))*sintheta, 2*(y-zdc(2))*v*costheta]*zdot(1:4,1);
gammay = 1;

Vvz = (v-zdc(3))^2;
Vvdotz = [0, 0, 2*(v-zdc(3)), 0]*zdot(1:4,1);
gammav = 1;


% HOCBF 7
dpvar a7_1
alpha7_2 = a7_1*psi7_1;
prog2 = sosdecvar(prog2, a7_1);
prog2 = sosineq(prog2,a7_1);
CBFcert2 = dotpsi7_1 + alpha7_2 - s2_0 - s2_1*psi7_0 - s2_2*(v-Zmin(3)) - s2_3*(theta-Zmin(4)) - s2_4*(Zmax(3)-v) - s2_5*(Zmax(4)-theta) - s2_6*psi1_0 - s2_7*psi2_0 - s2_8*psi3_0 - s2_9*psi4_0 - s2_10*psi5_0 - s2_11*psi6_0;
prog2 = soseq(prog2, CBFcert2);

dpvar  rhov
prog2 = sosdecvar(prog2,rhov);
prog2 = sosineq(prog2,rhov);
CLlFcertv = -Vvdotz + rhov - sV_10 - sV_11*psi7_0 - sV_12*(v-Zmin(3)) - sV_13*(theta-Zmin(4)) - sV_14*(Zmax(3)-v) - sV_15*(Zmax(4)-theta) - sV_26*psi1_0 - sV_27*psi2_0 - sV_28*psi3_0 - sV_29*psi4_0 - sV_30*psi5_0 - sV_31*psi6_0;
prog2 = soseq(prog2, CLlFcertv);
                         
dpvar  rhox
prog2 = sosdecvar(prog2,rhox);
prog2 = sosineq(prog2,rhox);
CLlFcertx = -Vxdotz + rhox - sV_0 - sV_1*psi7_0 - sV_2*(v-Zmin(3)) - sV_3*(theta-Zmin(4)) - sV_4*(Zmax(3)-v) - sV_5*(Zmax(4)-theta) - sV_32*psi1_0 - sV_33*psi2_0 - sV_34*psi3_0 - sV_35*psi4_0 - sV_36*psi5_0 - sV_37*psi6_0;
prog2 = soseq(prog2, CLlFcertx);

dpvar  rhoy
prog2 = sosdecvar(prog2,rhoy);
prog2 = sosineq(prog2,rhoy);
CLlFcerty = -Vydotz + rhoy - sV_20 - sV_21*psi7_0 - sV_22*(v-Zmin(3)) - sV_23*(theta-Zmin(4)) - sV_24*(Zmax(3)-v) - sV_25*(Zmax(4)-theta) - sV_38*psi1_0 - sV_39*psi2_0 - sV_40*psi3_0 - sV_41*psi4_0 - sV_42*psi5_0 - sV_43*psi6_0;
prog2 = soseq(prog2, CLlFcerty);


inputreq = -C*uz+c;
inputcert1 = inputreq(1) - su_0 - su_1*(v-Zmin(3)) - su_2*(theta-Zmin(4)) - su_3*(Zmax(3)-v) - su_4*(Zmax(4)-theta) - su_36*psi1_0 - su_37*psi2_0 - su_38*psi3_0 - su_39*psi4_0 - su_40*psi5_0 - su_41*psi6_0 - su_42*psi7_0;
prog2 = soseq(prog2, inputcert1);
inputcert2 = inputreq(2) - su_33 - su_9*(v-Zmin(3)) - su_10*(theta-Zmin(4)) - su_11*(Zmax(3)-v) - su_12*(Zmax(4)-theta) - su_43*psi1_0 - su_44*psi2_0 - su_45*psi3_0 - su_46*psi4_0 - su_47*psi5_0 - su_48*psi6_0 - su_49*psi7_0;
prog2 = soseq(prog2, inputcert2);
inputcert3 = inputreq(3) - su_34 -  su_17*(v-Zmin(3)) - su_18*(theta-Zmin(4)) - su_19*(Zmax(3)-v) - su_20*(Zmax(4)-theta) - su_50*psi1_0 - su_51*psi2_0 - su_52*psi3_0 - su_53*psi4_0 - su_54*psi5_0 - su_55*psi6_0 - su_56*psi7_0;
prog2 = soseq(prog2, inputcert3);
inputcert4 = inputreq(4) - su_35 - su_25*(v-Zmin(3)) - su_26*(theta-Zmin(4)) - su_27*(Zmax(3)-v) - su_28*(Zmax(4)-theta) - su_57*psi1_0 - su_58*psi2_0 - su_59*psi3_0 - su_60*psi4_0 - su_61*psi5_0 - su_62*psi6_0 - su_63*psi7_0;
prog2 = soseq(prog2, inputcert4);



% Solve
prog2 = sossetobj(prog2,-a7_1)
prog2 = sossolve(prog2)

fprintf('Feasibility Ratio: %g\n',prog2.solinfo.info.feasratio);
fprintf('Numerical issues: %g\n',prog2.solinfo.info.numerr);
fprintf('Primal infeasible? %g\n',prog2.solinfo.info.pinf);

uzsol = sosgetsol(prog2,uz)
a7_1sol = double(sosgetsol(prog2,a7_1))
rhoxsol = double(sosgetsol(prog2,rhox))
rhovsol = double(sosgetsol(prog2,rhov))


%% HOCBF 1
close all
t0 = 0;
tf = 200;
dt = 0.1;
tvec = t0:dt:tf;

% Initial state
x0MC = 5:10:45;
y0MC = 5:10:45;


costheta = @(z) 1-z(4)^2/2 ;
sintheta = @(z) z(4) - z(4)^3/factorial(3) ;



% HOCBF 1
psi1_0 =  @(z,xo,yo) (z(1)-xo)^2 + (z(2)-yo)^2 - 7^2;

dotpsi1_0 = @(z,xo,yo) 2*(z(1)-xo)*z(3)*costheta(z) + 2*(z(2)-yo)*z(3)*sintheta(z);
beta1 = @(z,xo,yo) b1_1sol*psi1_0(z,xo,yo);
alpha1 = @(z,xo,yo) 1*real(sqrt(2*beta1(z,xo,yo)));

psi1_1 = @(z,xo,yo) 2*z(3)*(z(1)-xo)*costheta(z) + 2*z(3)*(z(2)-yo)*sintheta(z) + alpha1(z,xo,yo);
alpha2 = @(z,xo,yo) 1*a1_1sol*psi1_1(z,xo,yo);

A1max = @(z,xo,yo) [2*z(3)*(z(1)-xo)*sintheta(z)-2*z(3)*(z(2)-yo)*costheta(z), -2*(z(1)-xo)*costheta(z)-2*(z(2)-yo)*sintheta(z), 0, 0, 0];
b1max = @(z,xo,yo) z(3)*costheta(z)*[2*z(3)*costheta(z)+(4*b1_1sol*(z(1)-xo))/(2*alpha1(z,xo,yo))] + z(3)*sintheta(z)*[2*z(3)*sintheta(z)+(4*b1_1sol*(z(2)-yo))/(2*alpha1(z,xo,yo))] + alpha2(z,xo,yo);


% HOCBF 2
psi2_0 = @(z,xo,yo) (z(1)-xo)^2 + (z(2)-yo)^2 - r2^2;

dotpsi2_0 = @(z,xo,yo) 2*(z(1)-xo)*z(3)*costheta(z) + 2*(z(2)-yo)*z(3)*sintheta(z);
beta2 = @(z,xo,yo) b2_1sol*psi2_0(z,xo,yo);
alpha2_1 = @(z,xo,yo) 1*real(sqrt(2*beta2(z,xo,yo)));

psi2_1 = @(z,xo,yo) 2*z(3)*(z(1)-xo)*costheta(z) + 2*z(3)*(z(2)-yo)*sintheta(z) + alpha2_1(z,xo,yo);
alpha2_2 = @(z,xo,yo) 1*a2_1sol*psi2_1(z,xo,yo);

A2max = @(z,xo,yo) [2*z(3)*(z(1)-xo)*sintheta(z)-2*z(3)*(z(2)-yo)*costheta(z), -2*(z(1)-xo)*costheta(z)-2*(z(2)-yo)*sintheta(z), 0, 0, 0];
b2max = @(z,xo,yo) z(3)*costheta(z)*[2*z(3)*costheta(z)+(4*b2_1sol*(z(1)-xo))/(2*alpha2_1(z,xo,yo))] + z(3)*sintheta(z)*[2*z(3)*sintheta(z)+(4*b2_1sol*(z(2)-yo))/(2*alpha2_1(z,xo,yo))] + alpha2_2(z,xo,yo);


% HOCBF 3
psi3_0 = @(z,xo,yo) (z(1)-xo)^2 + (z(2)-yo)^2 - r3^2;

dotpsi3_0 = @(z,xo,yo) 2*(z(1)-xo)*z(3)*costheta(z) + 2*(z(2)-yo)*z(3)*sintheta(z);
beta3 = @(z,xo,yo) b3_1sol*psi3_0(z,xo,yo);
alpha3_1 = @(z,xo,yo) 1*real(sqrt(2*beta3(z,xo,yo)));

psi3_1 = @(z,xo,yo) 2*z(3)*(z(1)-xo)*costheta(z) + 2*z(3)*(z(2)-yo)*sintheta(z) + alpha3_1(z,xo,yo);
alpha3_2 = @(z,xo,yo) 1*a3_1sol*psi3_1(z,xo,yo);

A3max = @(z,xo,yo) [2*z(3)*(z(1)-xo)*sintheta(z)-2*z(3)*(z(2)-yo)*costheta(z), -2*(z(1)-xo)*costheta(z)-2*(z(2)-yo)*sintheta(z), 0, 0, 0];
b3max = @(z,xo,yo) z(3)*costheta(z)*[2*z(3)*costheta(z)+(4*b3_1sol*(z(1)-xo))/(2*alpha3_1(z,xo,yo))] + z(3)*sintheta(z)*[2*z(3)*sintheta(z)+(4*b3_1sol*(z(2)-yo))/(2*alpha3_1(z,xo,yo))] + alpha3_2(z,xo,yo);


% HOCBF 4
psi4_0 = @(z,xo) z(1)-xo;

dotpsi4_0 = @(z,xo) z(3)*costheta(z);
beta4 = @(z,xo) b4_1sol*psi4_0(z,xo);
alpha4_1 = @(z,xo) 1*real(sqrt(2*beta4(z,xo)));

psi4_1 = @(z,xo) z(3)*costheta(z) + alpha4_1(z,xo);
alpha4_2 = @(z,xo) 1*a4_1sol*psi4_1(z,xo);

A4max = @(z,xo) [z(3)*sintheta(z), -costheta(z), 0, 0, 0];
b4max = @(z,xo,yo) a4_1sol*z(3)*costheta(z) + alpha4_2(z,xo);


% HOCBF 5
psi5_0 = @(z,xo) xo-z(1);

dotpsi5_0 = @(z,xo) -z(3)*costheta(z);
beta5 = @(z,xo) b5_1sol*psi5_0(z,xo);
alpha5_1 = @(z,xo) 1*real(sqrt(2*beta5(z,xo)));

psi5_1 = @(z,xo) -z(3)*costheta(z) + alpha5_1(z,xo);
alpha5_2 = @(z,xo) 1*a5_1sol*psi5_1(z,xo);

A5max = @(z,xo) [-z(3)*sintheta(z), costheta(z), 0, 0, 0];
b5max = @(z,xo,yo) -a5_1sol*z(3)*costheta(z) + alpha5_2(z,xo);


% HOCBF 6
psi6_0 = @(z,yo) z(2)-yo;

dotpsi6_0 = @(z,yo) z(3)*sintheta(z);
beta6 = @(z,yo) b6_1sol*psi6_0(z,yo);
alpha6_1 = @(z,yo) 1*real(sqrt(2*beta6(z,yo)));

psi6_1 = @(z,yo) z(3)*sintheta(z) + alpha6_1(z,yo);
alpha6_2 = @(z,yo) 1*a6_1sol*psi6_1(z,yo);

A6max = @(z,yo) [-z(3)*costheta(z), -sintheta(z), 0, 0, 0];
b6max = @(z,yo) a6_1sol*z(3)*sintheta(z) + alpha6_2(z,yo);


% HOCBF 7
psi7_0 = @(z,yo) yo-z(2);

dotpsiy_0 = @(z,yo) -z(3)*sintheta(z);
beta7 = @(z,yo) b7_1sol*psi7_0(z,yo);
alpha7_1 = @(z,yo) 1*real(sqrt(2*beta7(z,yo)));

psi7_1 = @(z,yo) -z(3)*sintheta(z) + alpha7_1(z,yo);
alpha7_2 = @(z,yo) 1*a7_1sol*psi7_1(z,yo);

A7max = @(z,yo) [z(3)*costheta(z), sintheta(z), 0, 0, 0];
b7max = @(z,yo) -a7_1sol*z(3)*sintheta(z) + alpha7_2(z,yo);





% CLFs
AV2 = @(z,zdc) [0, 2*(z(3)-zdc(3)), 0, 0, -1];
bV2 = @(z,zdc) -.01*(z(3)-zdc(3))^2;

Atest1 = @(z,zdc) [-2*(z(1)-zdc(1))*z(3)*sintheta(z), 2*(z(1)-zdc(1))*costheta(z), -1, 0, 0];
btest1 = @(z,zdc) -2*z(3)^2*costheta(z)^2 - 1*2*(z(1)-zdc(1))*z(3)*costheta(z) - .1*(2*(z(1)-zdc(1))*z(3)*costheta(z) + 1*(z(1)-zdc(1))^2);

Atest2 = @(z,zdc) [2*(z(2)-zdc(2))*z(3)*costheta(z), 2*(z(2)-zdc(2))*sintheta(z), 0, -1, 0];
btest2 = @(z,zdc) -2*z(3)^2*sintheta(z)^2 - 1*2*(z(2)-zdc(2))*z(3)*sintheta(z) - .1*(2*(z(2)-zdc(2))*z(3)*sintheta(z) + 1*(z(2)-zdc(2))^2);


% Input
lb = [-c(1); -c(3); 0; 0; 0];
ub = [c(1); c(3); inf; inf; inf];

% Objective
H = [2 0 0 0 0; 0 2 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
f = [0; 0; 2; 2; 1];
% Dynamics
zdyn = @(z,u) [z(3)*costheta(z); z(3)*sintheta(z);0;0] + [0 0 0 0 0; 0 0 0 0 0; 0 1 0 0 0; 1 0 0 0 0]*u;

options = optimoptions('quadprog','Display','off');
ll =1;
uvec1Matrix = zeros(24,length(tvec));
uvec2Matrix = zeros(24,length(tvec));
for ii = 1:length(x0MC)
    for kk = 1:length(y0MC)
        tendSQP = 1;
        tendnom = 1;
        flaggedSQP = 0;
        flaggednom = 0;
        clear zvecSQP znom unom uvecSQP psi2track psi1track psi0track 

        z01 = [x0MC(ii); y0MC(kk)];
        if (z01(1)-zdc(1))^2 + (z01(2)-zdc(2))^2 == 0 ||psi3_0([z01;0;0],xo3,yo3)<0 || psi2_0([z01;0;0],xo2,yo2)<0 || psi1_0([z01;0;0],xo1,yo1)<0
            continue 
        end
        z02 = [Zmin(3)/10 + (Zmax(3)/10-Zmin(3)/10)*rand(1,1); -3.14/2 + (3.14/2--3.14/2)*rand(1,1);];
        

        zvecSQP(1:4,1) = [z01; z02];
        uvecSQP(1:5,1) = [0; 0 ; 0; 0 ; 0];
        jj = 1;
        while (zvecSQP(1,end)-zdc(1))^2 + (zvecSQP(2,end)-zdc(2))^2 >= .1 && jj<(1/.1)*tf
            zcur = zvecSQP(:,jj);

            ASQP = [Atest2(zcur,zdc); Atest1(zcur,zdc); AV2(zcur,zdc); A1max(zcur,xo1,yo1); A2max(zcur,xo2,yo2); A3max(zcur,xo3,yo3); A4max(zcur,xo4); A5max(zcur,xo5); A6max(zcur,yo6); A7max(zcur,yo7)];
            bSQP = [btest2(zcur,zdc); btest1(zcur,zdc); bV2(zcur,zdc); b1max(zcur,xo1,yo1); b2max(zcur,xo2,yo2); b3max(zcur,xo3,yo3); b4max(zcur,xo4); b5max(zcur,xo5); b6max(zcur,yo6); b7max(zcur,yo7)];
            try 
            [un,fn] = quadprog(H,f,ASQP,bSQP,[],[],lb,ub);
            uvecSQP(:,jj+1) = un;
            zvecSQP(:,jj+1) = zcur + dt*zdyn(zcur,un);
            tendSQP = jj+1;
            jj = jj+1;
            catch
                jj = (1/.1)*tf;
                continue
            end
        end

        figure(2)
        set(gcf, 'Position', get(0, 'Screensize'));
        subplot(1,3,1)
        plot(zvecSQP(1,:),(zvecSQP(2,:)),'b')
        hold on
        obstx3 = r3.*cos(0:.1:2*pi)+xo3;
        obsty3 = r3.*sin(0:.1:2*pi)+yo3;
        plot(obstx3,obsty3,'r--','LineWidth',2)
        hold on
        
        plot(zvecSQP(1,1),zvecSQP(2,1),'gx','MarkerSize',10,'MarkerEdgeColor','k')
        hold on
        plot(zdc(1),zdc(2),'go','MarkerSize',9 ,'MarkerFaceColor','g','MarkerEdgeColor','k')
        hold on
        ylim([-5 70]);
        xlim([-15 60]);
        xlabel('Position, $x$','Interpreter','Latex')
        ylabel('Position, $y$', 'Interpreter','Latex')
        grayColor = [.7 .7 .7];
        xline(0,'Color','r','LineStyle','--','LineWidth',2)
        hold on
        xline(50,'Color','r','LineStyle','--','LineWidth',2)
        hold on
        yline(0,'Color','r','LineStyle','--','LineWidth',2)
        hold on
        yline(50,'Color','r','LineStyle','--','LineWidth',2)
        hold on
        obstx1 = r1.*cos(0:.1:2*pi)+xo1;
        obsty1 = r1.*sin(0:.1:2*pi)+yo1;
        plot(obstx1,obsty1,'r--','LineWidth',2)
        hold on
        obstx2 = r2.*cos(0:.1:2*pi)+xo2;
        obsty2 = r2.*sin(0:.1:2*pi)+yo2;
        plot(obstx2,obsty2,'r--','LineWidth',2)
        hold on
        
        uvec1Matrix(ll,1:length(uvecSQP(1,:))) = uvecSQP(1,:);
        uvec2Matrix(ll,1:length(uvecSQP(2,:))) = uvecSQP(2,:);
        ll = ll+1;

    end
end

subplot(1,3,1)
fill([0, 0, 50, 50],[50, 75, 50, 75],'red','FaceAlpha',0.2,'EdgeColor','none')
hold on
fill([0, 50, 0, 50],[50, 50, 75, 75],'red','FaceAlpha',0.2,'EdgeColor','none')
hold on
fill([0, 0, 50, 50],[0, -15, 0, -15],'red','FaceAlpha',0.2,'EdgeColor','none')
hold on
fill([0, 50, 0, 50],[0, 0, -15, -15],'red','FaceAlpha',0.2,'EdgeColor','none')
hold on
fill([50, 50, 70, 70],[-15, 80, -15, 80],'red','FaceAlpha',0.2,'EdgeColor','none')
hold on
fill([50, 70, 50, 70],[-15, -15, 80, 80],'red','FaceAlpha',0.2,'EdgeColor','none')
hold on
fill([-15, -15, 0, 0],[-15, 70, -15, 70],'red','FaceAlpha',0.2,'EdgeColor','none')
hold on
fill([-15, 0, -15, 0],[-15, -15, 70, 70],'red','FaceAlpha',0.2,'EdgeColor','none')
hold on
fill(obstx1,obsty1,'red','FaceAlpha',0.2,'EdgeColor','none')
hold on
fill(obstx2,obsty2,'red','FaceAlpha',0.2,'EdgeColor','none')
hold on
fill(obstx3,obsty3,'red','FaceAlpha',0.2,'EdgeColor','none')
hold on
lgd = legend('Trajectory','HOCBFs $j\in\{1,2,3,4,5,6,7\}$','Initial State', 'Desired State','State Constraints', 'Interpreter','Latex')
fontsize(lgd,12,'points')

uvec1Max = max(uvec1Matrix,[],1);
uvec1Min = min(uvec1Matrix,[],1);
subplot(1,3,2)
plot(tvec,uvec1Max,'Color',[0.4940 0.1840 0.5560])
hold on
plot(tvec,uvec1Min,'Color','m')
hold on
yline(lb(1),'Color',grayColor,'LineStyle','--','LineWidth',2)
shade(tvec,uvec1Max,tvec,uvec1Min,'FillColor',[0.4940 0.1840 0.5560],'FillType',[1 0; 0 2])
hold on
plot(tvec,uvec1Max,'Color',[0.4940 0.1840 0.5560])
hold on
plot(tvec,uvec1Min,'Color','m')
hold on
xlabel('Time','Interpreter','Latex')
ylabel('$u_1$','Interpreter','Latex')
ylim([-1.5, 1.5])
yline(lb(1),'Color',grayColor,'LineStyle','--','LineWidth',2)
yline(ub(1),'Color',grayColor,'LineStyle','--','LineWidth',2)
lgd = legend('Maximum Input over Trials','Minimum Input over Trials','Input Bounds','Interpreter','Latex')
fontsize(lgd,12,'points')

uvec2Max = max(uvec2Matrix,[],1);
uvec2Min = min(uvec2Matrix,[],1);
subplot(1,3,3)
plot(tvec,uvec2Max,'Color',[0.4660 0.6740 0.1880])
hold on
plot(tvec,uvec2Min,'Color',[0.9290 0.6940 0.1250])
hold on
yline(lb(2),'Color',grayColor,'LineStyle','--','LineWidth',2)
shade(tvec,uvec2Max,tvec,uvec2Min,'FillColor',[0.4660 0.6740 0.1880],'FillType',[1 0; 0 2])
hold on
plot(tvec,uvec2Max,'Color',[0.4660 0.6740 0.1880])
hold on
plot(tvec,uvec2Min,'Color',[0.9290 0.6940 0.1250])
hold on
xlabel('Time','Interpreter','Latex')
ylabel('$u_2$','Interpreter','Latex')
ylim([-2.5, 2.5])
yline(lb(2),'Color',grayColor,'LineStyle','--','LineWidth',2)
yline(ub(2),'Color',grayColor,'LineStyle','--','LineWidth',2)
lgd = legend('Maximum Input over Trials','Minimum Input over Trials','Input Bounds','Interpreter','Latex')
fontsize(lgd,12,'points')