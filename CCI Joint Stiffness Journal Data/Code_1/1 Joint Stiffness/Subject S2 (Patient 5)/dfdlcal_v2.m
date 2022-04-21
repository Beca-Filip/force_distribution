%% To compute the dfdl values with instaneous values
% lmt, vmt, a and muscle parameters
function dfdl_sym = dfdlcal_v2(Params)
% extract muscle parameters from input structure 'Params'
Fmax_sym = Params.Fmax_sym;
Lmt_sym = Params.Lmt_sym;
alphaAngle = Params.alphaAngle;
lmo_sym = Params.lmo_sym;
lts_sym = Params.lts_sym;
VmaxFactor_sym = Params.VmaxFactor_sym;
Vmt_sym = Params.Vmt_sym;
a_sym = Params.a_sym;

% define parameters in the models
b11 = 0.814483478343008;
b21 = 1.055033428970575;
b31 = 0.162384573599574;
b41 = 0.063303448465465;
b12 = 0.433004984392647;
b22 = 0.716775413397760;
b32 = -0.029947116970696;
b42 = 0.200356847296188;
b13 = 0.1;
b23 = 1;
b33 = 0.5.*sqrt(0.5);
b43 = 0;

d1 = -8.665833008200368;
d2 = 7.623504943764548;
d3 = 3.565407127548163;
d4 = 0.968466940650055;
d5 = -0.394453556647392;
d6 = 8.169548497993912;

e1 = 0.232000804702701;
e2 = 12.4385351684494;
e3 = 1.32947047766830;

%% sybolic expression of dfdl (derived from 'dfdlSymbolicExpression.m')
dfdl_sym = Fmax_sym.*a_sym.*cos(alphaAngle).*(d1 + d2.*atan(d3 + d4.*atan(d5 + (Vmt_sym.*d6)./(VmaxFactor_sym.*lmo_sym.*cos(alphaAngle)))))...
    .*(b11.*exp(-(b21 - (Lmt_sym - lts_sym)./(lmo_sym.*cos(alphaAngle))).^2./(2.*(b31 + (b41.*(Lmt_sym - lts_sym))./(lmo_sym.*cos(alphaAngle))).^2))...
    .*((b21 - (Lmt_sym - lts_sym)./(lmo_sym.*cos(alphaAngle)))./(lmo_sym.*cos(alphaAngle).*(b31 + (b41.*(Lmt_sym - lts_sym))./(lmo_sym.*cos(alphaAngle))).^2)...
    + (b41.*(b21 - (Lmt_sym - lts_sym)./(lmo_sym.*cos(alphaAngle))).^2)./(lmo_sym.*cos(alphaAngle).*(b31 + (b41.*(Lmt_sym - lts_sym))./(lmo_sym.*cos(alphaAngle))).^3))...
    + b12.*exp(-(b22 - (Lmt_sym - lts_sym)./(lmo_sym.*cos(alphaAngle))).^2./(2.*(b32 + (b42.*(Lmt_sym - lts_sym))./(lmo_sym.*cos(alphaAngle))).^2)).*((b22 - (Lmt_sym - lts_sym)...
    ./(lmo_sym.*cos(alphaAngle)))./(lmo_sym.*cos(alphaAngle).*(b32 + (b42.*(Lmt_sym - lts_sym))./(lmo_sym.*cos(alphaAngle))).^2) + (b42.*(b22 - (Lmt_sym - lts_sym)./(lmo_sym.*cos(alphaAngle))).^2)...
    ./(lmo_sym.*cos(alphaAngle).*(b32 + (b42.*(Lmt_sym - lts_sym))./(lmo_sym.*cos(alphaAngle))).^3)) + b13.*exp(-(b23 - (Lmt_sym - lts_sym)./(lmo_sym.*cos(alphaAngle))).^2./(2.*(b33 + (b43.*(Lmt_sym - lts_sym))...
    ./(lmo_sym.*cos(alphaAngle))).^2)).*((b23 - (Lmt_sym - lts_sym)./(lmo_sym.*cos(alphaAngle)))./(lmo_sym.*cos(alphaAngle).*(b33 + (b43.*(Lmt_sym - lts_sym))./(lmo_sym.*cos(alphaAngle))).^2) + (b43.*(b23 - (Lmt_sym - lts_sym)...
    ./(lmo_sym.*cos(alphaAngle))).^2)./(lmo_sym.*cos(alphaAngle).*(b33 + (b43.*(Lmt_sym - lts_sym))./(lmo_sym.*cos(alphaAngle))).^3))) + (Fmax_sym.*e1.*e2.*exp(-e2.*(e3 - (Lmt_sym - lts_sym)./(lmo_sym.*cos(alphaAngle)))))./(lmo_sym.*cos(alphaAngle)...
    .*(exp(-e2.*(e3 - (Lmt_sym - lts_sym)./(lmo_sym.*cos(alphaAngle)))) + 1));

dfdl_sym = dfdl_sym';

end
