function Jtot = cost_function(alpha, f, fmin, fmax, pcsa, vmt, M, fpassive, f0, m, r)

% J = zeros(1, 17, 'like', casadi.MX);
% 
% J(1) = cost_function1(f);
% J(2) = cost_function2(f);
% J(3) = cost_function3(f);
% J(4) = cost_function4(f);
% J(5) = cost_function5(f, fmin, fmax);
% J(6) = cost_function6(f, fmin, fmax);
% J(7) = cost_function7(f, fmin, fmax);
% J(8) = cost_function8(f, fmin, fmax);
% J(9) = cost_function9(f, pcsa);
% J(10) = cost_function10(f, pcsa);
% J(11) = cost_function11(f, pcsa);
% J(12) = cost_function12(f, pcsa);
% J(13) = cost_function13(f, vmt);
% J(14) = cost_function14(f, M);
% J(15) = cost_function15(f, fpassive, f0, pcsa, m);
% J(16) = cost_function16(f, f0, r);
% J(17) = cost_function17(f, fmax);
J = cost_function_set(f, fmin, fmax, pcsa, vmt, M, fpassive, f0, m, r);


Jtot = zeros(1, 1, 'like', casadi.MX);

for ii = 1 :  length(J)
    Jtot = Jtot + alpha(ii) * J(ii);
end

end