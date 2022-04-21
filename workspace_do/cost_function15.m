function J = cost_function15(f, fpassive, f0, pcsa, m)

factivenormalizedf0 = (f-fpassive)./f0;
stressactive = (f-fpassive)./pcsa;

n = length(f);
J = sqrt(sum(0.5*m.*(factivenormalizedf0 + stressactive.^2)) / 2 / n);

end