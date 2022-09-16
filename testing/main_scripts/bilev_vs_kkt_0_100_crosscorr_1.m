close all
clc

mm = 15;
tt = doc_sample_list(26);


x = data.f(mm, :, tt, doc_trial_list, doc_speed_list, doc_leg_list);
y = data_pred_bilev.f(mm, :, tt, doc_trial_list, doc_speed_list, doc_leg_list);
x = x(:)
y = y(:)

pcc = pearson_correlation_coefficient(x, y)

figure
hold on
plot(y);
plot(x);
xlabel('Realisation');
ylabel('Value');
legend('y', 'x');

figure
scatter(x, y, 10, 'x');