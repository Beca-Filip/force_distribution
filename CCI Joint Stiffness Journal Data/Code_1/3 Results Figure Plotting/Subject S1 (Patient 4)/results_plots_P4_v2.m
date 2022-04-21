% generating plots for results section
clear, clc, close all

% load data
load('P4.mat')

global fs1 fs2
fs1 = 26; % fontsize
fs2 = 24;


% figure 1, joint stiffness, and example CCI
figure_1_v1(P4)

% figure 2, r between joint stiffness and CCI
figure_2_v1(P4,'CCI_1')
figure_2_v1(P4,'CCI_2')

% table 2, CCI_1 and CCI_2 comparison
T2 = table_2(P4);

function figure_1_v1(Inpt)
% generate index and labels for antagonistic pairs for plots
CCI_type = 'CCI_1';
[CCI_1_idx.L_HipFE,~] = sort_CCI_index_v1('L_HipFE',Inpt,CCI_type);
[CCI_1_idx.L_KneeFE,~] = sort_CCI_index_v1('L_KneeFE',Inpt,CCI_type);
[CCI_1_idx.L_AnkleFE,~] = sort_CCI_index_v1('L_AnkleFE',Inpt,CCI_type);
[CCI_1_idx.R_HipFE,~] = sort_CCI_index_v1('R_HipFE',Inpt,CCI_type);
[CCI_1_idx.R_KneeFE,~] = sort_CCI_index_v1('R_KneeFE',Inpt,CCI_type);
[CCI_1_idx.R_AnkleFE,~] = sort_CCI_index_v1('R_AnkleFE',Inpt,CCI_type);
CCI_type = 'CCI_2';
[CCI_2_idx.L_HipFE,~] = sort_CCI_index_v1('L_HipFE',Inpt,CCI_type);
[CCI_2_idx.L_KneeFE,~] = sort_CCI_index_v1('L_KneeFE',Inpt,CCI_type);
[CCI_2_idx.L_AnkleFE,~] = sort_CCI_index_v1('L_AnkleFE',Inpt,CCI_type);
[CCI_2_idx.R_HipFE,~] = sort_CCI_index_v1('R_HipFE',Inpt,CCI_type);
[CCI_2_idx.R_KneeFE,~] = sort_CCI_index_v1('R_KneeFE',Inpt,CCI_type);
[CCI_2_idx.R_AnkleFE,~] = sort_CCI_index_v1('R_AnkleFE',Inpt,CCI_type);
% 6 rows * 3 columns, rows represent plot type, columns represent DoF
fig = 1;
plot_figure_1(fig,Inpt,CCI_1_idx,CCI_2_idx)
fig = 2;
plot_figure_1(fig,Inpt,CCI_1_idx,CCI_2_idx)    
end

function plot_figure_1(fig,Inpt,CCI_1_idx,CCI_2_idx)
global fs1 fs2

figure(fig)
% x axis
x_axis = [0:1:100]';
% plot the figure for either the non-paretic side or the paretic side
    switch fig
        case 1 % non-paretic side
           for i  = 1:3 % 3 rows
                switch i
                    case 1 % non-paretic side joint stiffness
                        for j = 1:3 %  number of DoFs
                            switch j
                                case 1 % HipFE
                                    K_mean = Inpt.stiffness.K_L_HipFE.mean; 
                                    K_sd = Inpt.stiffness.K_L_HipFE.sd;  
                                    ttl = 'Hip';
                                    ylb = {'Joint stiffness','(N-m/rad)'};
                                case 2 % KneeFE
                                    K_mean = Inpt.stiffness.K_L_KneeFE.mean;
                                    K_sd = Inpt.stiffness.K_L_KneeFE.sd;
                                    ttl = 'Knee';
                                case 3 % AnkleFE
                                    K_mean = Inpt.stiffness.K_L_AnkleFE.mean;
                                    K_sd = Inpt.stiffness.K_L_AnkleFE.sd;
                                    ttl = 'Ankle';
                            end
                            subplot(3,3,(i-1)*3+j)
                            hold on
                            plot(x_axis,K_mean,'b','linewidth',2)
                            fill1 = fill([x_axis;flipud(x_axis)],[K_mean-K_sd;flipud(K_mean+K_sd)],'b','linestyle','none');
                            set(fill1,'facealpha',.2) 
                            ylim([-10 230])
                            yticks([0 100 200])
                            a = get(gca,'XTickLabel');  
                            set(gca,'XTickLabel',a,'fontsize',fs2)
                            set(gca,'XTickLabelMode','auto')
                            title(ttl,'FontSize',fs1)
                            if j == 1
                                ylabel(ylb,'FontSize',fs1)
                            end
                        end
                    case 2 % non-paretic side best emg-based CCI_1
                        for j = 1:3 %  number of DoFs
                            switch j
                                case 1 % HipFE
                                    idx_best = best_in_class(Inpt,CCI_1_idx,'CCI_1','L_HipFE');
                                    switch idx_best
                                        case 3
                                            CCI = Inpt.CCI_e_delayed;
                                            lgd_txt = '(EMG delayed)';
                                        case 4
                                            CCI = Inpt.CCI_e_calibrated;
                                            lgd_txt = '(EMG calibrated)';
                                    end
                                    CCI = CCI.L_HipFE(:,:,CCI_1_idx.L_HipFE(idx_best));
                                    ylb = {'CCI_1',''};
                                    lgd = [char(Inpt.Antagon_Pairs.HipFE_label(1,CCI_1_idx.L_HipFE(idx_best))),' vs ',...
                                        char(Inpt.Antagon_Pairs.HipFE_label(2,CCI_1_idx.L_HipFE(idx_best))),lgd_txt];
                                    ylm = [0 1];
                                    ytks = 0:0.5:1;
                                case 2 % KneeFE
                                    idx_best = best_in_class(Inpt,CCI_1_idx,'CCI_1','L_KneeFE');
                                    switch idx_best
                                        case 3
                                            CCI = Inpt.CCI_e_delayed;
                                            lgd_txt = '(EMG delayed)';
                                        case 4
                                            CCI = Inpt.CCI_e_calibrated;
                                            lgd_txt = '(EMG calibrated)';
                                    end
                                    CCI = CCI.L_KneeFE(:,:,CCI_1_idx.L_KneeFE(3));
                                    lgd = [char(Inpt.Antagon_Pairs.KneeFE_label(1,CCI_1_idx.L_KneeFE(idx_best))),' vs ',...
                                        char(Inpt.Antagon_Pairs.KneeFE_label(2,CCI_1_idx.L_KneeFE(idx_best))),lgd_txt];
                                    ylm = [0 1];
                                    ytks = 0:0.5:1;
                                case 3 % AnkleFE
                                    idx_best = best_in_class(Inpt,CCI_1_idx,'CCI_1','L_AnkleFE');
                                    switch idx_best
                                        case 3
                                            CCI = Inpt.CCI_e_delayed;
                                            lgd_txt = '(EMG delayed)';
                                        case 4
                                            CCI = Inpt.CCI_e_calibrated;
                                            lgd_txt = '(EMG calibrated)';
                                    end
                                    CCI = CCI.L_AnkleFE(:,:,CCI_1_idx.L_AnkleFE(idx_best));
                                    lgd = [char(Inpt.Antagon_Pairs.AnkleFE_label(1,CCI_1_idx.L_AnkleFE(idx_best))),' vs ',...
                                        char(Inpt.Antagon_Pairs.AnkleFE_label(2,CCI_1_idx.L_AnkleFE(idx_best))),lgd_txt];
                                    ylm = [0 0.2];
                                    ytks = 0:0.1:0.2;
                            end
                            % find mean and std
                            CCI_mean = mean(CCI,2);
                            CCI_sd = std(CCI,0,2);
                            subplot(3,3,(i-1)*3+j)
                            hold on
                            plot(x_axis,CCI_mean,'b','linewidth',2)
                            fill1 = fill([x_axis;flipud(x_axis)],[CCI_mean-CCI_sd;flipud(CCI_mean+CCI_sd)],'b','linestyle','none');
                            set(fill1,'facealpha',.2)
                            a = get(gca,'XTickLabel');
                            xticks([0:20:100])
                            set(gca,'XTickLabel',a,'fontsize',fs2)
                            set(gca,'xticklabel',[])
                            ylim(ylm)
                            yticks(ytks)
                            LGD = legend(lgd);
                            LGD.FontSize = fs2;
                            if j == 1
                                ylabel(ylb,'FontSize',fs1)
                            end                    
                        end                       
                    case 3 % non-paretic side best emg-based CCI_2
                        for j = 1:3
                            switch j
                                case 1 % HipFE
                                    idx_best = best_in_class(Inpt,CCI_2_idx,'CCI_2','L_HipFE');
                                    switch idx_best
                                        case 3
                                            CCI = Inpt.CCI_2_e_delayed;
                                            lgd_txt = '(EMG delayed)';
                                        case 4
                                            CCI = Inpt.CCI_2_e_calibrated;
                                            lgd_txt = '(EMG calibrated)';
                                    end
                                    CCI = CCI.L_HipFE(:,:,CCI_2_idx.L_HipFE(idx_best));
                                    ylb = {'CCI_2',''};
                                    lgd = [char(Inpt.Antagon_Pairs.HipFE_label(1,CCI_2_idx.L_HipFE(idx_best))),' vs ',...
                                        char(Inpt.Antagon_Pairs.HipFE_label(2,CCI_2_idx.L_HipFE(idx_best))),lgd_txt];
                                    ylm = [0 1.2];
                                    ytks = 0:0.4:1.2;
                                case 2 % KneeFE
                                    idx_best = best_in_class(Inpt,CCI_2_idx,'CCI_2','L_KneeFE');
                                    switch idx_best
                                        case 3
                                            CCI = Inpt.CCI_2_e_delayed;
                                            lgd_txt = '(EMG delayed)';
                                        case 4
                                            CCI = Inpt.CCI_2_e_calibrated;
                                            lgd_txt = '(EMG calibrated)';
                                    end
                                    CCI = CCI.L_KneeFE(:,:,CCI_2_idx.L_KneeFE(idx_best));
                                    lgd = [char(Inpt.Antagon_Pairs.KneeFE_label(1,CCI_2_idx.L_KneeFE(idx_best))),' vs ',...
                                        char(Inpt.Antagon_Pairs.KneeFE_label(2,CCI_2_idx.L_KneeFE(idx_best))),lgd_txt];
                                    ylm = [0 1.2];
                                    ytks = 0:0.4:1.2;
                                case 3 % AnkleFE
                                    idx_best = best_in_class(Inpt,CCI_2_idx,'CCI_2','L_AnkleFE');
                                    switch idx_best
                                        case 3
                                            CCI = Inpt.CCI_2_e_delayed;
                                            lgd_txt = '(EMG delayed)';
                                        case 4
                                            CCI = Inpt.CCI_2_e_calibrated;
                                            lgd_txt = '(EMG calibrated)';
                                    end
                                    CCI = CCI.L_AnkleFE(:,:,CCI_2_idx.L_AnkleFE(idx_best));
                                    lgd = [char(Inpt.Antagon_Pairs.AnkleFE_label(1,CCI_2_idx.L_AnkleFE(idx_best))),' vs ',...
                                        char(Inpt.Antagon_Pairs.AnkleFE_label(2,CCI_2_idx.L_AnkleFE(idx_best))),lgd_txt];
                                    ylm = [0 1.2];
                                    ytks = 0:0.4:1.2;
                            end
                            % find mean and std
                            CCI_mean = mean(CCI,2);
                            CCI_sd = std(CCI,0,2);
                            subplot(3,3,(i-1)*3+j)
                            hold on
                            plot(x_axis,CCI_mean,'b','linewidth',2)
                            fill1 = fill([x_axis;flipud(x_axis)],[CCI_mean-CCI_sd;flipud(CCI_mean+CCI_sd)],'b','linestyle','none');
                            set(fill1,'facealpha',.2)
                            a = get(gca,'XTickLabel');
                            xticks([0:20:100])
                            set(gca,'XTickLabel',a,'fontsize',fs2)
                            set(gca,'xticklabel',[])
                            ylim(ylm)
                            yticks(ytks)
                            LGD = legend(lgd);
                            LGD.FontSize = fs2;
                            if j == 1
                                ylabel(ylb,'FontSize',fs1)
                            end
                        end   
                end
            end 
        case 2 % paretic side
            for i = 1:3 % 3 rows
                switch i
                    case 1 % paretic side joint stiffness
                        for j = 1:3 % number of DoFs
                            switch j
                                case 1 % HipFE
                                    K_mean = Inpt.stiffness.K_R_HipFE.mean;
                                    K_sd = Inpt.stiffness.K_R_HipFE.sd;
                                    ttl = 'Hip';
                                    ylb = {'Joint stiffness','(N-m/rad)'};
                                    ylm = [-10 230];
                                case 2 % KneeFE 
                                    K_mean = Inpt.stiffness.K_R_KneeFE.mean;
                                    K_sd = Inpt.stiffness.K_R_KneeFE.sd;
                                    ttl = 'Knee';
                                    ylm = [-10 230];
                                case 3 % AnkleFE
                                    K_mean = Inpt.stiffness.K_R_AnkleFE.mean;
                                    K_sd = Inpt.stiffness.K_R_AnkleFE.sd;
                                    ttl = 'Ankle';
                                    ylm = [-10 300];
                            end
                            subplot(3,3,(i-1)*3+j)
                            hold on   
                            plot(x_axis,K_mean,'r','linewidth',2)
                            fill2 = fill([x_axis;flipud(x_axis)],[K_mean-K_sd;flipud(K_mean+K_sd)],'r','linestyle','none');
                            set(fill2,'facealpha',.2)
                            ylim(ylm)
                            yticks([0 100 200 300])
                            a = get(gca,'XTickLabel');  
                            set(gca,'XTickLabel',a,'fontsize',fs2)
                            set(gca,'XTickLabelMode','auto')
                            title(ttl,'FontSize',fs1)
                            if j == 1
                                ylabel(ylb,'FontSize',fs1)
                            end
                        end
                    case 2 % paretic side best emg-based CCI_1
                        for j = 1:3
                            switch j
                                case 1 % HipFE
                                    idx_best = best_in_class(Inpt,CCI_1_idx,'CCI_1','R_HipFE');
                                    switch idx_best
                                        case 3
                                            CCI = Inpt.CCI_e_delayed;
                                            lgd_txt = '(EMG delayed)';
                                        case 4
                                            CCI = Inpt.CCI_e_calibrated;
                                            lgd_txt = '(EMG calibrated)';
                                    end
                                    CCI = CCI.R_HipFE(:,:,CCI_1_idx.R_HipFE(idx_best));
                                    ylb = {'CCI_1',''};
                                    lgd = [char(Inpt.Antagon_Pairs.HipFE_label(1,CCI_1_idx.R_HipFE(idx_best))),' vs ',...
                                        char(Inpt.Antagon_Pairs.HipFE_label(2,CCI_1_idx.R_HipFE(idx_best))),lgd_txt];
                                    ylm = [-0.05 1.2];
                                    ytks = [0:0.4:1.2];
                                case 2 % KneeFE
                                    idx_best = best_in_class(Inpt,CCI_1_idx,'CCI_1','R_KneeFE');
                                    switch idx_best
                                        case 3
                                            CCI = Inpt.CCI_e_delayed;
                                            lgd_txt = '(EMG delayed)';
                                        case 4
                                            CCI = Inpt.CCI_e_calibrated;
                                            lgd_txt = '(EMG calibrated)';
                                    end
                                    CCI = CCI.R_KneeFE(:,:,CCI_1_idx.R_KneeFE(idx_best));
                                    lgd = [char(Inpt.Antagon_Pairs.KneeFE_label(1,CCI_1_idx.R_KneeFE(idx_best))),' vs ',...
                                        char(Inpt.Antagon_Pairs.KneeFE_label(2,CCI_1_idx.R_KneeFE(idx_best))),lgd_txt];
                                    ylm = [-0.05 1.2];
                                    ytks = [0:0.4:1.2];
                                case 3 % AnkleFE
                                    idx_best = best_in_class(Inpt,CCI_1_idx,'CCI_1','R_AnkleFE');
                                    switch idx_best
                                        case 3
                                            CCI = Inpt.CCI_e_delayed;
                                            lgd_txt = '(EMG delayed)';
                                        case 4
                                            CCI = Inpt.CCI_e_calibrated;
                                            lgd_txt = '(EMG calibrated)';
                                    end
                                    CCI = CCI.R_AnkleFE(:,:,CCI_1_idx.R_AnkleFE(idx_best));
                                    lgd = [char(Inpt.Antagon_Pairs.AnkleFE_label(1,CCI_1_idx.R_AnkleFE(idx_best))),' vs ',...
                                        char(Inpt.Antagon_Pairs.AnkleFE_label(2,CCI_1_idx.R_AnkleFE(idx_best))),lgd_txt];
                                    ylm = [-0.05 1.2];
                                    ytks = [0:0.4:1.2];
                            end
                            % find mean and std
                            CCI_mean = mean(CCI,2);
                            CCI_sd = std(CCI,0,2);
                            subplot(3,3,(i-1)*3+j)
                            hold on
                            plot(x_axis,CCI_mean,'r','linewidth',2)
                            fill1 = fill([x_axis;flipud(x_axis)],[CCI_mean-CCI_sd;flipud(CCI_mean+CCI_sd)],'r','linestyle','none');
                            set(fill1,'facealpha',.2)
                            ylim(ylm)
                            yticks(ytks)
                            a = get(gca,'XTickLabel');
                            xticks([0:20:100])
                            set(gca,'XTickLabel',a,'fontsize',fs2)
                            set(gca,'xticklabel',[])
                            LGD = legend(lgd);
                            LGD.FontSize = fs2;
                            if j == 1
                                ylabel(ylb,'FontSize',fs1)
                            end                    
                        end
                    case 3 % paretic side best emg-based CCI_2
                        for j = 1:3
                            switch j
                                case 1 % HipFE
                                    idx_best = best_in_class(Inpt,CCI_2_idx,'CCI_2','R_HipFE');
                                    switch idx_best
                                        case 3
                                            CCI = Inpt.CCI_2_e_delayed;
                                            lgd_txt = '(EMG delayed)';
                                        case 4
                                            CCI = Inpt.CCI_2_e_calibrated;
                                            lgd_txt = '(EMG calibrated)';
                                    end
                                    CCI = CCI.R_HipFE(:,:,CCI_2_idx.R_HipFE(idx_best));
                                    ylb = {'CCI_2',''};
                                    lgd = [char(Inpt.Antagon_Pairs.HipFE_label(1,CCI_2_idx.R_HipFE(idx_best))),' vs ',...
                                        char(Inpt.Antagon_Pairs.HipFE_label(2,CCI_2_idx.R_HipFE(idx_best))),lgd_txt];
                                    ylm = [-0.05 1.2];
                                    ytks = [0:0.4:1.2];
                                case 2 % KneeFE
                                    idx_best = best_in_class(Inpt,CCI_2_idx,'CCI_2','R_KneeFE');
                                    switch idx_best
                                        case 3
                                            CCI = Inpt.CCI_2_e_delayed;
                                            lgd_txt = '(EMG delayed)';
                                        case 4
                                            CCI = Inpt.CCI_2_e_calibrated;
                                            lgd_txt = '(EMG calibrated)';
                                    end
                                    CCI = CCI.R_KneeFE(:,:,CCI_2_idx.R_KneeFE(idx_best));
                                    lgd = [char(Inpt.Antagon_Pairs.KneeFE_label(1,CCI_2_idx.R_KneeFE(idx_best))),' vs ',...
                                        char(Inpt.Antagon_Pairs.KneeFE_label(2,CCI_2_idx.R_KneeFE(idx_best))),lgd_txt];
                                    ylm = [-0.05 1.2];
                                    ytks = [0:0.4:1.2];
                                case 3 % AnkleFE
                                    idx_best = best_in_class(Inpt,CCI_2_idx,'CCI_2','R_AnkleFE');
                                    switch idx_best
                                        case 3
                                            CCI = Inpt.CCI_2_e_delayed;
                                            lgd_txt = '(EMG delayed)';
                                        case 4
                                            CCI = Inpt.CCI_2_e_calibrated;
                                            lgd_txt = '(EMG calibrated)';
                                    end
                                    CCI = CCI.R_AnkleFE(:,:,CCI_2_idx.R_AnkleFE(idx_best));
                                    lgd = [char(Inpt.Antagon_Pairs.AnkleFE_label(1,CCI_2_idx.R_AnkleFE(idx_best))),' vs ',...
                                        char(Inpt.Antagon_Pairs.AnkleFE_label(2,CCI_2_idx.R_AnkleFE(idx_best))),lgd_txt];
                                    ylm = [-0.05 1.2];
                                    ytks = [0:0.4:1.2];
                            end
                            % find mean and std
                            CCI_mean = mean(CCI,2);
                            CCI_sd = std(CCI,0,2);
                            subplot(3,3,(i-1)*3+j)
                            hold on
                            plot(x_axis,CCI_mean,'r','linewidth',2)
                            fill1 = fill([x_axis;flipud(x_axis)],[CCI_mean-CCI_sd;flipud(CCI_mean+CCI_sd)],'r','linestyle','none');
                            set(fill1,'facealpha',.2)
                            ylim(ylm)
                            yticks(ytks)
                            a = get(gca,'XTickLabel');
                            xticks([0:20:100])
                            set(gca,'XTickLabel',a,'fontsize',fs2)
                            set(gca,'xticklabel',[])
                            LGD = legend(lgd);
                            LGD.FontSize = fs2;
                            if j == 1
                                ylabel(ylb,'FontSize',fs1)
                            end                    
                        end
                end
            end
    end            
end

function [CCI_idx,CCI_label] = sort_CCI_index_v1(Joint,Inpt,CCI_type)
% remove antagonistic pairs with insignificant contributor(s)
% insignificant contributors contribute less than 2% of joint stiffness
    switch Joint
        case 'L_HipFE'
            % remove insignificant contributors of joint stiffness
            insignificant_agon = Inpt.Contribution_L_HipFE.agon_idx(find(Inpt.Contribution_L_HipFE.agon_amount<0.02));
            insignificant_antagon = Inpt.Contribution_L_HipFE.antagon_idx(find(Inpt.Contribution_L_HipFE.antagon_amount<0.02));
            AP = Inpt.Antagon_Pairs.HipFE_idx; % antagonistic pairs
            labels = Inpt.Antagon_Pairs.HipFE_label;
            switch CCI_type
                case 'CCI_1'
                    r_e_basic_mean = Inpt.r_e_raw.L_HipFE_mean;
                    r_e_scaled_mean = Inpt.r_e_scaled.L_HipFE_mean;
                    r_e_delayed_mean = Inpt.r_e_delayed.L_HipFE_mean;
                    r_e_calibrated_mean = Inpt.r_e_calibrated.L_HipFE_mean;
                    r_m_mus_mean = Inpt.r_m_mus.L_HipFE_mean; 
                case 'CCI_2'
                    r_e_basic_mean = Inpt.r_e_raw_2.L_HipFE_mean;
                    r_e_scaled_mean = Inpt.r_e_scaled_2.L_HipFE_mean;
                    r_e_delayed_mean = Inpt.r_e_delayed_2.L_HipFE_mean;
                    r_e_calibrated_mean = Inpt.r_e_calibrated_2.L_HipFE_mean;
                    r_m_mus_mean = Inpt.r_m_mus_2.L_HipFE_mean;                
            end          
        case 'L_KneeFE'
            insignificant_agon = Inpt.Contribution_L_KneeFE.agon_idx(find(Inpt.Contribution_L_KneeFE.agon_amount<0.02));
            insignificant_antagon = Inpt.Contribution_L_KneeFE.antagon_idx(find(Inpt.Contribution_L_KneeFE.antagon_amount<0.02));
            AP = Inpt.Antagon_Pairs.KneeFE_idx; % antagonistic pairs
            labels = Inpt.Antagon_Pairs.KneeFE_label;
            switch CCI_type
                case 'CCI_1'
                    r_e_basic_mean = Inpt.r_e_raw.L_KneeFE_mean;
                    r_e_scaled_mean = Inpt.r_e_scaled.L_KneeFE_mean;
                    r_e_delayed_mean = Inpt.r_e_delayed.L_KneeFE_mean;
                    r_e_calibrated_mean = Inpt.r_e_calibrated.L_KneeFE_mean;
                    r_m_mus_mean = Inpt.r_m_mus.L_KneeFE_mean;
                case 'CCI_2'
                    r_e_basic_mean = Inpt.r_e_raw_2.L_KneeFE_mean;
                    r_e_scaled_mean = Inpt.r_e_scaled_2.L_KneeFE_mean;
                    r_e_delayed_mean = Inpt.r_e_delayed_2.L_KneeFE_mean;
                    r_e_calibrated_mean = Inpt.r_e_calibrated_2.L_KneeFE_mean;
                    r_m_mus_mean = Inpt.r_m_mus_2.L_KneeFE_mean;                    
            end          
        case 'L_AnkleFE'
            insignificant_agon = Inpt.Contribution_L_AnkleFE.agon_idx(find(Inpt.Contribution_L_AnkleFE.agon_amount<0.02));
            insignificant_antagon = Inpt.Contribution_L_AnkleFE.antagon_idx(find(Inpt.Contribution_L_AnkleFE.antagon_amount<0.02));
            AP = Inpt.Antagon_Pairs.AnkleFE_idx; % antagonistic pairs
            labels = Inpt.Antagon_Pairs.AnkleFE_label;
            switch CCI_type
                case 'CCI_1'
                    r_e_basic_mean = Inpt.r_e_raw.L_AnkleFE_mean;
                    r_e_scaled_mean = Inpt.r_e_scaled.L_AnkleFE_mean;
                    r_e_delayed_mean = Inpt.r_e_delayed.L_AnkleFE_mean;
                    r_e_calibrated_mean = Inpt.r_e_calibrated.L_AnkleFE_mean;
                    r_m_mus_mean = Inpt.r_m_mus.L_AnkleFE_mean;
                case 'CCI_2'
                    r_e_basic_mean = Inpt.r_e_raw_2.L_AnkleFE_mean;
                    r_e_scaled_mean = Inpt.r_e_scaled_2.L_AnkleFE_mean;
                    r_e_delayed_mean = Inpt.r_e_delayed_2.L_AnkleFE_mean;
                    r_e_calibrated_mean = Inpt.r_e_calibrated_2.L_AnkleFE_mean;
                    r_m_mus_mean = Inpt.r_m_mus_2.L_AnkleFE_mean;
            end
        case 'R_HipFE'
            insignificant_agon = Inpt.Contribution_R_HipFE.agon_idx(find(Inpt.Contribution_R_HipFE.agon_amount<0.02));
            insignificant_antagon = Inpt.Contribution_R_HipFE.antagon_idx(find(Inpt.Contribution_R_HipFE.antagon_amount<0.02));
            AP = Inpt.Antagon_Pairs.HipFE_idx; % antagonistic pairs
            labels = Inpt.Antagon_Pairs.HipFE_label;
            switch CCI_type
                case 'CCI_1'
                    r_e_basic_mean = Inpt.r_e_raw.R_HipFE_mean;
                    r_e_scaled_mean = Inpt.r_e_scaled.R_HipFE_mean;
                    r_e_delayed_mean = Inpt.r_e_delayed.R_HipFE_mean;
                    r_e_calibrated_mean = Inpt.r_e_calibrated.R_HipFE_mean;
                    r_m_mus_mean = Inpt.r_m_mus.R_HipFE_mean; 
                case 'CCI_2'
                    r_e_basic_mean = Inpt.r_e_raw_2.R_HipFE_mean;
                    r_e_scaled_mean = Inpt.r_e_scaled_2.R_HipFE_mean;
                    r_e_delayed_mean = Inpt.r_e_delayed_2.R_HipFE_mean;
                    r_e_calibrated_mean = Inpt.r_e_calibrated_2.R_HipFE_mean;
                    r_m_mus_mean = Inpt.r_m_mus_2.R_HipFE_mean;
            end
        case 'R_KneeFE'
            insignificant_agon = Inpt.Contribution_R_KneeFE.agon_idx(find(Inpt.Contribution_R_KneeFE.agon_amount<0.02));
            insignificant_antagon = Inpt.Contribution_R_KneeFE.antagon_idx(find(Inpt.Contribution_R_KneeFE.antagon_amount<0.02));
            AP = Inpt.Antagon_Pairs.KneeFE_idx; % antagonistic pairs
            labels = Inpt.Antagon_Pairs.KneeFE_label;
            switch CCI_type
                case 'CCI_1'
                    r_e_basic_mean = Inpt.r_e_raw.R_KneeFE_mean;
                    r_e_scaled_mean = Inpt.r_e_scaled.R_KneeFE_mean;
                    r_e_delayed_mean = Inpt.r_e_delayed.R_KneeFE_mean;
                    r_e_calibrated_mean = Inpt.r_e_calibrated.R_KneeFE_mean;
                    r_m_mus_mean = Inpt.r_m_mus.R_KneeFE_mean;
                case 'CCI_2'
                    r_e_basic_mean = Inpt.r_e_raw_2.R_KneeFE_mean;
                    r_e_scaled_mean = Inpt.r_e_scaled_2.R_KneeFE_mean;
                    r_e_delayed_mean = Inpt.r_e_delayed_2.R_KneeFE_mean;
                    r_e_calibrated_mean = Inpt.r_e_calibrated_2.R_KneeFE_mean;
                    r_m_mus_mean = Inpt.r_m_mus_2.R_KneeFE_mean;
            end
        case 'R_AnkleFE'
            insignificant_agon = Inpt.Contribution_R_AnkleFE.agon_idx(find(Inpt.Contribution_R_AnkleFE.agon_amount<0.003));
            insignificant_antagon = Inpt.Contribution_R_AnkleFE.antagon_idx(find(Inpt.Contribution_R_AnkleFE.antagon_amount<0.02));
            AP = Inpt.Antagon_Pairs.AnkleFE_idx; % antagonistic pairs
            labels = Inpt.Antagon_Pairs.AnkleFE_label;
            switch CCI_type
                case 'CCI_1'
                    r_e_basic_mean = Inpt.r_e_raw.R_AnkleFE_mean;
                    r_e_scaled_mean = Inpt.r_e_scaled.R_AnkleFE_mean;
                    r_e_delayed_mean = Inpt.r_e_delayed.R_AnkleFE_mean;
                    r_e_calibrated_mean = Inpt.r_e_calibrated.R_AnkleFE_mean;
                    r_m_mus_mean = Inpt.r_m_mus.R_AnkleFE_mean;
                case 'CCI_2'
                    r_e_basic_mean = Inpt.r_e_raw_2.R_AnkleFE_mean;
                    r_e_scaled_mean = Inpt.r_e_scaled_2.R_AnkleFE_mean;
                    r_e_delayed_mean = Inpt.r_e_delayed_2.R_AnkleFE_mean;
                    r_e_calibrated_mean = Inpt.r_e_calibrated_2.R_AnkleFE_mean;
                    r_m_mus_mean = Inpt.r_m_mus_2.R_AnkleFE_mean;
            end
    end
    % initialize
    CCI_idx = [1:size(AP,2)];
    % remove insignificant contributors
    insignificant = [insignificant_agon,insignificant_antagon];
    col_remove_1 = [];
    for i = 1:size(insignificant,2)
        [~,col] = find(AP == insignificant(i));
        col_remove_1 = [col_remove_1;col];
    end
    col_remove_1 = unique(col_remove_1);
    
    fine_wire_remove = 1; % 0 fine-wire included
    if fine_wire_remove
        % remove the muscles whose EMG were measured by fine wire
        fine_wire_mus = [16 17 29 34 35];
        col_remove_2 = [];
        for i = 1:size(fine_wire_mus,2)
            [~,col] = find(AP == fine_wire_mus(i));
            col_remove_2 = [col_remove_2;col];
        end
        col_remove_2 = unique(col_remove_2);    
    else
        col_remove_2 = [];
    end
    col_remove = unique([col_remove_1;col_remove_2]);
    CCI_idx(:,col_remove) = [];
    % identify antagonistic pairs with highest cor between CCI and K_joint
    % EMG basic
    [~,idx1] = max(r_e_basic_mean(1,CCI_idx));
    idx_e_basic = CCI_idx(idx1);
    % EMG scaled
    [~,idx2] = max(r_e_scaled_mean(1,CCI_idx));
    idx_e_scaled = CCI_idx(idx2);
    % EMG delayed
    [~,idx3] = max(r_e_delayed_mean(1,CCI_idx));
    idx_e_delayed = CCI_idx(idx3);
    % EMG calibrated
    [~,idx4] = max(r_e_calibrated_mean(1,CCI_idx));
    idx_e_calibrated = CCI_idx(idx4);
    % muscle moment
    [~,idx5] = max(r_m_mus_mean(1,CCI_idx));
    idx_m_mus = CCI_idx(idx5);
    % put together
    CCI_idx = [idx_e_basic,idx_e_scaled,idx_e_delayed,idx_e_calibrated,idx_m_mus];
    % create pair label
    for i = 1:size(CCI_idx,2) % loop through the pair
        CCI_label{1,i} = ['\begin{tabular}{c}',labels{1,CCI_idx(i)},'\\','vs','\\',labels{2,CCI_idx(i)},'\end{tabular}'];
    end
end

function idx_best = best_in_class(Inpt,CCI_idx,CCI_type,DoF)
    switch CCI_type
        case 'CCI_1'
            r_basic = Inpt.r_e_raw;
            r_scaled = Inpt.r_e_scaled;
            r_delayed = Inpt.r_e_delayed;
            r_calibrated = Inpt.r_e_calibrated;
        case 'CCI_2'
            r_basic = Inpt.r_e_raw_2;
            r_scaled = Inpt.r_e_scaled_2;
            r_delayed = Inpt.r_e_delayed_2;
            r_calibrated = Inpt.r_e_calibrated_2;
    end
    switch DoF
        case 'L_HipFE'
            r_basic_mean = r_basic.L_HipFE_mean;
            r_scaled_mean = r_scaled.L_HipFE_mean;
            r_delayed_mean = r_delayed.L_HipFE_mean;
            r_calibrated_mean = r_calibrated.L_HipFE_mean;
            idx = CCI_idx.L_HipFE;
        case 'L_KneeFE'
            r_basic_mean = r_basic.L_KneeFE_mean;
            r_scaled_mean = r_scaled.L_KneeFE_mean;
            r_delayed_mean = r_delayed.L_KneeFE_mean;
            r_calibrated_mean = r_calibrated.L_KneeFE_mean;
            idx = CCI_idx.L_KneeFE;
        case 'L_AnkleFE'
            r_basic_mean = r_basic.L_AnkleFE_mean;
            r_scaled_mean = r_scaled.L_AnkleFE_mean;
            r_delayed_mean = r_delayed.L_AnkleFE_mean;
            r_calibrated_mean = r_calibrated.L_AnkleFE_mean;
            idx = CCI_idx.L_AnkleFE;
        case 'R_HipFE'
            r_basic_mean = r_basic.R_HipFE_mean;
            r_scaled_mean = r_scaled.R_HipFE_mean;
            r_delayed_mean = r_delayed.R_HipFE_mean;
            r_calibrated_mean = r_calibrated.R_HipFE_mean;
            idx = CCI_idx.R_HipFE;
        case 'R_KneeFE'
            r_basic_mean = r_basic.R_KneeFE_mean;
            r_scaled_mean = r_scaled.R_KneeFE_mean;
            r_delayed_mean = r_delayed.R_KneeFE_mean;
            r_calibrated_mean = r_calibrated.R_KneeFE_mean;
            idx = CCI_idx.R_KneeFE;
        case 'R_AnkleFE'
            r_basic_mean = r_basic.R_AnkleFE_mean;
            r_scaled_mean = r_scaled.R_AnkleFE_mean;
            r_delayed_mean = r_delayed.R_AnkleFE_mean;
            r_calibrated_mean = r_calibrated.R_AnkleFE_mean;
            idx = CCI_idx.R_AnkleFE;
    end
    r = [r_basic_mean(idx(1)),r_scaled_mean(idx(2)),r_delayed_mean(idx(3)),r_calibrated_mean(idx(4))];
    [~,idx_best] = max(r);
end

function figure_2_v1(Inpt,CCI_type)
figure(2)
global fs1 fs2

% number of antagonistic pairs to be plotted
nAP = 4;
% making the plots - bar plots with error bars
    switch CCI_type
        case 'CCI_1' % CCI plot
            % generate index and labels for antagonistic pairs for plots
            [CCI_idx.L_HipFE,CCI_label.L_HipFE] = sort_CCI_index_v1('L_HipFE',Inpt,CCI_type);
            [CCI_idx.L_KneeFE,CCI_label.L_KneeFE] = sort_CCI_index_v1('L_KneeFE',Inpt,CCI_type);
            [CCI_idx.L_AnkleFE,CCI_label.L_AnkleFE] = sort_CCI_index_v1('L_AnkleFE',Inpt,CCI_type);
            [CCI_idx.R_HipFE,CCI_label.R_HipFE] = sort_CCI_index_v1('R_HipFE',Inpt,CCI_type);
            [CCI_idx.R_KneeFE,CCI_label.R_KneeFE] = sort_CCI_index_v1('R_KneeFE',Inpt,CCI_type);
            [CCI_idx.R_AnkleFE,CCI_label.R_AnkleFE] = sort_CCI_index_v1('R_AnkleFE',Inpt,CCI_type);
            for i = 1:2 % two row of plots, non-paretic and paretic side
                switch i % two rows
                    case 1 % non-paretic side
                        for j = 1:3 % number of DoFs
                            switch j
                                case 1 % HipFE
                                    idx = CCI_idx.L_HipFE(1:nAP);
                                    idx_label = 1:nAP;
                                    r_e_raw = Inpt.r_e_raw.L_HipFE(:,idx);
                                    r_e_scaled = Inpt.r_e_scaled.L_HipFE(:,idx);
                                    r_e_delayed = Inpt.r_e_delayed.L_HipFE(:,idx);
                                    r_e_calibrated = Inpt.r_e_calibrated.L_HipFE(:,idx);
                                    xtk = CCI_label.L_HipFE(:,idx_label);
                                    ttl = 'Non-paretic side';
                                    ylb = {'Hip'};
                                case 2 % KneeFE
                                    idx = CCI_idx.L_KneeFE(1:nAP);
                                    idx_label = 1:nAP;
                                    r_e_raw = Inpt.r_e_raw.L_KneeFE(:,idx);
                                    r_e_scaled = Inpt.r_e_scaled.L_KneeFE(:,idx);
                                    r_e_delayed = Inpt.r_e_delayed.L_KneeFE(:,idx);
                                    r_e_calibrated = Inpt.r_e_calibrated.L_KneeFE(:,idx);
                                    xtk = CCI_label.L_KneeFE(:,idx_label);
                                    ylb = {'Knee'};
                                case 3 % AnkleFE
                                    idx = CCI_idx.L_AnkleFE(1:nAP);
                                    idx_label = 1:nAP;
                                    r_e_raw = Inpt.r_e_raw.L_AnkleFE(:,idx);
                                    r_e_scaled = Inpt.r_e_scaled.L_AnkleFE(:,idx);
                                    r_e_delayed = Inpt.r_e_delayed.L_AnkleFE(:,idx);
                                    r_e_calibrated = Inpt.r_e_calibrated.L_AnkleFE(:,idx);
                                    xtk = CCI_label.L_AnkleFE(:,idx_label);
                                    ylb = {'Ankle'};
                            end
                            subplot(3,2,(j-1)*2+i)
                            hold on
                            % x-axis
                            x_raw = [1 7 13 19];
                            x_scaled = 8;
                            x_delayed = 14;
                            x_calibrated = 20;
                            % bar plots
                            bar(x_raw,mean(r_e_raw),1/6)
                            bar(x_scaled,mean(r_e_scaled(:,2)),1)
                            bar(x_delayed,mean(r_e_delayed(:,3)),1)
                            bar(x_calibrated,mean(r_e_calibrated(:,4)),1)
                            errorbar_plot(x_raw,r_e_raw)
                            errorbar_plot(x_scaled,r_e_scaled(:,2))
                            errorbar_plot(x_delayed,r_e_delayed(:,3))
                            errorbar_plot(x_calibrated,r_e_calibrated(:,4))
                            % plot edits
                            xlim([0 21])
                            ylim([-0.5 1])
                            % xticks
                            tks = 1.5:6:19.5;
                            xticklabels({xtk{1,1},xtk{1,2},xtk{1,3},xtk{1,4}})
                            set(gca,'xtick', tks, 'XTickLabel', xtk, 'TickLabelInterpreter','latex','fontsize',fs2);
                            % title
                            if j == 1
                                title(ttl,'FontSize',fs1)
                            end
                            % ylabel
                            ylabel(ylb,'FontSize',fs1)
                            % legend
                            if j == 1
                                LGD = legend('EMG basic','EMG scaled','EMG delayed','EMG fully cal.');
                                LGD.FontSize = fs2;
                                LGD.Location = 'EastOutside';
                            end
                        end
                    case 2 % paretic side
                        for j = 1:3 % number of DoFs
                            switch j
                                case 1 % HipFE
                                    idx = CCI_idx.R_HipFE(1:nAP);
                                    idx_label = 1:nAP;
                                    r_e_raw = Inpt.r_e_raw.R_HipFE(:,idx);
                                    r_e_scaled = Inpt.r_e_scaled.R_HipFE(:,idx);
                                    r_e_delayed = Inpt.r_e_delayed.R_HipFE(:,idx);
                                    r_e_calibrated = Inpt.r_e_calibrated.R_HipFE(:,idx);
                                    xtk = CCI_label.R_HipFE(:,idx_label);
                                    ttl = 'Paretic side';
                                case 2 % KneeFE
                                    idx = CCI_idx.R_KneeFE(1:nAP);
                                    idx_label = 1:nAP;
                                    r_e_raw = Inpt.r_e_raw.R_KneeFE(:,idx);
                                    r_e_scaled = Inpt.r_e_scaled.R_KneeFE(:,idx);
                                    r_e_delayed = Inpt.r_e_delayed.R_KneeFE(:,idx);
                                    r_e_calibrated = Inpt.r_e_calibrated.R_KneeFE(:,idx);
                                    xtk = CCI_label.R_KneeFE(:,idx_label);
                                case 3 % AnkleFE
                                    idx = CCI_idx.R_AnkleFE(1:nAP);
                                    idx_label = 1:nAP;
                                    r_e_raw = Inpt.r_e_raw.R_AnkleFE(:,idx);
                                    r_e_scaled = Inpt.r_e_scaled.R_AnkleFE(:,idx);
                                    r_e_delayed = Inpt.r_e_delayed.R_AnkleFE(:,idx);
                                    r_e_calibrated = Inpt.r_e_calibrated.R_AnkleFE(:,idx);
                                    xtk = CCI_label.R_AnkleFE(:,idx_label);
                            end
                            subplot(3,2,(j-1)*2+i)
                            hold on
                            % x-axis
                            x_raw = [1 7 13 19];
                            x_scaled = 8;
                            x_delayed = 14;
                            x_calibrated = 20;
                            % bar plots
                            bar(x_raw,mean(r_e_raw),1/6)
                            bar(x_scaled,mean(r_e_scaled(:,2)),1)
                            bar(x_delayed,mean(r_e_delayed(:,3)),1)
                            bar(x_calibrated,mean(r_e_calibrated(:,4)),1)
                            errorbar_plot(x_raw,r_e_raw)
                            errorbar_plot(x_scaled,r_e_scaled(:,2))
                            errorbar_plot(x_delayed,r_e_delayed(:,3))
                            errorbar_plot(x_calibrated,r_e_calibrated(:,4))
                            % plot edits
                            xlim([0 21])
                            ylim([-0.5 1])
                            % xticks
                            tks = 1.5:6:19.5;
                            xticklabels({xtk{1,1},xtk{1,2},xtk{1,3},xtk{1,4}})
                            set(gca,'xtick', tks, 'XTickLabel', xtk, 'TickLabelInterpreter','latex','fontsize',fs2);
                            % title
                            if j == 1
                                title(ttl,'FontSize',fs1)
                            end
                        end
                end
            end
        case 'CCI_2' % CCI2 plot
            [CCI_idx.L_HipFE,CCI_label.L_HipFE] = sort_CCI_index_v1('L_HipFE',Inpt,CCI_type);
            [CCI_idx.L_KneeFE,CCI_label.L_KneeFE] = sort_CCI_index_v1('L_KneeFE',Inpt,CCI_type);
            [CCI_idx.L_AnkleFE,CCI_label.L_AnkleFE] = sort_CCI_index_v1('L_AnkleFE',Inpt,CCI_type);
            [CCI_idx.R_HipFE,CCI_label.R_HipFE] = sort_CCI_index_v1('R_HipFE',Inpt,CCI_type);
            [CCI_idx.R_KneeFE,CCI_label.R_KneeFE] = sort_CCI_index_v1('R_KneeFE',Inpt,CCI_type);
            [CCI_idx.R_AnkleFE,CCI_label.R_AnkleFE] = sort_CCI_index_v1('R_AnkleFE',Inpt,CCI_type);
            for i = 1:2 % two row of plots, non-paretic and paretic side
                switch i % two rows
                    case 1 % non-paretic side
                        for j = 1:3 % number of DoFs
                            switch j
                                case 1 % HipFE
                                    idx = CCI_idx.L_HipFE(1:nAP);
                                    idx_label = 1:nAP;
                                    r_e_raw = Inpt.r_e_raw_2.L_HipFE(:,idx);
                                    r_e_scaled = Inpt.r_e_scaled_2.L_HipFE(:,idx);
                                    r_e_delayed = Inpt.r_e_delayed_2.L_HipFE(:,idx);
                                    r_e_calibrated = Inpt.r_e_calibrated_2.L_HipFE(:,idx);
                                    xtk = CCI_label.L_HipFE(:,idx_label);
                                    ttl = {'Non-paretic side'};
                                    ylb = {'Hip'};
                                case 2 % KneeFE
                                    idx = CCI_idx.L_KneeFE(1:nAP);
                                    idx_label = 1:nAP;
                                    r_e_raw = Inpt.r_e_raw_2.L_KneeFE(:,idx);
                                    r_e_scaled = Inpt.r_e_scaled_2.L_KneeFE(:,idx);
                                    r_e_delayed = Inpt.r_e_delayed_2.L_KneeFE(:,idx);
                                    r_e_calibrated = Inpt.r_e_calibrated_2.L_KneeFE(:,idx);
                                    xtk = CCI_label.L_KneeFE(:,idx_label);
                                    ylb = {'Knee'};
                                case 3 % AnkleFE
                                    idx = CCI_idx.L_AnkleFE(1:nAP);
                                    idx_label = 1:nAP;
                                    r_e_raw = Inpt.r_e_raw_2.L_AnkleFE(:,idx);
                                    r_e_scaled = Inpt.r_e_scaled_2.L_AnkleFE(:,idx);
                                    r_e_delayed = Inpt.r_e_delayed_2.L_AnkleFE(:,idx);
                                    r_e_calibrated = Inpt.r_e_calibrated_2.L_AnkleFE(:,idx);
                                    xtk = CCI_label.L_AnkleFE(:,idx_label);
                                    ylb = {'Ankle'};
                            end
                            subplot(3,2,(j-1)*2+i)
                            hold on
                            % x-axis
                            x_raw = [1 7 13 19];
                            x_scaled = 8;
                            x_delayed = 14;
                            x_calibrated = 20;
                            % bar plots
                            bar(x_raw,mean(r_e_raw),1/6)
                            bar(x_scaled,mean(r_e_scaled(:,2)),1)
                            bar(x_delayed,mean(r_e_delayed(:,3)),1)
                            bar(x_calibrated,mean(r_e_calibrated(:,4)),1)
                            errorbar_plot(x_raw,r_e_raw)
                            errorbar_plot(x_scaled,r_e_scaled(:,2))
                            errorbar_plot(x_delayed,r_e_delayed(:,3))
                            errorbar_plot(x_calibrated,r_e_calibrated(:,4))
                            % plot edits
                            xlim([0 21])
                            ylim([-0.5 1])
                            % xticks
                            tks = 1.5:6:19.5;
                            xticklabels({xtk{1,1},xtk{1,2},xtk{1,3},xtk{1,4}})
                            set(gca,'xtick', tks, 'XTickLabel', xtk, 'TickLabelInterpreter','latex','fontsize',fs2);
                            % title
                            if j == 1
                                title(ttl,'FontSize',fs1)
                            end
                            % ylabel
                            ylabel(ylb,'FontSize',fs1)
                            % legend
                            if j == 1
                                LGD = legend('EMG basic','EMG scaled','EMG delayed','EMG fully cal.');
                                LGD.FontSize = fs2;
                                LGD.Location = 'EastOutside';
                            end
                        end
                    case 2 % paretic side
                        for j = 1:3 % number of DoFs
                            switch j
                                case 1 % HipFE
                                    idx = CCI_idx.R_HipFE(1:nAP);
                                    idx_label = 1:nAP;
                                    r_e_raw = Inpt.r_e_raw_2.R_HipFE(:,idx);
                                    r_e_scaled = Inpt.r_e_scaled_2.R_HipFE(:,idx);
                                    r_e_delayed = Inpt.r_e_delayed_2.R_HipFE(:,idx);
                                    r_e_calibrated = Inpt.r_e_calibrated_2.R_HipFE(:,idx);
                                    xtk = CCI_label.R_HipFE(:,idx_label);
                                    ttl = 'Paretic side';
                                case 2 % KneeFE
                                    idx = CCI_idx.R_KneeFE(1:nAP);
                                    idx_label = 1:nAP;
                                    r_e_raw = Inpt.r_e_raw_2.R_KneeFE(:,idx);
                                    r_e_scaled = Inpt.r_e_scaled_2.R_KneeFE(:,idx);
                                    r_e_delayed = Inpt.r_e_delayed_2.R_KneeFE(:,idx);
                                    r_e_calibrated = Inpt.r_e_calibrated_2.R_KneeFE(:,idx);
                                    xtk = CCI_label.R_KneeFE(:,idx_label);
                                case 3 % AnkleFE
                                    idx = CCI_idx.R_AnkleFE(1:nAP);
                                    idx_label = 1:nAP;
                                    r_e_raw = Inpt.r_e_raw_2.R_AnkleFE(:,idx);
                                    r_e_scaled = Inpt.r_e_scaled_2.R_AnkleFE(:,idx);
                                    r_e_delayed = Inpt.r_e_delayed_2.R_AnkleFE(:,idx);
                                    r_e_calibrated = Inpt.r_e_calibrated_2.R_AnkleFE(:,idx);
                                    xtk = CCI_label.R_AnkleFE(:,idx_label);
                            end
                            subplot(3,2,(j-1)*2+i)
                            hold on
                            % x-axis
                            x_raw = [1 7 13 19];
                            x_scaled = 8;
                            x_delayed = 14;
                            x_calibrated = 20;
                            % bar plots
                            bar(x_raw,mean(r_e_raw),1/6)
                            bar(x_scaled,mean(r_e_scaled(:,2)),1)
                            bar(x_delayed,mean(r_e_delayed(:,3)),1)
                            bar(x_calibrated,mean(r_e_calibrated(:,4)),1)
                            errorbar_plot(x_raw,r_e_raw)
                            errorbar_plot(x_scaled,r_e_scaled(:,2))
                            errorbar_plot(x_delayed,r_e_delayed(:,3))
                            errorbar_plot(x_calibrated,r_e_calibrated(:,4))
                            % plot edits
                            xlim([0 21])
                            ylim([-0.5 1])
                            % xticks
                            tks = 1.5:6:19.5;
                            xticklabels({xtk{1,1},xtk{1,2},xtk{1,3},xtk{1,4}})
                            set(gca,'xtick', tks, 'XTickLabel', xtk, 'TickLabelInterpreter','latex','fontsize',fs2);
                            % title
                            if j == 1
                                title(ttl,'FontSize',fs1)
                            end
                        end
                end
            end
    end
end

function errorbar_plot(x,data)
% initialization
erh = zeros(1,size(data,2));
erl = zeros(1,size(data,2));
% look at the sign of mean value
data_mean = mean(data);
% calculate standard deviation
data_sd = std(data);
% give the sign to sd
idx_pos = find(data_mean>0);
idx_neg = find(data_mean<0);
% assign value
erh(1,idx_pos) = data_sd(idx_pos);
erl(1,idx_neg) = data_sd(idx_neg);

er = errorbar(x,mean(data),erl,erh);
er.Color = [0 0 0];
er.LineStyle = 'none';
end

function T2 = table_2(Inpt)
global fs2
% CCI_1
CCI_type = 'CCI_1';
[CCI_1_idx.L_HipFE,~] = sort_CCI_index_v1('L_HipFE',Inpt,CCI_type);
[CCI_1_idx.L_KneeFE,~] = sort_CCI_index_v1('L_KneeFE',Inpt,CCI_type);
[CCI_1_idx.L_AnkleFE,~] = sort_CCI_index_v1('L_AnkleFE',Inpt,CCI_type);
[CCI_1_idx.R_HipFE,~] = sort_CCI_index_v1('R_HipFE',Inpt,CCI_type);
[CCI_1_idx.R_KneeFE,~] = sort_CCI_index_v1('R_KneeFE',Inpt,CCI_type);
[CCI_1_idx.R_AnkleFE,~] = sort_CCI_index_v1('R_AnkleFE',Inpt,CCI_type);
% CCI_2
CCI_type = 'CCI_2';
[CCI_2_idx.L_HipFE,~] = sort_CCI_index_v1('L_HipFE',Inpt,CCI_type);
[CCI_2_idx.L_KneeFE,~] = sort_CCI_index_v1('L_KneeFE',Inpt,CCI_type);
[CCI_2_idx.L_AnkleFE,~] = sort_CCI_index_v1('L_AnkleFE',Inpt,CCI_type);
[CCI_2_idx.R_HipFE,~] = sort_CCI_index_v1('R_HipFE',Inpt,CCI_type);
[CCI_2_idx.R_KneeFE,~] = sort_CCI_index_v1('R_KneeFE',Inpt,CCI_type);
[CCI_2_idx.R_AnkleFE,~] = sort_CCI_index_v1('R_AnkleFE',Inpt,CCI_type);
% compute mean difference, significance
% initialization
T2.mean_diff = zeros(4,6);
T2.stats = zeros(4,6);
    for i = 1:4 % types of CCI input quantity
        switch i
            case 1 % EMG basic
                r1 = Inpt.r_e_raw;
                r2 = Inpt.r_e_raw_2;
            case 2 % EMG scaled
                r1 = Inpt.r_e_scaled;
                r2 = Inpt.r_e_scaled_2;
            case 3 % EMG delayed
                r1 = Inpt.r_e_delayed;
                r2 = Inpt.r_e_delayed_2;
            case 4 % EMG calibrated
                r1 = Inpt.r_e_calibrated;
                r2 = Inpt.r_e_calibrated_2;
        end
        for j = 1:6 % DoF
            switch j
                case 1 % L_HipFE
                    r1_DoF = r1.L_HipFE;
                    r2_DoF = r2.L_HipFE;
                    CCI1_idx_DoF = CCI_1_idx.L_HipFE;
                    CCI2_idx_DoF = CCI_2_idx.L_HipFE;
                case 2 % L_KneeFE
                    r1_DoF = r1.L_KneeFE;
                    r2_DoF = r2.L_KneeFE;
                    CCI1_idx_DoF = CCI_1_idx.L_KneeFE;
                    CCI2_idx_DoF = CCI_2_idx.L_KneeFE;
                case 3 % L_AnkleFE
                    r1_DoF = r1.L_AnkleFE;
                    r2_DoF = r2.L_AnkleFE;
                    CCI1_idx_DoF = CCI_1_idx.L_AnkleFE;
                    CCI2_idx_DoF = CCI_2_idx.L_AnkleFE;
                case 4 % R_HipFE
                    r1_DoF = r1.R_HipFE;
                    r2_DoF = r2.R_HipFE;
                    CCI1_idx_DoF = CCI_1_idx.R_HipFE;
                    CCI2_idx_DoF = CCI_2_idx.R_HipFE;
                case 5 % R_KneeFE
                    r1_DoF = r1.R_KneeFE;
                    r2_DoF = r2.R_KneeFE;
                    CCI1_idx_DoF = CCI_1_idx.R_KneeFE;
                    CCI2_idx_DoF = CCI_2_idx.R_KneeFE;
                case 6 % R_AnkleFE
                    r1_DoF = r1.R_AnkleFE;
                    r2_DoF = r2.R_AnkleFE;
                    CCI1_idx_DoF = CCI_1_idx.R_AnkleFE;
                    CCI2_idx_DoF = CCI_2_idx.R_AnkleFE;
            end
            r1_DoF_best = r1_DoF(:,CCI1_idx_DoF(i));
            r2_DoF_best = r2_DoF(:,CCI2_idx_DoF(i));
            T2.mean_diff(i,j) = mean(r1_DoF_best)-mean(r2_DoF_best);
            [~,T2.stats(i,j)] = ranksum(r1_DoF_best,r2_DoF_best);    
        end
    end
    cdata = round(T2.mean_diff,2);
    xvalues = {'Hip (NP)','Knee (NP)','Ankle (NP)','Hip (P)','Knee (P)','Ankle (P)'};
    yvalues = {'EMG-basic','EMG-scaled','EMG-delayed','EMG-calibrated'};
    h = heatmap(xvalues,yvalues,cdata,'Colormap',jet);
    h.FontSize = fs2;
end
