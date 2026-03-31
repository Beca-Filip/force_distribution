function musc_traj_rmse = muscle_traj_and_trial_rmse(fref, fpred)
%MUSCLE_TRAJ_AND_TRIAL_RMSE computes the muscle trajectory RMSE for all
%muscles and all trials. Individual muscles are along dimension 1, and 
%trials along dimension 4, while samples are along dimension 3.

musc_traj_rmse = zeros(size(fref, 1), size(fref, 4));
% Loop across trials and muscles
for tt =  1 : size(fref, 4)
    for mm = 1 : size(fref, 1)
        
        musc_traj_rmse(mm, tt) = rmse(fref(mm, :, :, tt), fpred(mm, :, :, tt));
        
    end
end
end