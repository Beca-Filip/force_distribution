function musc_traj_cc = muscle_traj_and_trial_cc(fref, fpred)
%MUSCLE_TRAJ_AND_TRIAL_CC computes the muscle trajectory CC for all
%muscles and all trials. Individual muscles are along dimension 1, and 
%trials along dimension 4, while samples are along dimension 3.

musc_traj_cc = zeros(size(fref, 1), size(fref, 4));
% Loop across trials and muscles
for tt =  1 : size(fref, 4)
    for mm = 1 : size(fref, 1)
        
        musc_traj_cc(mm, tt) = corr2(squeeze(fref(mm, :, :, tt)), squeeze(fpred(mm, :, :, tt)));
        
    end
end

end