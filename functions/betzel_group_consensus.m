HC = importdata("GOC_allsubs_HC_commit_conmat_array.mat") ;
dist = importdata("GOC_dist.mat") ;
hemi= importdata("GOC_hemiid.mat");

[G,Gc] = fcn_group_bins(HC,dist,hemi,41)

dlmwrite('GOC_allsubs_HC_bin_distance_based_consensus.txt',G)
dlmwrite('GOC_allsubs_HC_bin_traditional_consistency_based_consensus.txt',Gc)

