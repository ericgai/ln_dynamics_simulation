import numpy as np
from numba import njit
import pickle
from simulation import ln_cells

@njit
def augment_sample(th_mAPC_dist,th_gfp_hi_mask,bound_dist,freq_table):
    assert th_mAPC_dist.shape[0]==len(th_gfp_hi_mask)
    cdf=np.cumsum(freq_table)
    rands=np.random.rand(th_mAPC_dist.shape[1])
    bound_th=np.searchsorted(cdf,rands)
    dist_mask=th_mAPC_dist<bound_dist
    for i in np.nonzero(bound_th)[0]:
        n_act=bound_th[i]
        prox_th_idx_arr=np.nonzero(dist_mask[:,i]*~th_gfp_hi_mask)[0]
        if len(prox_th_idx_arr)<n_act:
            n_act=len(prox_th_idx_arr)
        th_gfp_hi_mask[np.random.choice(prox_th_idx_arr,size=n_act,replace=False)]=1
    return th_gfp_hi_mask

# the function below does not work as numba does not allow calls to jitclass with arguments inside an njit function
# @njit
# def interacting_freq_table(n_runs=100,n_bins=10,k1=0.125,k_in_T=0.082,diffusivity=10.,contact_dist=12,r_paracortex=300,k_off=0.125,n_tot_Tconv_ln=3.e5,n_apc=1500,ag_specific_frequency=5.e-5,n_tot_Tconv=2.e7,recruit_fact=1,mode=1):
#     dist_cont=np.zeros((n_runs,n_bins))
#     for i in range(n_runs):
#         simulation=ln_cells(k1=k1,k_in_T=k_in_T,diffusivity=diffusivity,contact_dist=contact_dist,r_paracortex=r_paracortex,k_off=k_off,n_tot_Tconv_ln=n_tot_Tconv_ln,n_apc=n_apc,ag_specific_frequency=ag_specific_frequency,n_tot_Tconv=n_tot_Tconv,recruit_fact=recruit_fact,mode=mode)
#         simulation.step_until(48)
#         hist,_=np.histogram(simulation.interacting,bins=np.arange(n_bins+1))
#         dist_cont[i]=hist
#     agg_hist=dist_cont.sum(axis=0)
#     return agg_hist/agg_hist.sum()



# with open('240730_processed_export.pickle','rb') as f:
#     exports=pickle.load(f)

# sample='1R2'
# th_gfp_hi_mask=np.array(exports[f'{sample}_Th']['Intensity_Mean_GFP']<1.0)
# print(th_gfp_hi_mask.sum())
# th_gfp_hi_mask_new=augment_sample(exports[f'{sample}_Th_mAPC_dist'],th_gfp_hi_mask,10,np.array([0.3,0.3,0.4]))
# print(th_gfp_hi_mask_new.sum())