import numpy as np
import numba
import matplotlib.pyplot as plt

class ln_cells:
    rng=np.random.default_rng(0)

    def __init__(self,k1,k_in_T,diffusivity,contact_dist,r_paracortex,k_off,n_tot_Tconv_ln,n_apc,ag_specific_frequency,n_tot_Tconv,mode='vanHeijst'):
        self.n_ln_apc=0
        self.n_naive_free=0 # assume no ag specific T cells at time zero; realistic number is likely less than ten
        self.mode=mode
        self.n_active=n_tot_Tconv_ln*ag_specific_frequency  
        self.interacting=np.zeros(n_apc) # container for how many Tconv are interacting with each APC
        self.diffusivity=diffusivity
        self.contact_dist=contact_dist
        self.r_paracortex=r_paracortex
        self.n_tot_Tconv_ln=n_tot_Tconv_ln
        self.n_tot_Tconv=n_tot_Tconv
        self.n_naive_circ=np.round(n_tot_Tconv*ag_specific_frequency).astype(int)
        self.k1=k1
        if mode=='Mandl':
            self.k2=k_in_T * n_tot_Tconv_ln / n_tot_Tconv
        elif mode=='vanHeijst':
            self.k2=0.042
        else:
            raise('Unknown estimation mode for k_2')
        self.k_on=3*diffusivity*contact_dist/r_paracortex**3
        self.k_off=k_off
        self.n_tissue_apc=n_apc
        self.ag_specific_frequency=ag_specific_frequency
        self.t=0
    
    def step(self):
        r1,r2=self.rng.random(size=2) #draw random numbers
        propen_on=self.n_naive_free*self.n_ln_apc*self.k_on
        propen_off=np.sum(self.interacting)*self.k_off
        propen_apc_in=self.n_tissue_apc*self.k1
        propen_tconv_in=self.k2*self.n_naive_circ
        propen_sum=propen_on+propen_off+propen_apc_in+propen_tconv_in #calculate process propensities and their sum
        # print(propen_on,propen_off,propen_apc_in,propen_sum)
        tau=np.log(1/r1)/propen_sum
        self.t+=tau
        if r2<(propen_on/propen_sum):
            r2_norm=r2/(propen_on/propen_sum)
            apc_idx=np.floor(r2_norm*self.n_ln_apc).astype(int)
            self.interacting[apc_idx]+=1
            self.n_naive_free-=1
        elif r2<(propen_on+propen_off)/propen_sum:
            r2_norm=(r2-(propen_on/propen_sum))/(propen_off/propen_sum)
            interaction_idx=np.floor(r2_norm*np.sum(self.interacting))
            self.interacting[np.searchsorted(np.cumsum(self.interacting),interaction_idx)]-=1
            self.n_active+=1
        elif r2<(propen_on+propen_off+propen_apc_in)/propen_sum:
            self.n_tissue_apc-=1
            self.n_ln_apc+=1
        else:
            self.n_naive_circ-=1
            self.n_naive_free+=1

    def leap(self, steps):
        for _ in range(steps):
            self.step()

# time step in units of hrs
fig,ax=plt.subplots()
for mode in ['Mandl','vanHeijst']:
    simulation=ln_cells(k1=0.125,k_in_T=0.082,diffusivity=10,contact_dist=12,r_paracortex=300,k_off=0.125,n_tot_Tconv_ln=3.e5,ag_specific_frequency=5e-5,n_apc=1500,n_tot_Tconv=2.e7,mode=mode)
    n_steps=5000
    container=np.zeros(n_steps)
    t=np.zeros(n_steps)
    for i in range(n_steps):
        simulation.step()
        container[i]=np.sum(simulation.interacting)
        # container[i]=simulation.n_naive_free
        # container[i]=simulation.n_active
        t[i]=simulation.t


    ax.plot(t,container,label=mode)

ax.set_xlabel('time (hrs)')
ax.set_ylabel('number of APC-bound T$_\mathrm{conv}$')
ax.legend()

plt.show()
