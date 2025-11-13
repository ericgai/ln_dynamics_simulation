import numpy as np
from numba.experimental import jitclass 
from numba import njit,int8,int16, int32, float64
import matplotlib.pyplot as plt
import pickle

spec = [
    ('n_ln_apc',int16),
    ('n_active',int32),
    ('n_naive_free',int32),
    ('interacting',int16[:]),
    ('n_tot_Tconv_ln',int32),
    ('n_tot_Tconv',int32),
    ('n_naive_circ',int32),
    ('k1',float64),
    ('k2',float64),
    ('k_on',float64),
    ('k_off',float64),
    ('n_tissue_apc',int16),
    ('ag_specific_frequency',float64),
    ('t',float64)
]

@jitclass(spec)
class ln_cells:
    np.random.seed(0)
    def __init__(self,k1=0.125,k_in_T=0.082,diffusivity=10.,contact_dist=12,r_paracortex=300,k_off=0.125,n_tot_Tconv_ln=3.e5,n_apc=1500,ag_specific_frequency=5.e-5,n_tot_Tconv=2.e7,recruit_fact=1,mode=1):
        self.n_ln_apc= 0
        self.n_naive_free=n_tot_Tconv_ln*ag_specific_frequency 
        self.n_active= 0  
        self.interacting=np.zeros(n_apc,dtype=int16) # container for how many Tconv are interacting with each APC
        # self.diffusivity=diffusivity
        # self.contact_dist=contact_dist
        # self.r_paracortex=r_paracortex
        self.n_tot_Tconv_ln=n_tot_Tconv_ln
        self.n_tot_Tconv=n_tot_Tconv
        self.n_naive_circ=int32(np.round(n_tot_Tconv*ag_specific_frequency))
        self.k1=k1
        if mode==1: #based on Mandl, detailed balance of lymphocyte ingress and egress
            self.k2=k_in_T * n_tot_Tconv_ln / n_tot_Tconv * recruit_fact
        elif mode==2: #based on van Heijst, arguments of complete recruitment
            self.k2=0.042
        else:
            raise KeyError('Unknown estimation mode for k_2')
        self.k_on=3*diffusivity*contact_dist/r_paracortex**3
        self.k_off=k_off
        self.n_tissue_apc=n_apc
        self.ag_specific_frequency=ag_specific_frequency
        self.t=0
    
    def step(self):
        r1,r2=np.random.rand(2) #draw random numbers
        propen_on=self.n_naive_free*self.n_ln_apc*self.k_on
        propen_off=np.sum(self.interacting)*self.k_off
        propen_apc_in=self.n_tissue_apc*self.k1
        propen_tconv_in=self.k2*self.n_naive_circ
        propen_sum=propen_on+propen_off+propen_apc_in+propen_tconv_in #calculate process propensities and their sum
        if propen_sum==0.: #all possible kinetic processes have been realized, do not continue stepping
            return None
        # print(propen_on,propen_off,propen_apc_in,propen_sum)
        tau=np.log(1/r1)/propen_sum
        self.t+=tau
        if r2<(propen_on/propen_sum):
            r2_norm=r2/(propen_on/propen_sum)
            apc_idx=int16(np.floor(r2_norm*self.n_ln_apc))
            self.interacting[apc_idx]+=1
            self.n_naive_free-=1
        elif r2<(propen_on+propen_off)/propen_sum:
            r2_norm=(r2-(propen_on/propen_sum))/(propen_off/propen_sum)
            interaction_idx=np.floor(r2_norm*np.sum(self.interacting))
            self.interacting[np.searchsorted(np.cumsum(self.interacting),interaction_idx+1)]-=1
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
    
    def step_until(self,t):
        while self.t<t:
            self.step()

# 
# time step in units of hrs
# fig,ax=plt.subplots()
# for mode in ['Mandl','vanHeijst']:
#     simulation=ln_cells(mode=mode)
#     n_steps=5000
#     container=np.zeros(n_steps)
#     t=np.zeros(n_steps)
#     for i in range(n_steps):
#         simulation.step()
#         container[i]=np.sum(simulation.interacting)
#         # container[i]=simulation.n_naive_free
#         # container[i]=simulation.n_active
#         t[i]=simulation.t

#     ax.plot(t,container,label=mode)


# simulation=ln_cells(n_apc=100,diffusivity=1000,mode='Mandl',recruit_fact=5)
# n_steps=2500
# container=np.zeros((3,n_steps))
# t=np.zeros(n_steps)
# for i in range(n_steps):
#     simulation.step()
#     container[0,i]=np.sum(simulation.interacting)
#     container[1,i]=simulation.n_naive_free
#     container[2,i]=simulation.n_active
#     t[i]=simulation.t

# import timeit
# init_stmt='''\
# simulation=ln_cells(n_apc=100,diffusivity=1000,recruit_fact=5)
# simulation.step()'''

# loop_stmt='''\
# simulation=ln_cells(n_apc=100,diffusivity=1000,recruit_fact=5)
# for _ in range(1000):
#     simulation.step()'''

# leap_stmt='''\
# simulation=ln_cells(n_apc=100,diffusivity=1000,recruit_fact=5)
# simulation.leap(1000)'''

# leap_stmt2='''\
# simulation=ln_cells(n_apc=100,diffusivity=1000,recruit_fact=5)
# simulation.leap_while(65)'''

# timeit.timeit(init_stmt,setup='from __main__ import ln_cells',number=100)

# print(timeit.timeit(leap_stmt,setup='from __main__ import ln_cells',number=100))

# print(timeit.timeit(leap_stmt2,setup='from __main__ import ln_cells',number=100))

# t_cont=np.zeros(100)
# for i in range(100):
#     simulation=ln_cells(n_apc=100,diffusivity=1000,recruit_fact=5)
#     simulation.leap(1000)
#     t_cont[i]=simulation.t

# plt.hist(t_cont)
# plt.show()


# ax.set_title('migratory APCs')

# ax.plot(t,container[1],label='naive')
# ax.plot(t,container[0],label='APC-bound')
# ax.plot(t,container[2],label='disengaged')

# ax.set_xlabel('time (hrs)')
# ax.set_ylabel('T$_\mathrm{conv}$ count')
# ax.legend()

# plt.show()

# while simulation.t<48.:
#     simulation.step()
# with open('simulation_48h_interacting.pickle','wb') as f:
#     pickle.dump(simulation.interacting,f)

# ax.hist(simulation.interacting,bins=np.arange(10))
# ax.set_xlabel('number of proximal activated CD4 T$_\mathrm{conv}$')
# ax.set_ylabel('count')

# plt.show()
