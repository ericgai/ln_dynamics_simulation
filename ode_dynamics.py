import numpy as np

class ln_dynamics:

    def __init__(self,e_sp=1,n_ln=30,e_ln=1.5,fac_e_i=1,fac_e_j=1,ti=12.,tj=12.,tsp=6.,tln=12.,f0=[1.,0.,0.,0.],row_labels=['spleen','LN_i','LN_j','LN_others']):
        self.e_sp=e_sp
        self.n_ln=n_ln
        self.e_ln=e_ln
        self.fac_e_i=fac_e_i
        self.fac_e_j=fac_e_j
        self.ti=ti
        self.tj=tj
        self.tsp=tsp
        self.tln=tln
        self.f=f0
        self.row_labels=row_labels

    def update_params(self):
        self.e_oln=self.e_ln*(1-2/self.n_ln)
        self.e_i=self.fac_e_i*self.e_ln*(1/self.n_ln)
        self.e_j=self.fac_e_j*self.e_ln*(1/self.n_ln)
        self.M=np.array([[-self.e_i-self.e_j-self.e_oln-self.e_sp-1/self.tsp,1/self.ti-1/self.tsp,1/self.tj-1/self.tsp,1/self.tln-1/self.tsp],
                    [self.e_i,-1/self.ti,0,0],
                    [self.e_j,0,-1/self.tj,0],
                    [self.e_oln,0,0,-1/self.tln]])
        self.k=np.array([1/self.tsp,0,0,0])
        self.lam, self.V = np.linalg.eig(self.M)
        self.w= - np.linalg.pinv(self.M) @ self.k
    
    def equilibrate(self):
        self.update_params()
        self.f=self.w

    def evolve(self,t):
        self.update_params()
        self.update_y()
        try:
            traj=np.array([(self.V @ (self.y*np.exp(self.lam*ti)) + self.w) for ti in t])
            self.f=traj[-1]
        except TypeError:
            traj=self.V @ (self.y*np.exp(self.lam*t)) + self.w
            self.f=traj
        self.update_y()
        return traj
    
    def update_y(self):
        self.y=np.linalg.pinv(self.V) @ (self.f - self.w)
