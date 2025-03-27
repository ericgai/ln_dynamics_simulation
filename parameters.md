## disengagement from APC $k_{off}$  
Miller *et al.* 2004 J Exp Med  
T-DC stable interactions for 10-15 hrs after contact, remain localized to DC 24 hrs after signal initiation  
model as exponential decay  
dwell time distribution   
$p_{dwell}(t)=\frac{1}{\tau} e^{-t/\tau}$  
integrate to yield survival function  
$p_{survive}(t)=\int^{\infty}_t dt' p_{dwell}(t')=e^{-t/\tau}$  

assume 5% survival at 24h:  
$\tau=-\frac{24 hr}{\ln(0.05)}=8hr$   
$k_{off}=\frac{1}{\tau}=0.125 hr^{-1}$  

## T cell arrival $k_2$

Mandl *et al.* Proc. Natl. Acad. Sci. USA 2012  
mean dwell time of CD4 T<sub>conv</sub> in a peripheral lymph node is 12.2h  
thus egress rate per cell is $\frac{1}{12.\mathrm{hr}}=0.082\mathrm{hr}^{-1}$ 

thus aggregate egress rate of CD4+ T<sub>conv</sub> is $0.082N \mathrm{hr}^{-1}$
where $N$ is the total number of CD4+ T<sub>conv</sub> in the lymph node 

which should be equal to the ingress rate  

this can be treated as a mass action term involving the entire circulating T cell pool, so each circulating T cell ingresses with rate $0.082\frac{N}{N_{tot}} \mathrm{hr}^{-1}$

at steady state, fraction $f=\frac{N_{0}}{N_{tot}}$ of these cells are antigen specific, then ingress rate of antigen specific cells is $0.082fN \mathrm{hr}^{-1}$, where $N_0$ is the total number of Ag specific cells

during a productive response, the antigen specific T cells are retained in the LN, therefore depleting the circulating population, aggregate Ag specific ingress rate is then $0.082\frac{NN_{circ}}{N_{tot}} \mathrm{hr}^{-1}$ where $N_{circ}$ is the remaining circulating pool 

precursor frequency is one in 1e6 naive cells. with 50 peptides frequency becomes one in 2e4  

estimate a total of 5e5 T cells in a resting popliteal LN (1e4/18um section * 50 sections)
70% are CD4+, of which 90% are T<sub>conv</sub>, then total of 3e5 CD4+ T<sub>conv</sub>  

alternatively, consider the observation in van Heijst *et al.* *Science* 2009 that most Ag specific precursors are recruited within 72 hours. 

we have $\tau=-\frac{72 \mathrm{hr}}{\ln(0.05)}=24\mathrm{hr}$   
then $k_2=\frac{1}{24 \mathrm{hr}}=0.042 \mathrm{hr}^{-1}$

## APC arrival $k_1$  
from cleared imaging, we have about 1500 APCs migrate in during the first 24hrs, assume this is 95% of the pool of  APCs   
model this as first order depletion of labelled tissue resident APCs  
$k_{1}=-\frac{\ln(0.05)}{24 \mathrm{hr}}=0.125 \mathrm{hr}^{-1}$  

## APC-T<sub>conv</sub> encounter $k_{on}$ 
Celli *et al.* Blood 2012  
rate constant for encounters, derived by Day and Lythe, in the limit where movement space (radius $R$) is much larger than contact radius $b$   
$\alpha=\frac{3Db}{R^3}$  
quoting parameters from this work  
$D=10\mu \mathrm{m^2/s}$, $b=12\mu\mathrm{m}$  
approximate measurement from images, $R=300\mu\mathrm{m}$ for an inflamed popliteal LN


at steady state there are about 3e5 total CD4 Tconv in the pLN. about 4e4 new CD4 Tconv enter every hour. There are about ~2e7 CD4 Tconv in a mouse, then it would take ~250h or ten days for all T cells to enter the LN once on average  
for each T cell to have visited the LN with 95% probability, it would take a month!

there must be a marked increase in the entry rate of T cells during inflammation   
but the paracortex only expands 2x in volume, and the exit rate would also be drastically reduced. at equilibrium this should lead to a much more significant increase in paracortex T cell mass $N=k_{in}/k_{out}$