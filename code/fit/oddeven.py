# Odd/even box depth at fitted median ephemeris on fitdata.npz (in = |dt|<0.8*T14/2; out = |dt|>1.5*T14/2)
import json, numpy as np
H='work/planet-hunt/deep/375939912/fit/'
r=json.load(open(H+'result.json')); D=np.load(H+'fitdata.npz')
P=r['P_d'][0]; T0=r['T0_BJD_TDB'][0]-2457000; h=r['T14_h'][0]/48
n=np.round((D['t']-T0)/P).astype(int); ph=D['t']-T0-n*P; o={}
for nm,par in (('odd',1),('even',0)):
    m=(n%2==par); i=m&(np.abs(ph)<0.8*h); x=m&(np.abs(ph)>1.5*h); w=1/D['e']**2
    d=np.average(D['f'][x],weights=w[x])-np.average(D['f'][i],weights=w[i])
    s=np.std(D['f'][x])*np.sqrt(1/i.sum()+1/x.sum()); o[nm]=[d*1e6,s*1e6,int(i.sum()),len(set(n[i]))]
o['diff_sigma']=abs(o['odd'][0]-o['even'][0])/np.hypot(o['odd'][1],o['even'][1]); print(o)
json.dump(o,open(H+'oddeven.json','w'),indent=1)
