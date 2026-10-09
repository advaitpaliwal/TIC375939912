# Which neighbours are eligible NEB/NTP hosts in TRICERATOPS (calc_depths tdepth>0), before vs after exclusion.
import sys, pickle, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
S='work/planet-hunt/deep/375939912/fpp/'
tic,P,T0,dur=375939912,4.300336864213205,2798.8583594995516,0.06
z=np.load('work/planet-hunt/hunt5b/lc/375939912_TESS-SPOC_masked.npz'); t,f=z['t'],z['f']
ph=t-T0-np.round((t-T0)/P)*P; w=np.abs(ph)<0.25; edges=np.arange(-0.25,0.25+1e-9,10/1440); idx=np.digitize(ph[w],edges)
tb=np.array([ph[w][idx==i].mean() for i in np.unique(idx)]); fb=np.array([f[w][idx==i].mean() for i in np.unique(idx)])
oot=np.abs(tb)>dur; depth=float(np.median(fb[oot])-fb[np.abs(tb)<dur/2*0.8].mean())
tg=pickle.load(open(S+'trc/out/target_375939912_from_hunt5b.pkl','rb')); tg.calc_depths(tdepth=depth,all_ap_pixels=tg._aps)
x=pd.read_csv(S+'exclusions.csv',dtype={'tic':str,'gaia':str})
st=tg.stars[['ID','Tmag','sep (arcsec)','fluxratio','tdepth']].copy(); st['ID']=st.ID.astype(str)
st=st.merge(x[['tic','gaia','delta_chi2_sum','excluded']],left_on='ID',right_on='tic',how='left').drop(columns='tic')
nb=st[st.ID!=str(tic)]
print(f'depth {depth*1e6:.1f} ppm; non-target catalogue stars {len(nb)}; flux in aperture (fluxratio>0) {int((nb.fluxratio>0).sum())}; eligible NEB/NTP hosts (tdepth>0) {int((nb.tdepth>0).sum())}')
print(f'after exclusion: stars retained {int((~nb.excluded.fillna(False)).sum())}; eligible retained {int(((nb.tdepth>0)&~nb.excluded.fillna(False)).sum())}')
print('eligible hosts:'); print(nb[nb.tdepth>0].to_string(index=False))
print('max fluxratio among retained non-eligible: %.3e'%nb[~nb.excluded.fillna(False)].fluxratio.max())
st.to_csv(S+'eligibility.csv',index=False)
