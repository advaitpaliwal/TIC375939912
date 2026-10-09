# Summarise TRICERATOPS runs for TIC 375939912: mean +/- sample SD, scenario sums, and nearby-star scenario accounting.
import json, glob, numpy as np, pandas as pd, pickle
F='work/planet-hunt/deep/375939912/fpp/trc/out/'
d=pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(F+'res_*.json'))])
d['exptime_s']=d.get('exptime_s',pd.Series(1800.0,index=d.index)).fillna(1800.0)
for (m,e),g in d.groupby(['mode','exptime_s']):
    tag=f'{m}' if e==1800 else f'{m}_exp{int(e)}'
    P=[pd.read_csv(F+f'probs_{tag}_{k}.csv') for k in g['run']]
    s=pd.concat([p.groupby('scenario')['prob'].sum() for p in P],axis=1).fillna(0).mean(axis=1)
    p0=P[0]; nb=p0[p0['ID'].astype(str)!='375939912']
    print(f"{m} exptime={int(e)}s n={len(g)} FPP {g.FPP.mean():.5f} +/- {g.FPP.std(ddof=1):.5f} [{g.FPP.min():.5f},{g.FPP.max():.5f}]  NFPP {g.NFPP.mean():.6f} +/- {g.NFPP.std(ddof=1):.6f}  n_stars {g.n_stars.iloc[0]} depth {g.depth_ppm.iloc[0]:.1f}ppm")
    print('   scenarios(mean):',' '.join(f'{k}={v:.5f}' for k,v in s.items() if v>0))
    print(f'   non-target star rows in probs: {nb.ID.nunique()} stars, {len(nb)} rows; rows with prob>0: {(nb.prob>0).sum()}; stars with nonzero NEB/NTP prob: {nb[nb.prob>0].ID.nunique()}')
tg=pickle.load(open(F+'target_375939912_from_hunt5b.pkl','rb'))
