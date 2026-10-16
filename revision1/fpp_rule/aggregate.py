# Aggregate rule-based-exclusion runs (same stats as ../fpp_fitdepth/aggregate.py) + pixel-rule check of the dropped TICs.
import json, glob, os, numpy as np, pandas as pd
H=os.path.dirname(os.path.abspath(__file__))
R=[json.load(open(f)) for f in sorted(glob.glob(H+'/out/res_sel_*.json'))]
F=np.array([r['FPP'] for r in R]); N=np.array([r['NFPP'] for r in R])
P=pd.concat([pd.read_csv(H+f"/out/probs_sel_{r['run']}.csv",dtype={'ID':str}).assign(run=r['run']) for r in R])
sc=(P.groupby(['run','scenario']).prob.sum().unstack().mean()).sort_values(ascending=False)
ps=P[P.scenario.isin(['NTP','NEB','NEBx2P'])].groupby(['run','ID']).prob.sum().unstack(fill_value=0)
share=(ps.div(ps.sum(1),axis=0)).mean().sort_values(ascending=False)
# pixel rule: TIC -> Gaia DR3 matched by G mag and separation (pickle stars vs dchi2 table)
GA={'375939920':1737322487506264320,'375939917':1737322418786784256}
d=pd.read_csv('~/Scratch/planet-hunt/prerev/375939912/loc/dchi2_sum_persector.csv')
rule={t:dict(gaia=g,per_sector={int(r.sector):r.delta_chi2 for r in d[d.source==g].itertuples()},sum=float(d[d.source==g].delta_chi2.sum()),
     met=bool(d[d.source==g].delta_chi2.sum()>=25 and (d[d.source==g].delta_chi2>-25).all())) for t,g in GA.items()}
out=dict(n=len(R),seeds=[r['seed'] for r in R],FPP_mean=F.mean(),FPP_sd=F.std(ddof=1),FPP_min=F.min(),FPP_max=F.max(),
 n_le_0p035=int((F<=0.035).sum()),n_lt_0p015=int((F<0.015).sum()),NFPP_mean=N.mean(),NFPP_sd=N.std(ddof=1),NFPP_max=N.max(),
 n_stars=R[0]['n_stars'],n_before=R[0]['n_before'],dropped=R[0]['dropped'],pixel_rule=rule,tdepth_ppm=R[0]['tdepth_ppm'],boxdepth_lc_ppm=R[0]['boxdepth_lc_ppm'],
 err_ppm=R[0]['err_ppm'],nbins=R[0]['nbins'],P=R[0]['P'],T0_btjd=R[0]['T0_btjd'],scenarios=sc.to_dict(),nfpp_share_top=share.head(6).to_dict(),
 runs=[dict(run=r['run'],seed=r['seed'],FPP=r['FPP'],NFPP=r['NFPP']) for r in R])
json.dump(out,open(H+'/summary.json','w'),indent=1,default=float); print(json.dumps(out,indent=1,default=float))
