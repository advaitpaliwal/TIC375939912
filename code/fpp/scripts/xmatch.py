# Map TRICERATOPS TIC stars to localization Gaia DR3 sources (hunt5b/localize/375939912/delta_chi2_sum.csv) by position.
import pickle, numpy as np, pandas as pd
S='work/planet-hunt/deep/375939912/fpp/'; LOC='work/planet-hunt/hunt5b/localize/375939912/'
tg=pickle.load(open(S+'trc/out/target_375939912_from_hunt5b.pkl','rb')); st=tg.stars
g=pd.read_csv(LOC+'gaia_field.csv',dtype={'Source':str}); d=pd.read_csv(LOC+'delta_chi2_sum.csv',dtype={'source':str})
g=g.merge(d[['source','delta_chi2_sum','min_sector_delta']],left_on='Source',right_on='source',how='left')
dec0=np.deg2rad(g.DE_ICRS); ra=g.RA_ICRS-16*g.pmRA.fillna(0)/3.6e6/np.cos(dec0); de=g.DE_ICRS-16*g.pmDE.fillna(0)/3.6e6  # 2016.0 -> 2000.0
rows=[]
for _,s in st.iterrows():
    sep=np.hypot((ra-s.ra)*np.cos(np.deg2rad(s.dec)),de-s.dec)*3600; j=int(np.argmin(sep.values)); m=sep.iloc[j]<1.5; r=g.iloc[j]
    rows.append(dict(tic=str(s.ID),Tmag=s.Tmag,trc_sep=s['sep (arcsec)'],gaia=r.Source if m else '',match_sep=round(float(sep.iloc[j]),3),
        Gmag=r.Gmag if m else np.nan,delta_chi2_sum=r.delta_chi2_sum if m else np.nan,min_sector_delta=r.min_sector_delta if m else np.nan))
o=pd.DataFrame(rows); o['excluded']=(o.delta_chi2_sum>=25)&(o.tic!='375939912')
o.to_csv(S+'exclusions.csv',index=False)
print(o[o.gaia!=''].to_string()); print('excluded',int(o.excluded.sum()),'of',len(o)-1,'neighbours; localized',len(d),'matched',int((o.gaia!='').sum())-1)
print('unmatched localized sources:',sorted(set(d.source)-set(o.gaia)))
