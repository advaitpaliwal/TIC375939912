import json,glob,numpy as np,pandas as pd,warnings; warnings.filterwarnings('ignore')
rows=[json.load(open(f)) for f in glob.glob('res/*.json')]
d=pd.DataFrame(rows); print('done',len(d),'errors',d.err.notna().sum() if 'err' in d else 0)
d=d[d.SDE.notna()]; d['dur_h']=d.dur*24
stars=pd.read_csv('pool_all.csv').set_index('ID')
g=pd.read_csv('gaia_flags.csv').set_index('tic') if __import__('os').path.exists('gaia_flags.csv') else None
out=[]
for _,r in d[(d.SDE>=7)&(d.ntr>=3)].iterrows():
    tic=int(r.tic); z=np.load(f'lc/{tic}.npz'); t,f,sid=z['t'],z['f'],z['s']; P,T0=r.P,r.T0
    s=stars.loc[tic]; G=6.674e-11
    a=(G*s.mass*1.989e30*(P*86400)**2/(4*np.pi**2))**(1/3); Dexp=P*24/np.pi*np.arcsin(min(1,s.rad*6.957e8/a)); D=max(0.85*Dexp/24,0.03)
    ph=((t-T0+0.5*P)%P)-0.5*P; n=np.round((t-T0)/P).astype(int)
    it=np.abs(ph)<D/2*0.8; oo=(np.abs(ph)>D)&(np.abs(ph)<3*D)
    def dep(m):
        A=f[m&it]; B=f[m&oo]
        if len(A)<2 or len(B)<5: return np.nan,np.nan
        return np.median(B)-np.mean(A), np.std(B)/np.sqrt(len(A))
    dd,de=dep(np.ones_like(f,bool)); o=dep(n%2==1); e=dep(n%2==0)
    ss=[(int(k),)+dep(sid==k) for k in np.unique(sid)]
    ss=[x for x in ss if x[1]==x[1]]
    dv=np.array([x[1:] for x in ss]); chi=float(np.sum((dv[:,0]-dd)**2/dv[:,1]**2)/max(len(dv)-1,1)) if len(dv)>1 else np.nan
    wing=(np.abs(ph)>D/2)&(np.abs(ph)<1.5*D); far=(np.abs(ph)>3*D)&(np.abs(ph)<6*D)
    wdep=(np.median(f[far])-np.mean(f[wing]))/dd if wing.any() and far.any() and dd>0 else np.nan
    ph2=((t-T0)%P)-0.5*P; sd=(np.median(f[(np.abs(ph2)>D)&(np.abs(ph2)<3*D)])-np.mean(f[np.abs(ph2)<D/2]))/de
    snrs=[x[1]/x[2] for x in ss]
    row=dict(tic=tic,SDE=r.SDE,P=P,T0=T0,dppm=dd*1e6,dsnr=dd/de,oe_sig=abs(o[0]-e[0])/np.hypot(o[1],e[1]),nsec_tr=len(ss),sec_snr_min=min(snrs) if snrs else np.nan,sec_snrs=[round(x,1) for x in snrs],dur_h=r.dur_h,dur_exp_h=Dexp,Rp=np.sqrt(max(dd,0))*s.rad*109.1,ntr=r.ntr,chi=chi,wing_frac=wdep,sec_sig=sd,Tmag=s.Tmag,Teff=s.Teff)
    if g is not None and tic in g.index: row.update(ruwe=g.loc[tic,'ruwe'],nb=g.loc[tic,'n_bright_30'])
    out.append(row)
o=pd.DataFrame(out).sort_values('SDE',ascending=False)
o['pass']=(o.SDE>9)&(o.dsnr>7)&(o.chi<3)&(o.wing_frac.abs()<0.5)&(o.sec_sig.abs()<3)&~o.P.between(12,16)&(o.Rp<6)&(o.oe_sig<3)&(o.nsec_tr>=2)&(o.sec_snr_min>1)
o['pass'] &= o.P.ge(1)
o.to_csv('autovet.csv',index=False)
pd.set_option('display.width',250)
print(o.head(30).round(2).drop(columns=['sec_snrs']).to_string())
print('PASS'); print(o[o['pass']].round(2).to_string())
