# Reuses the hunt event-local linear-baseline estimator, with requested +/-0.5d baselines.
from pathlib import Path
import sys,json,subprocess,warnings
warnings.filterwarnings('ignore')
import numpy as np,pandas as pd
from astropy.io import fits
root=Path('work/planet-hunt/hunt5b')
tic=int(sys.argv[1]);r=json.load(open(root/f'res/{tic}.json'));P=r['P'];T0=r['T0'];secs=[int(x[0]) for x in r['used']]
star=pd.read_csv(root/'pool_all.csv').set_index('ID').loc[tic]
a=(6.674e-11*star.mass*1.989e30*(P*86400)**2/(4*np.pi**2))**(1/3)
D=max(.85*P/np.pi*np.arcsin(min(1,star.rad*6.957e8/a)),.03) # identical window to inherited autovet/vet
ts=f'{tic:016d}';path='/'.join(ts[i:i+4] for i in range(0,16,4));rows=[];issues=[]
for s in secs:
 for prod in ['spoc','qlp']:
  f=root/f'lc/{tic}/{prod}_s{s}.fits'
  if not f.exists():
   url=(f'https://archive.stsci.edu/hlsps/tess-spoc/s{s:04d}/target/{path}/hlsp_tess-spoc_tess_phot_{ts}-s{s:04d}_tess_v1_lc.fits' if prod=='spoc' else f'https://archive.stsci.edu/hlsps/qlp/s{s:04d}/{path}/hlsp_qlp_tess_ffi_s{s:04d}-{ts}_tess_v01_llc.fits')
   subprocess.run(['curl','-sf','-m','60','-o',str(f),url])
  try:
   with fits.open(f) as h:
    data=h[1].data;names=data.names;corr='PDCSAP_FLUX' if prod=='spoc' else next(c for c in ['KSPSAP_FLUX','DET_FLUX'] if c in names)
    for kind,col in [('raw','SAP_FLUX'),('corrected',corr)]:
     t=np.asarray(data['TIME'],float);v=np.asarray(data[col],float);q=np.asarray(data['QUALITY']);ok=np.isfinite(t)&np.isfinite(v)&(q==0);t=t[ok];v=v[ok];v=v/np.median(v);cad=np.median(np.diff(t))
     for epoch in range(int(np.floor((t.min()-T0)/P)),int(np.ceil((t.max()-T0)/P))+1):
      tc=T0+epoch*P;dt=t-tc;inn=np.abs(dt)<.4*D;oot=(np.abs(dt)>D)&(np.abs(dt)<.5)
      if inn.sum()<2 or (dt[oot]<0).sum()<3 or (dt[oot]>0).sum()<3:continue
      X=np.column_stack([np.ones(oot.sum()),dt[oot]]);beta=np.linalg.lstsq(X,v[oot],rcond=None)[0];res=v[oot]-X@beta;var=np.sum(res**2)/max(len(res)-2,1);xbar=np.array([1,np.mean(dt[inn])]);err=np.sqrt(var*(1/inn.sum()+xbar@np.linalg.inv(X.T@X)@xbar));depth=np.mean(beta[0]+beta[1]*dt[inn]-v[inn]);cov=((np.abs(dt)<.5*D).sum()*cad/D)
      rows.append(dict(tic=tic,sector=s,product=prod,kind=kind,column=col,epoch=epoch,tc=tc,depth_ppm=depth*1e6,err_ppm=err*1e6,coverage=cov,n_in=int(inn.sum())))
  except Exception as e:issues.append(dict(sector=s,product=prod,error=repr(e)))
e=pd.DataFrame(rows);e.to_csv(root/f'raw_events_{tic}.csv',index=False)
def combined(x):
 x=x[x.coverage>=.7];v=x.depth_ppm.values;se=x.err_ppm.values
 if not len(v):return dict(depth_ppm=None,err_ppm=None,n_events=0)
 w=1/se**2;d=np.sum(w*v)/w.sum();formal=1/np.sqrt(w.sum());scatter=np.std(v,ddof=1)/np.sqrt(len(v)) if len(v)>1 else 0
 return dict(depth_ppm=float(d),err_ppm=float(max(formal,scatter)),n_events=len(v))
out=dict(tic=tic,P=P,T0=T0,duration_d=D,baseline_halfwidth_d=.5,threshold_sigma=3,issues=issues,measurements={},gate_a={},gate_b={},stopped_at=None)
if e.empty:out['stopped_at']='a';out['gate_a']={'pass':False,'reasons':['no event measurements']}
else:
 for s in secs:
  for prod in ['spoc','qlp']:
   for kind in ['raw','corrected']:out['measurements'][f'{prod}_S{s}_{kind}']=combined(e[(e.sector==s)&(e['product']==prod)&(e.kind==kind)])
 reasons=[]
 for s in secs:
  for prod in ['spoc','qlp']:
   a=out['measurements'][f'{prod}_S{s}_raw'];b=out['measurements'][f'{prod}_S{s}_corrected'];key=f'{prod}_S{s}'
   if not a['n_events'] or not b['n_events']:reasons.append(key+': missing raw/corrected events');continue
   z=abs(a['depth_ppm']-b['depth_ppm'])/np.hypot(a['err_ppm'],b['err_ppm']);out['gate_a'][key+'_disagreement_sigma']=float(z)
   if a['depth_ppm']<=0 or b['depth_ppm']<=0:reasons.append(key+': nonpositive raw/corrected depth')
   if z>3:reasons.append(key+': raw/corrected disagreement >3sigma')
 out['gate_a'].update({'pass':not reasons,'reasons':reasons})
 if reasons:out['stopped_at']='a'
 else:
  reasons=[]
  for s in secs:
   a=out['measurements'][f'spoc_S{s}_corrected'];b=out['measurements'][f'qlp_S{s}_corrected']
   for prod,x in [('spoc',a),('qlp',b)]:
    if x['n_events']<3:reasons.append(f'{prod}_S{s}: fewer than 3 covered events')
   z=abs(a['depth_ppm']-b['depth_ppm'])/np.hypot(a['err_ppm'],b['err_ppm']);out['gate_b'][f'S{s}_disagreement_sigma']=float(z)
   if z>3:reasons.append(f'S{s}: pipeline disagreement >3sigma')
  out['gate_b'].update({'pass':not reasons,'reasons':reasons})
  if reasons:out['stopped_at']='b'
json.dump(out,open(root/f'event_gates_{tic}.json','w'),indent=2);print(json.dumps(out),flush=True)
