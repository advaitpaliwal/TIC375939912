# Detrend exactly as deep/375939912/fit/prep.py (no +-3 T14 window cut); saves lc_fitdetrend.npz and checks against fitdata.npz.
import json, glob, os, numpy as np
from astropy.io import fits
from wotan import flatten, transit_mask
HERE=os.path.dirname(os.path.abspath(__file__))
V=json.load(open('~/Scratch/planet-hunt/hunt5b/vet_375939912.json')); P0,T00,D0=V['P'],V['T0'],V['dur_h']/24
T,F,W=[],[],[]
for fn in sorted(glob.glob('~/Scratch/planet-hunt/deep/375939912/sectors/data/mastDownload/HLSP/hlsp_tess-spoc*/*lc.fits')):
    d=fits.open(fn)[1].data
    m=(d['QUALITY']==0)&np.isfinite(d['PDCSAP_FLUX'])&np.isfinite(d['TIME'])
    t,f=d['TIME'][m],d['PDCSAP_FLUX'][m]; f=f/np.nanmedian(f)
    msk=transit_mask(time=t,period=P0,duration=2.5*D0,T0=T00)
    tr=flatten(t[~msk],f[~msk],method='biweight',window_length=0.75,return_trend=True)[1]
    g=np.isfinite(tr); f=f/np.interp(t,t[~msk][g],tr[g])
    r=f-1; sig=1.4826*np.median(np.abs(r[~msk]-np.median(r[~msk]))); w=(r<5*sig)&(r>-10*sig)
    ph=((t-T00+0.5*P0)%P0)-0.5*P0
    T.append(t[w]); F.append(f[w]); W.append(np.abs(ph[w])<3*D0)
t,f,win=map(np.concatenate,(T,F,W))
np.savez(HERE+'/lc_fitdetrend.npz',t=t,f=f)
z=np.load('~/Scratch/planet-hunt/deep/375939912/fit/fitdata.npz')
print('n',len(t),'in window',win.sum(),'fitdata n',len(z['t']),'max|df| in window',np.max(np.abs(f[win]-z['f'])) if win.sum()==len(z['t']) else 'mismatch')
