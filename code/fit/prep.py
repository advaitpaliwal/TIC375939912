# Build fitdata.npz from TESS-SPOC HLSP FFI LCs (S55, S82): PDCSAP, quality==0, wotan biweight (transit masked), +-3 T14 windows.
import numpy as np, json, glob
from astropy.io import fits
from wotan import flatten, transit_mask
HERE='work/planet-hunt/deep/375939912/fit'
V=json.load(open('work/planet-hunt/hunt5b/vet_375939912.json')); P0,T00,D0=V['P'],V['T0'],V['dur_h']/24
T,F,E,S,X=[],[],[],[],[]
for fn in sorted(glob.glob('work/planet-hunt/deep/375939912/sectors/data/mastDownload/HLSP/hlsp_tess-spoc*/*lc.fits')):
    h=fits.open(fn); d=h[1].data; s=h[0].header['SECTOR']; hd=h[1].header
    exp=round(hd.get('TIMEDEL')*86400); print(fn.split('/')[-1], s, 'TIMEDEL s', exp, 'TIMESYS', hd.get('TIMESYS'), 'TIMEREF', hd.get('TIMEREF'), 'BJDREFI', hd.get('BJDREFI'))
    m=(d['QUALITY']==0)&np.isfinite(d['PDCSAP_FLUX'])&np.isfinite(d['TIME'])
    t,f,e=d['TIME'][m],d['PDCSAP_FLUX'][m],d['PDCSAP_FLUX_ERR'][m]; e=e/np.nanmedian(f); f=f/np.nanmedian(f)
    msk=transit_mask(time=t,period=P0,duration=2.5*D0,T0=T00)
    tr=flatten(t[~msk],f[~msk],method='biweight',window_length=0.75,return_trend=True)[1]
    g=np.isfinite(tr); fl=np.interp(t,t[~msk][g],tr[g]); f,e=f/fl,e/fl
    r=f-1; oot=~msk; sig=1.4826*np.median(np.abs(r[oot]-np.median(r[oot])))
    ph=((t-T00+0.5*P0)%P0)-0.5*P0; w=(r<5*sig)&(r>-10*sig)&(np.abs(ph)<3*D0)
    T+= [t[w]];F+=[f[w]];E+=[e[w]];S+=[np.full(w.sum(),s)];X+=[np.full(w.sum(),exp)]
    print(s,'kept',w.sum(),'oot rms ppm',round(sig*1e6),'median err ppm',round(np.median(e)*1e6))
np.savez(HERE+'/fitdata.npz',t=np.concatenate(T),f=np.concatenate(F),e=np.concatenate(E),s=np.concatenate(S),x=np.concatenate(X))
