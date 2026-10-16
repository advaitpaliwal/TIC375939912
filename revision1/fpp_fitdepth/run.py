# TRICERATOPS+MOLUSC for TIC 375939912 at the FITTED depth, NO neighbours removed.
# Copy of prerev/375939912/fpp/run.py (mode full). Changes ONLY: (1) input light curve = TESS-SPOC PDCSAP detrended exactly as
# deep/375939912/fit/prep.py (wotan biweight 0.75 d, transits masked 2.5xD0, trend interpolated, 5/10 sigma clip) but without
# the +-3 T14 window cut; (2) folded on the fitted ephemeris (deep/375939912/fit/result.json); (3) tdepth = fitted depth 867.986 ppm.
# Same: pickle, Gaia background, MOLUSC file, 10-min bins over |phase|<0.25 d, err = std of OOT bins, filt, N=1e6, exptime=30 min,
# nsamples=20, seeds 1000+k, parallel.
import sys, json, pickle, os, time, types, glob, numpy as np, warnings; warnings.filterwarnings('ignore')
S='~/Scratch/planet-hunt/deep/375939912/fpp/'; HERE=os.path.dirname(os.path.abspath(__file__))
FIT='~/Scratch/planet-hunt/deep/375939912/fit/'
sys.path.insert(0,S+'scripts')
def lc():
    z=np.load(HERE+'/lc_fitdetrend.npz'); return z['t'],z['f']  # built by build_lc.py (wotan env)
def main():
    import gaia_patch
    import triceratops.triceratops as tr, triceratops.marginal_likelihoods as ml
    class _F:
        def __init__(s,d): s.d=d
        def __getitem__(s,k):
            if isinstance(k,str): return types.SimpleNamespace(values=s.d[k].to_numpy(copy=True))
            return _F(s.d[np.asarray(k)])
    _rc=ml.read_csv; ml.read_csv=lambda f,*a,**kw: _F(_rc(f,*a,**kw)) if 'molusc' in str(f) else _rc(f,*a,**kw)
    k=int(sys.argv[1])
    R=json.load(open(FIT+'result.json')); P=R['P_d'][0]; T0=R['T0_BJD_TDB'][0]-2457000; depth=R['depth_ppm'][0]*1e-6
    dur=0.06
    t,f=lc()
    tg=pickle.load(open(S+'trc/out/target_375939912_from_hunt5b.pkl','rb'))
    assert tg.trilegal_fname=='375939912_gaia_background.csv'; tg.trilegal_fname=S+'trc/375939912_gaia_background.csv'
    ph=t-T0-np.round((t-T0)/P)*P; w=np.abs(ph)<0.25
    edges=np.arange(-0.25,0.25+1e-9,10/1440); idx=np.digitize(ph[w],edges)
    tb=np.array([ph[w][idx==i].mean() for i in np.unique(idx)]); fb=np.array([f[w][idx==i].mean() for i in np.unique(idx)])
    oot=np.abs(tb)>dur; err=float(np.std(fb[oot])); boxdepth=float(np.median(fb[oot])-fb[np.abs(tb)<dur/2*0.8].mean())
    tg.calc_depths(tdepth=depth,all_ap_pixels=tg._aps)
    mf=S+'molusc_out/molusc_for_trc_N1e6.csv'
    np.random.seed(1000+k); t0=time.time()
    tg.calc_probs(time=tb,flux_0=fb,flux_err_0=err,P_orb=P,contrast_curve_file=None,filt='TESS',N=1000000,parallel=True,
                  exptime=30/1440,nsamples=20,verbose=0,molusc_file=mf)
    tag=f'fit_{k}'; tg.probs.to_csv(HERE+f'/out/probs_{tag}.csv',index=False)
    json.dump(dict(mode='fitdepth_full',run=k,seed=1000+k,dropped=[],n_stars=len(tg.stars),tdepth_ppm=depth*1e6,boxdepth_lc_ppm=boxdepth*1e6,
                   err_ppm=err*1e6,nbins=len(tb),P=P,T0_btjd=T0,molusc=mf,FPP=float(tg.FPP),NFPP=float(tg.NFPP),sec=time.time()-t0),open(HERE+f'/out/res_{tag}.json','w'))
    print(tag,tg.FPP,tg.NFPP,round(time.time()-t0),flush=True)
if __name__=='__main__': main()
