# TRICERATOPS FPP for TIC 375939912; copied from hunt5b/trc/run_fpp.py (same LC, binning, depth, aperture pickle)
# + neighbour removal from tic240488148-localize/trc_rerun.py + molusc_file branch patch from tic240488148-companions/trc/run_fpp.py.
# usage (cwd = fpp/trc): run_fpp.py MODE RUN ; MODE in {base, loc, loc_molusc}
import sys, json, pickle, os, time, types, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
S='work/planet-hunt/deep/375939912/fpp/'; sys.path.insert(0,S+'scripts')
def main():
    import gaia_patch
    import triceratops.triceratops as tr, triceratops.marginal_likelihoods as ml
    class _F:
        def __init__(s,d): s.d=d
        def __getitem__(s,k):
            if isinstance(k,str): return types.SimpleNamespace(values=s.d[k].to_numpy(copy=True))
            return _F(s.d[np.asarray(k)])
    _rc=ml.read_csv; ml.read_csv=lambda f,*a,**kw: _F(_rc(f,*a,**kw)) if 'molusc' in str(f) else _rc(f,*a,**kw)
    mode,k=sys.argv[1],int(sys.argv[2]); EXP=float(os.environ.get('EXPTIME_S','1800'))
    tic,P,T0,dur=375939912,4.300336864213205,2798.8583594995516,0.06
    L='work/planet-hunt/hunt5b/lc/'
    z=np.load(L+f'{tic}_TESS-SPOC_masked.npz'); t,f=z['t'],z['f']
    tg=pickle.load(open('out/target_375939912_from_hunt5b.pkl','rb'))
    drop=[]
    if mode in ('loc','loc_molusc'):
        drop=[str(x) for x in pd.read_csv(S+'exclusions.csv',dtype={'tic':str}).query('excluded')['tic']]
        keep=~tg.stars['ID'].astype(str).isin(drop).values
        tg.pix_coords=[pc[keep] for pc in tg.pix_coords]; tg.remove_star(drop); tg.stars=tg.stars.reset_index(drop=True)
    ph=t-T0-np.round((t-T0)/P)*P; w=np.abs(ph)<0.25
    edges=np.arange(-0.25,0.25+1e-9,10/1440); idx=np.digitize(ph[w],edges)
    tb=np.array([ph[w][idx==i].mean() for i in np.unique(idx)]); fb=np.array([f[w][idx==i].mean() for i in np.unique(idx)])
    oot=np.abs(tb)>dur; err=float(np.std(fb[oot])); depth=float(np.median(fb[oot])-fb[np.abs(tb)<dur/2*0.8].mean())
    tg.calc_depths(tdepth=depth,all_ap_pixels=tg._aps)
    mf=S+'molusc_out/molusc_for_trc_N1e6.csv' if mode=='loc_molusc' else None
    np.random.seed(1000+k); t0=time.time()
    tg.calc_probs(time=tb,flux_0=fb,flux_err_0=err,P_orb=P,contrast_curve_file=None,filt='TESS',N=1000000,parallel=True,
                  exptime=EXP/86400,nsamples=20,verbose=0,molusc_file=mf)
    tag=f'{mode}_{k}' if EXP==1800 else f'{mode}_exp{int(EXP)}_{k}'; tg.probs.to_csv(f'out/probs_{tag}.csv',index=False)
    json.dump(dict(mode=mode,run=k,seed=1000+k,dropped=drop,n_stars=len(tg.stars),depth_ppm=depth*1e6,err_ppm=err*1e6,nbins=len(tb),molusc=mf,
                   exptime_s=EXP,FPP=float(tg.FPP),NFPP=float(tg.NFPP),sec=time.time()-t0),open(f'out/res_{tag}.json','w'))
    print(tag,tg.FPP,tg.NFPP,round(time.time()-t0),flush=True)
if __name__=='__main__': main()
