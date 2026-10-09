# Glue (adapted from submission/fit/fit.py): batman + emcee fit for TIC 375939912 candidate.
# Input: fitdata.npz here with t (BTJD), f (detrended ~1), e, s (sector), x (exposure s), written by prep.py.
import numpy as np, emcee, batman, json, sys, os
from scipy.optimize import minimize
HERE = 'work/planet-hunt/deep/375939912/fit'
np.random.seed(42)
V = json.load(open('work/planet-hunt/hunt5b/vet_375939912.json'))
P0, T00, DUR0 = V['P'], V['T0'], V['dur_h']/24

# Stellar prior: TIC v8.2 via ExoFOP (_sources/tic_exofop.json); dens in g/cm3
st = [s for s in json.load(open(HERE+'/_sources/tic_exofop.json'))['stellar_parameters'] if s.get('snotes') == 'TIC v8.2'][0]
RHO, E_RHO = float(st['dens']), float(st['dens_e'])
TEFF, LOGG = float(st['teff']), float(st['logg'])
# Claret 2017 TESS quadratic LD, ATLAS LSM, Z=0, xi=2, bilinear in Teff (4500,4750), logg (4.5,5.0)
g = {}
for ln in open(HERE+'/_sources/claret_raw.tsv'):
    try: lg, T_, Z, xi, a_, b_ = map(float, ln.split('\t')[:6])
    except Exception: continue
    if Z == 0 and xi == 2: g[(T_, lg)] = (a_, b_)
x_ = (TEFF-4500)/250; y_ = (LOGG-4.5)/0.5
def bil(i): return (1-x_)*(1-y_)*g[(4500, 4.5)][i] + x_*(1-y_)*g[(4750, 4.5)][i] + (1-x_)*y_*g[(4500, 5.0)][i] + x_*y_*g[(4750, 5.0)][i]
U1, U2 = bil(0), bil(1); E_U = 0.1

D = np.load(HERE+'/fitdata.npz')
T, F, E, S, X = D['t'], D['f'], D['e'], D['s'], D['x']
EXPS = sorted(set(np.round(X).astype(int).tolist()))
grp = np.searchsorted(EXPS, np.round(X).astype(int))
NG = len(EXPS)
TREF_N = int(np.round((np.median(T) - T00)/P0)); TREF0 = T00 + TREF_N*P0
p = batman.TransitParams(); p.t0=0.; p.per=P0; p.rp=0.028; p.a=12.; p.inc=89.; p.ecc=0.; p.w=90.; p.limb_dark='quadratic'; p.u=[U1, U2]
MODS = [batman.TransitModel(p, T[grp == i], supersample_factor=max(1, int(np.ceil(EXPS[i]/120))), exp_time=EXPS[i]/86400.) for i in range(NG)]

def a_rs(rho, P): G = 6.674e-8; return ((G*rho*(P*86400)**2)/(3*np.pi))**(1/3)
names = ['T0', 'P', 'k', 'b', 'rho', 'q1', 'q2', 'f0'] + ['lns%d' % e for e in EXPS]
def unpack(th):
    T0, P, k, b, rho, q1, q2, f0 = th[:8]
    a = a_rs(rho, P); inc = np.degrees(np.arccos(min(b/a, 1)))
    return T0, P, k, b, rho, a, inc, 2*np.sqrt(q1)*q2, np.sqrt(q1)*(1-2*q2), f0
def model(th):
    T0, P, k, b, rho, a, inc, u1, u2, f0 = unpack(th)
    p.t0=T0; p.per=P; p.rp=k; p.a=a; p.inc=inc; p.u=[u1, u2]
    out = np.empty_like(F)
    for i in range(NG): out[grp == i] = MODS[i].light_curve(p)
    return out*f0
def lnprob(th):
    T0, P, k, b, rho, q1, q2, f0 = th[:8]; ls = np.asarray(th[8:])
    if not (0 < k < 0.2 and 0 <= b < 1+k and 0.1 < rho < 30 and 0 < q1 < 1 and 0 < q2 < 1 and np.all((ls > -14) & (ls < -4))): return -np.inf
    if abs(T0 - TREF0) > 0.1 or abs(P - P0) > 0.01: return -np.inf
    _, _, _, _, _, a, inc, u1, u2, _ = unpack(th)
    if b >= a: return -np.inf
    lp = -0.5*((rho-RHO)/E_RHO)**2 - 0.5*((u1-U1)/E_U)**2 - 0.5*((u2-U2)/E_U)**2
    var = E**2 + np.exp(2*ls)[grp]
    r = F - model(th)
    return lp - 0.5*np.sum(r**2/var + np.log(var))

if __name__ == '__main__':
    print('LD', U1, U2, 'rho', RHO, E_RHO, 'exps', EXPS, 'N', len(T), 'TREF_N', TREF_N, flush=True)
    q1_0 = (U1+U2)**2; q2_0 = U1/(2*(U1+U2))
    best = np.array([TREF0, P0, np.sqrt(V['depth_all'][0]), 0.3, RHO, q1_0, q2_0, 1.0] + [-9]*NG)
    nll = lambda x: -lnprob(x)
    for _ in range(4):
        r = minimize(nll, best, method='Nelder-Mead', options=dict(maxiter=20000, maxfev=20000, xatol=1e-9, fatol=1e-3)); best = r.x
    print('MAP', dict(zip(names, best)), -r.fun, flush=True)
    nd = len(best); nw = 4*nd
    sc = np.array([5e-4, 2e-6, 1e-3, 0.05, 0.2, 0.02, 0.02, 1e-5] + [0.05]*NG)
    pos = []
    while len(pos) < nw:
        c = best + sc*np.random.randn(nd)
        if np.isfinite(lnprob(c)): pos.append(c)
    nsteps = int(sys.argv[1]) if len(sys.argv) > 1 else 15000
    sam = emcee.EnsembleSampler(nw, nd, lnprob)
    sam.run_mcmc(np.array(pos), nsteps, progress=False)
    tau = sam.get_autocorr_time(tol=0)
    print('tau', dict(zip(names, np.round(tau, 1))))
    burn = nsteps//3
    ch = sam.get_chain(discard=burn, thin=10, flat=True)
    pc = sam.get_chain()[burn:]; W_ = pc.var(0, ddof=1).mean(0); B_ = pc.mean(0).var(0, ddof=1)*pc.shape[0]
    n_ = pc.shape[0]; rhat = np.sqrt(((n_-1)/n_*W_ + B_/n_)/W_)
    np.savez(HERE+'/chain.npz', chain=ch, tau=tau, names=names, nsteps=nsteps, burn=burn, nw=nw, rhat=rhat,
             acc=sam.acceptance_fraction, TREF_N=TREF_N, U1=U1, U2=U2, best=best)
    print('rhat', dict(zip(names, np.round(rhat, 4))), 'acc', np.mean(sam.acceptance_fraction), 'n', len(ch))
