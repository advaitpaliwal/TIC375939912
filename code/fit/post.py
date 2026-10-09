# Post-process chain.npz -> result.json + folded model-over-data PNG (TIC 375939912)
import numpy as np, sys, json, batman
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.argv = ['x']; sys.path.insert(0, 'work/planet-hunt/deep/375939912/fit'); import fit as Fm
HERE = Fm.HERE; np.random.seed(1)
c = np.load(HERE+'/chain.npz', allow_pickle=True); ch = c['chain']; N = len(ch)
T0, P, k, b, rho, q1, q2, f0 = ch.T[:8]
a = Fm.a_rs(rho, P); inc = np.degrees(np.arccos(b/a))
u1 = 2*np.sqrt(q1)*q2; u2 = np.sqrt(q1)*(1-2*q2)
T14 = P/np.pi*np.arcsin(np.sqrt((1+k)**2-b**2)/(a*np.sin(np.radians(inc))))*24
def bp(i):
    p = batman.TransitParams(); p.t0=0; p.per=P[i]; p.rp=k[i]; p.a=a[i]; p.inc=inc[i]; p.ecc=0; p.w=90; p.limb_dark='quadratic'; p.u=[u1[i], u2[i]]; return p
sel = np.random.choice(N, min(N, 4000), replace=False)
dep = np.array([1-batman.TransitModel(bp(i), np.array([0.])).light_curve(bp(i))[0] for i in sel])*1e6
st = Fm.st
Rs = np.random.normal(float(st['srad']), float(st['srad_e']), N); Ms = np.random.normal(float(st['mass']), float(st['mass_e']), N)
Teff = np.random.normal(float(st['teff']), float(st['teff_e']), N); L = np.random.normal(float(st['lum']), float(st['lum_e']), N)
Rp = k*Rs*109.076
sma = (Ms*(P/365.25)**2)**(1/3)
Teq = Teff*np.sqrt(1/(2*a))     # Bond albedo 0, full day-night redistribution (f=1/4)
S = L/sma**2                    # S_earth, TIC luminosity + Kepler III a with TIC mass
def q(x): m, lo, hi = np.percentile(x, [50, 15.865, 84.135]); return [float(m), float(hi-m), float(m-lo)]
out = {}
for nm, x in [('P_d', P), ('T0_BJD_TDB', T0+2457000.), ('depth_ppm', dep), ('T14_h', T14), ('k_RpRs', k), ('a_Rs', a), ('b', b), ('inc_deg', inc),
              ('rho_star_gcc', rho), ('u1', u1), ('u2', u2), ('Rp_Re', Rp), ('Teq_K', Teq), ('S_Se', S), ('a_au', sma)]:
    out[nm] = q(x)
for i, nm in enumerate(c['names'][8:]): out['jitter_ppm_'+nm[3:]+'s'] = q(np.exp(ch[:, 8+i])*1e6)
tau = c['tau']; ns = int(c['nsteps']); nw = int(c['nw']); burn = int(c['burn'])
out['convergence'] = dict(nsteps=ns, nwalkers=nw, burn=burn, tau=dict(zip(c['names'].tolist(), np.round(tau, 1).tolist())),
    nsteps_over_max_tau=float(ns/np.nanmax(tau)), ess_post_burn=dict(zip(c['names'].tolist(), np.round((ns-burn)*nw/tau).tolist())),
    rhat=dict(zip(c['names'].tolist(), np.round(c['rhat'], 4).tolist())), acc_mean=float(np.mean(c['acc'])), N_thinned=N)
out['b_95pct_upper'] = float(np.percentile(b, 95)); out['corr_T0_P'] = float(np.corrcoef(T0, P)[0, 1]); out['epoch_ref_n_from_vet_T0'] = int(c['TREF_N'])
out['ld_prior'] = [float(c['U1']), float(c['U2']), 0.1]; out['rho_prior'] = [Fm.RHO, Fm.E_RHO]
json.dump(out, open(HERE+'/result.json', 'w'), indent=1)
for kx, v in out.items(): print(kx, v)

med = np.median(ch, 0); Pm, Tm = np.median(P), np.median(T0)
ph = ((Fm.T - Tm + 0.5*Pm) % Pm - 0.5*Pm)*24
fn = Fm.F/np.median(f0); res = Fm.F - Fm.model(med)
fig, ax = plt.subplots(2, 1, figsize=(7.5, 6.2), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
for A in ax:
    for s in ('top', 'right'): A.spines[s].set_visible(False)
    A.grid(axis='y', c='0.92', lw=0.6)
tt = np.linspace(-0.25, 0.25, 3000)
pm = bp(int(np.argmin(np.abs(P-Pm)))); pm.per=Pm; pm.rp=np.median(k); pm.a=np.median(a); pm.inc=np.median(inc); pm.u=[np.median(u1), np.median(u2)]
for i in np.random.choice(N, 50):
    ax[0].plot(tt*24, (batman.TransitModel(bp(i), tt).light_curve(bp(i))-1)*1e6, c='#d9822b', alpha=0.06, lw=1)
ax[0].plot(ph, (fn-1)*1e6, '.', c='0.82', ms=1.5, zorder=0, label='TESS data (all cadences)')
w = 4*np.median(T14); bins = np.arange(-w, w+1e-9, 0.25); bc = (bins[1:]+bins[:-1])/2; idx = np.digitize(ph, bins)
def bstat(y):
    m = np.array([np.mean(y[idx == j]) if (idx == j).sum() > 3 else np.nan for j in range(1, len(bins))])
    e = np.array([np.std(y[idx == j])/np.sqrt((idx == j).sum()) if (idx == j).sum() > 3 else np.nan for j in range(1, len(bins))]); return m, e
mb, eb = bstat(fn); rb, re_ = bstat(res)
ax[0].errorbar(bc, (mb-1)*1e6, eb*1e6, fmt='o', c='0.15', ms=4, lw=1, label='15-min bins')
ax[0].plot(tt*24, (batman.TransitModel(pm, tt).light_curve(pm)-1)*1e6, c='#d9822b', lw=2, label='median model (unbinned)')
ax[0].set_ylim(-1.8*np.median(dep), 0.8*np.median(dep)); ax[0].set_ylabel('relative flux [ppm]'); ax[0].legend(fontsize=7, loc='lower left', frameon=False)
ax[0].set_title('TIC 375939912 candidate   P = %.5f d   depth %.0f ppm   Rp = %.2f R$_\\oplus$' % (Pm, np.median(dep), np.median(Rp)), fontsize=10)
ax[1].errorbar(bc, rb*1e6, re_*1e6, fmt='o', c='0.15', ms=3, lw=1); ax[1].axhline(0, c='#d9822b')
ax[1].set_ylabel('resid. [ppm]'); ax[1].set_xlabel('hours from mid-transit'); ax[1].set_xlim(-w, w)
plt.tight_layout(); plt.savefig(HERE+'/folded_model_TIC375939912.png', dpi=130)
