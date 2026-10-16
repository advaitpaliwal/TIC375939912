# TIC 375939912.01: a 2.06 Earth-radius, 4.30-day transiting planet candidate

Discovered by Advait Paliwal on 2026-10-09 in public TESS full-frame-image photometry.

## Summary

TIC 375939912 (TYC 525-1773-1) is a K dwarf at 84.5 pc. TESS full-frame-image photometry from Sectors 55 and 82 shows a transit-like signal with period 4.3003975 d, depth 868 ppm and duration 2.43 h, implying a radius of 2.06 R⊕. The signal appears in raw and corrected flux, with consistent depths in TESS-SPOC and QLP. Difference images in both sectors favour the target over each of the 23 Gaia sources within 60″. As of 2026-10-09, the star is not a TOI or CTOI.

With MOLUSC Gaia companion limits and no neighbours removed, TRICERATOPS gives a false-positive probability of **0.034 ± 0.002** and NFPP 0.028 at the fitted depth (Revision 1, below). This is above the 0.015 validation threshold. No high-resolution imaging exists. **Status: planet candidate.**

| Parameter | Value |
|---|---|
| Host | TIC 375939912 / TYC 525-1773-1 |
| Distance | 84.5 pc |
| Effective temperature | 4635 ± 121 K |
| Stellar radius | 0.708 ± 0.052 R☉ |
| Period | 4.3003975 ± 0.0000386 d |
| Reference mid-transit | BJD_TDB 2460542.81698 ± 0.0046 |
| Transit depth | 868 ± 76 ppm |
| Transit duration | 2.43 h |
| Candidate radius | 2.06 (+0.19 / −0.18) R⊕ |
| Equilibrium temperature | about 897 K |
| False-positive probability | 0.034 ± 0.002, NFPP 0.028 (no neighbours removed, fitted depth; Revision 1) |

The symmetric period and epoch errors above use the larger side of the fitted interval. Full asymmetric intervals and convergence diagnostics are in `data/fit/result.json`. The false-positive probability is conditional on the adopted model and constraints; its quoted scatter is across ten runs. The version 1.0 false-positive runs in `data/fpp/` used the search depth of 762.3 ppm; Revision 1 replaces them with runs at the fitted 868 ppm.

## Revision 1 (2026-10-16)

The research note was revised in response to the journal editor. The corrected numbers are:

- **False-positive probability, no neighbours removed:** FPP = 0.0337 ± 0.0018 (ten runs, range 0.0308–0.0368) and NFPP = 0.0281 ± 0.0008, with TRICERATOPS given the fitted depth of 868 ppm and the wotan-detrended light curve folded on the fitted ephemeris (`revision1/fpp_fitdepth/summary.json`). All 220 neighbours are kept. The version 1.0 value of 0.0075 ± 0.0019 depended on removing neighbours and on the shallower search depth, and is withdrawn.
- **NFPP carriers:** TIC 375939920 (59%, 43.8″) and TIC 375939917 (41%, 35.2″). The difference images disfavour them by summed Δχ² of 60.5 and 58.3, and by at least 21.0 in each sector.
- **Removing only those two** (following Giacalone et al. 2021): FPP = 0.0058 ± 0.0017 (range 0.0033–0.0086), NFPP = 7.4 × 10⁻⁶ (`revision1/fpp_rule/summary.json`). This rests on TESS pixels, not on ground-based photometry, so no validation is claimed.
- **MOLUSC:** only its Gaia RUWE and contrast tests were used (no imaging); 200,000 companions generated, 52% kept.
- **Figure 1** now has three panels: folded, unfolded Sector 82 with the trend, and four single transits (`revision1/figure1_combined_TIC375939912.png`).

Every change, with its source, is listed in [`revision1/CHANGES.md`](revision1/CHANGES.md). The revised note text is with the journal and is not included here.

## Contents

- [Research note](note/TIC375939912_note.md): methods, results, limitations and references.
- [Figure 1](note/figure1_TIC375939912.png): folded TESS-SPOC photometry, transit model and residuals.
- `revision1/`: Revision 1 figure, the new false-positive runs (scripts, per-run outputs and summaries) and `CHANGES.md`.
- [Ephemeris](ephemeris.csv): 104 predicted mid-transits from 2026-10-10 through 2027-12-31, with BJD_TDB, geocentric UTC and 1σ uncertainty.
- `code/fit/`: transit fit, preparation, posterior analysis and odd/even scripts.
- `code/fpp/`: companion-limit and false-positive scripts, without third-party packages.
- `code/search/`: search and event-vetting scripts.
- `data/fit/`: posterior chain, prepared fit data, fit results and folded plot.
- `data/fpp/`: TRICERATOPS and MOLUSC summaries, exclusions and eligibility tables.
- `data/sectors/`: sector results and coverage audit.
- `data/localize/`: neighbour tables, difference images and localization gate.
- `data/search/`: event-vetting output.

## Predicted transit times

Each ephemeris row propagates all 80,000 posterior samples as `T0 + epoch * P`. The reported BJD_TDB is the median, and `sigma_d` / `sigma_min` are the standard deviation in days / minutes. This retains the fitted T0–P correlation (0.7335), rather than treating the errors as independent. Timing uncertainty grows from about 15.1 to 20.6 minutes over the listed interval.

UTC times are geocentric: the barycentric light-travel delay toward the target is removed before converting TDB to UTC. They are not observatory-specific times. The UTC times carry the same quoted timing uncertainty. Ingress and egress use mid-transit ± 1.215 h and do not include duration uncertainty. Future UTC conversions assume the current leap-second table; regenerate after any leap-second announcement. The ephemeris assumes a constant period and does not include possible transit-timing variations.

Regenerate with NumPy and Astropy:

```sh
python code/ephemeris.py data/fit/chain.npz ephemeris.csv
```

## Reproduction and scope

The scripts are an archive of the analysis, not a self-contained pipeline. Absolute private home paths have been replaced by relative `work/planet-hunt/` paths. Adapt these to your own working directory; some scripts also expect intermediate search catalogues, stellar-prior files, original light curves and third-party tools that are not bundled. The included posterior is sufficient to regenerate the ephemeris independently.

The analysis uses NumPy, SciPy, Astropy, pandas, lightkurve, Transit Least Squares, wotan, batman, emcee, astroquery, Matplotlib, TRICERATOPS and [MOLUSC](https://github.com/woodml/MOLUSC). FITS light curves and TESScut cutouts are not included; the observations are public at [MAST](https://mast.stsci.edu/). No virtual environments or third-party source code are bundled.

## Credits

This work uses data collected by the TESS mission, funded by the NASA Science Mission Directorate, and obtained from MAST, including TESS-SPOC and QLP light curves and TESScut cutouts. It also uses data from the ESA Gaia mission, processed by the Gaia Data Processing and Analysis Consortium, and the TESS Input Catalog. Full references are in the research note.

## Licence and citation

Code in `code/` is under the [MIT licence](LICENSE). Text, figures and data are under [CC BY 4.0](LICENSE-CC-BY-4.0). Citation metadata is in [CITATION.cff](CITATION.cff). Version 1.0 records this candidate and its analysis; the repository and release dates are public timestamps, not a claim of peer review or validation.

## Contact

Advait Paliwal — advaitpaliwal@gmail.com
