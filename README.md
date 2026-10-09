# TIC 375939912.01: a 2.06 Earth-radius, 4.30-day transiting planet candidate

Discovered by Advait Paliwal on 2026-10-09 in public TESS full-frame-image photometry.

## Summary

TIC 375939912 (TYC 525-1773-1) is a K dwarf at 84.5 pc. TESS full-frame-image photometry from Sectors 55 and 82 shows a transit-like signal with period 4.3003975 d, depth 868 ppm and duration 2.43 h, implying a radius of 2.06 R⊕. The signal appears in raw and corrected flux, with consistent depths in TESS-SPOC and QLP. Difference images exclude all 23 Gaia neighbours within 60″ as its source. As of 2026-10-09, the star is not a TOI or CTOI.

With pixel exclusions and MOLUSC Gaia companion limits, TRICERATOPS gives a false-positive probability of **0.0075 ± 0.0019**, below the 0.015 statistical-validation threshold. No high-resolution imaging exists, and unresolved companions remain possible. **This is a strong planet candidate, not a validated or confirmed planet.**

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
| False-positive probability | 0.0075 ± 0.0019 with pixel exclusions and Gaia limits |

The symmetric period and epoch errors above use the larger side of the fitted interval. Full asymmetric intervals and convergence diagnostics are in `data/fit/result.json`. The false-positive probability is conditional on the adopted model and constraints; its quoted scatter is across ten runs. The false-positive runs used the search depth of 762.3 ppm, while the later transit fit gives 868 ppm; see `data/fpp/trc_summary.txt`.

## Contents

- [Research note](note/TIC375939912_note.md): methods, results, limitations and references.
- [Figure 1](note/figure1_TIC375939912.png): folded TESS-SPOC photometry, transit model and residuals.
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
