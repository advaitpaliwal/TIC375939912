Changes made in Revision 1 of the research note. The revised note text is with the journal and is not in this repository; the figure, new runs, scripts and source files it cites are in this folder. Paths in this file are the author's working paths.

# AAS82123 (TIC 375939912): changes from submission2/rnaas2.md to rnaas_rev1.md

Paths are relative to ~/Scratch/planet-hunt/. Word counts: abstract 151; body 997 (all sections except abstract and references, counting acknowledgments and the caption; 903 without them). Same convention as revision-82088 (990).

## Bar check (plain statement)

- At the fitted depth with NO neighbours removed: FPP = 0.0337 ± 0.0018 (10 runs, range 0.0308–0.0368), NFPP = 0.0281 ± 0.0008. Source: `revision-82123/fpp_fitdepth/summary.json`.
- **The mean FPP passes the ≤0.035 bar, but only just: 7 of 10 runs are ≤0.035 and 3 are above it (0.0351, 0.0354, 0.0368).** 0 of 10 runs are below 0.015, so there is no validation claim.
- NFPP = 0.028 is 3.8× the 0.0074 found at the 762 ppm search depth. It is above 10⁻³. It comes from two neighbours, TIC 375939920 (58.8%) and TIC 375939917 (41.2%) (`summary.json` key `nfpp_share_top`).
- FPP is now mostly NEB (0.0272) + NEBx2P (0.0009), then STP (0.0041). At 762 ppm it was NEB 0.0071, STP 0.0067, NEBx2P 0.0003 and smaller terms. In TRICERATOPS 1.1.0, FPP = 1 − TP − PTP − DTP, so PTP and DTP do not count toward it.

## The new FPP run (what was done)

- Script `revision-82123/fpp_fitdepth/run.py` is a copy of `prerev/375939912/fpp/run.py` (mode `full`, no neighbours removed, 221 stars). Unchanged: target pickle, Gaia background, MOLUSC file `deep/375939912/fpp/molusc_out/molusc_for_trc_N1e6.csv`, 10-min bins over |phase| < 0.25 d (72 bins), err = SD of out-of-transit bins, `filt='TESS'`, N = 10⁶, `exptime=30/1440`, `nsamples=20`, seeds 1000–1009, TRICERATOPS 1.1.0 (same venv).
- Three changes, all needed to run "at the fitted depth":
  1. `tdepth` = the fitted 867.986 ppm (`deep/375939912/fit/result.json` `depth_ppm`), not 762.32 ppm.
  2. Light curve = TESS-SPOC PDCSAP detrended as in `deep/375939912/fit/prep.py`: wotan biweight 0.75 d, transits masked over 2.5 × 1.965 h, trend interpolated, clipped 5σ up and 10σ down. Built by `fpp_fitdepth/build_lc.py` → `lc_fitdetrend.npz` (12,906 points). Inside the fit window it reproduces `fitdata.npz` to 6.6 × 10⁻⁵ (1,615 of 1,615 points). This replaces the unmasked lightkurve search detrend that caused the shallow depth.
  3. Folded on the fitted ephemeris (P = 4.3003975 d, T0 = BTJD 3542.81698) instead of the search ephemeris.
- Box depth of the new binned input: 856.4 ppm, consistent with the fitted 868. Bin error 212.9 ppm (was 202.5).
- Outputs: `fpp_fitdepth/out/res_fit_*.json`, `probs_fit_*.csv`, `run_fit_*.log`; summary `fpp_fitdepth/aggregate.py` → `summary.json`.

## Second number: rule-based exclusion (added after the fitted-depth run)

- Pixel rule (AAS82088R1): summed Δχ² ≥ 25 AND no sector preferring the neighbour by ≥ 25. Checked from `prerev/375939912/loc/dchi2_sum_persector.csv`. TIC → Gaia DR3 matched on G and separation (pickle `GAIAmag`/`sep` vs the loc table):
  - TIC 375939920 = Gaia 1737322487506264320 (G 13.40 vs 13.39; 43.8″ vs 44.0″): S55 +22.6, S82 +37.9, Σ 60.5. Rule met.
  - TIC 375939917 = Gaia 1737322418786784256 (G 15.81 vs 15.80; 35.2″ vs 35.1″): S55 +21.0, S82 +37.4, Σ 58.3. Rule met.
- Run: `revision-82123/fpp_rule/run.py` = `fpp_fitdepth/run.py` plus only the removal block copied from `prerev/375939912/fpp/run.py` (mode sel), dropping the two TICs in `fpp_rule/drop_sel.txt` (221 → 219 stars). Same `lc_fitdetrend.npz` (copied), fitted ephemeris, tdepth 867.986 ppm, MOLUSC file, N = 10⁶, exptime 30 min, nsamples 20, seeds 1000–1009, TRICERATOPS 1.1.0. Binned input identical (box 856.4 ppm, err 212.9 ppm, 72 bins). All 10 runs at once (`go.sh`).
- Result (`fpp_rule/summary.json`): FPP = 0.0058 ± 0.0017 (range 0.0033–0.0086), 10 of 10 runs < 0.015; NFPP = 7.4 × 10⁻⁶ ± 0.4 × 10⁻⁶. Scenarios: TP 0.841, PTP 0.147, DTP 0.0064, STP 0.0043, BEBx2P 0.0007, BEB 0.0006, NEB 7.4 × 10⁻⁶. Remaining NFPP shares: TIC 1935544897 64%, 375939914 31%, 1935544908 5%.
- Text (FPP section only; abstract unchanged, as in the AAS82162 rev1): "but I have not removed them" → per-sector Δχ² ≥ 21.0 stated, then "Following Giacalone et al. (2021)… removing only these two gives FPP = 0.006 ± 0.002 (range 0.003–0.009) and NFPP = 7 × 10⁻⁶. These exclusions rest on TESS pixels, not on ground-based photometry." Closing line → "Only with the pixel-based exclusions does the FPP fall below 0.015, so I make no validation claim (cf. Collier et al. 2026)." Removed for length: the earlier 762 ppm result (FPP 0.0148 ± 0.0018, NFPP 0.0074; row 14 above) and "as in my earlier runs". STP 0.004 holds in both runs (0.00415, 0.00426).
- docx rebuilt with `tools/mkdocx.py` (run via `uv run --with python-docx`); 1 image; no-removal FPP precedes the exclusion FPP.

## Changes, section by section

| # | Where | Was | Now | Source |
|---|---|---|---|---|
| 1 | Author block | the postal address | removed (name, "Independent researcher", email only) | — |
| 2 | Abstract, FPP | "Difference images exclude all 23 Gaia sources… With these exclusions… FPP = 0.0075 ± 0.0019… below the 0.015 validation threshold" | "favour the target over each of the 23… With Gaia companion limits and no neighbours removed… FPP = 0.034 ± 0.002 and NFPP = 0.028 at the fitted depth… above the 0.015 validation threshold" | `fpp_fitdepth/summary.json` |
| 3 | Abstract | "strong candidate, not a validated planet" | "awaits ground-based photometry and high-resolution imaging" (the word "validated" is gone) | rule |
| 4 | Abstract / Star | coverage stated only in Summary | "the only two sectors in which TESS has observed it"; "no further coverage is predicted through Sector 134" moved to Star | `deep/375939912/sectors/coverage_audit.json` |
| 5 | Search (Q4) | no search detrending stated | lightkurve Savitzky–Golay, 0.75 d window, no mask, 3σ upper clip, 30-min bins, TLS 1–12 d, SDE 24.5 | `hunt5b/run1.py` lines 35–38, 46; `hunt5b/res/375939912.json` SDE 24.4577 |
| 6 | Fit (Q4) | "after detrending with wotan" | biweight, 0.75 d window; ±2.46 h masked (2.5 × 1.965 h / 2) and interpolated; 5σ/10σ clip; 1,615 points within 5.9 h (3 × 1.965 h) | `deep/375939912/fit/prep.py` lines 11–17; `hunt5b/vet_375939912.json` dur_h 1.9652; `fitdata.npz` n = 1615 |
| 7 | Fit | a = 0.0468 ± 0.0019 AU | a = 0.0468 (+0.0018/−0.0019) AU | `fit/result.json` a_au [0.04682, 0.00180, 0.00194] |
| 8 | Fit | stellar-density comparison in Vetting | moved to the end of Fit (2.5 ± 0.6 vs 2.9 ± 0.8 g cm⁻³, unchanged) | `fit/result.json` rho_star_gcc, rho_prior |
| 9 | Vetting (Q2) | depths given without method | lightkurve 1 d flatten, transits masked; depth = median OOT − mean in-transit | `deep/375939912/sectors/results.json` key `method` |
| 10 | Vetting | "raw and corrected fluxes give consistent depths" | "agree to within 0.5σ in each sector and pipeline" (max 0.46σ) | `hunt5b/event_gates_375939912.json` gate_a |
| 11 | Vetting (Q2) | odd/even "0.48σ… 0.73σ… 2.3σ, with the sign reversing between sectors and pipelines" | σ defined as difference over quadrature-summed errors; the sign statement is replaced by "in QLP Sector 55 the odd transits are the deeper ones" (TESS-SPOC has even ≥ odd in both sectors, so the old wording was wrong) | `fit/oddeven.json` 0.478; `sectors/results_qlp_combined.json` 0.734; `sectors/results.json` comb_oe_sigma 2.3143 and per-sector odd/even |
| 12 | Vetting (Q5) | none | "the four in Figure 1 have S/N 2.4–6.7" | `prerev/375939912/fig/fig_stats.json` zoom snr 5.39, 2.44, 6.75, 5.21 |
| 13 | Vetting, localisation (Q1, Q2) | "fitted… with the TESS PRF… disfavoured by a summed Δχ² of at least 29.2, and none is preferred in either sector" | both sectors named (5 and 6 transits); 15 × 15 stamp; noise from 120 random-phase difference images; free amplitude + constant; Δχ² = χ²(neighbour) − χ²(target); Σ ≥ 29.2; each sector separately ≥ 6.5 | `hunt5b/localize/375939912/out/diff_s55.json`, `diff_s82.json` (n_events 5, 6); `hunt5b/localize/diff.py` lines 17–48; `prerev/375939912/loc/dchi2_sum.csv` (min sum 29.23, min per-sector 6.53, 23 rows, max sep 56.5″) |
| 14 | FPP (Q3) | "FPP = 0.049 ± 0.002 and NFPP = 0.0069 without constraints. Removing the 23 excluded neighbours… FPP = 0.0075 ± 0.0019 (range 0.0048–0.0109). All ten runs fall below the 0.015 threshold." | all 220 neighbours kept; FPP 0.034 ± 0.002 (0.031–0.037), NFPP 0.028; the two NFPP carriers named with their Σ Δχ² 60.5, 58.3 "but I have not removed them"; earlier 762 ppm no-removal result 0.0148 ± 0.0018 / 0.0074 also given | `fpp_fitdepth/summary.json`; `prerev/375939912/fpp/summary_full.json`; `prerev/375939912/fpp/perstar_nfpp_full.csv` |
| 15 | FPP | "Modelling the exposure time as 600 s or 200 s… gives 0.0068–0.0072" | removed (it was a with-removal result); the text now says exposure time 30 min, as before | `prerev/375939912/fpp/run.py` exptime |
| 16 | FPP (Q6) | "companion limits from Gaia astrometry and imaging via MOLUSC" | Gaia RUWE and contrast tests only; 200,000 generated, 52% kept; kept companions replace the default companion prior in bound-companion scenarios ("imaging" removed: it meant the Gaia contrast curve) | `prerev/375939912/methods/molusc_counts.json` (kept_frac 0.521615); PREREV.md §(f) |
| 17 | FPP | "Almost all of the remaining probability is a planet orbiting an unresolved bound companion. Gaia cannot exclude companions inside about 1 AU." | "Most of the rest of the FPP (0.004) is a planet transiting an unresolved bound companion" (NEB now dominates FPP; the 1 AU figure had no source file and was dropped) | `summary.json` scenarios STP 0.00415 |
| 18 | FPP | "strong candidate pending imaging" | "The FPP is above 0.015, so I make no validation claim (cf. Collier et al. 2026)" | |
| 19 | Summary | imaging and ground photometry, in that order | ground photometry resolving the two NFPP stars first, then imaging; "most useful next steps" removed | |
| 20 | Figure 1 | one panel (folded) | three panels: folded; unfolded Sector 82 with wotan trend; four Sector 82 transits | `prerev/375939912/fig/figure1_combined_TIC375939912.png` (copied unchanged), made by `prerev/tools/make_fig.py` |

Unchanged numbers checked against files: T 10.58 (ExoFOP 10.5847), Teff 4635 ± 121, R★ 0.708 ± 0.052 (srad 0.707854, e 0.0522), M★ 0.74 ± 0.09, 84.5 pc (84.456), RUWE 0.98 (`hunt5b/gaia_flags.csv` 0.976), 2,220 stars (`hunt5b/pool_counts.json`), every fit value (`deep/375939912/fit/result.json`), per-sector depths 956 ± 117, 769 ± 92, 919 ± 113, 903 ± 106 and injections 673–840 ppm, phase-0.5 −38 ± 77 ppm (`sectors/results.json`, comb_phase0.5 −37.7 ± 77.1).

## How each of the six AAS82088 questions is answered

1. Sectors in the difference analysis: both sectors that exist (S55, S82), each with its event count, and per-sector minimum Δχ² stated.
2. χ² definitions: localisation Δχ² defined in full; odd/even σ defined; per-sector depths stated in ppm ± error with method.
3. Neighbour removal: none. The lead number is the no-removal FPP at the fitted depth.
4. Detrending: the search (lightkurve SG) and the fit (wotan biweight) are described separately with windows and masks.
5. Unfolded light curve: Figure 1 middle panel, and single transits in the bottom panels.
6. MOLUSC: tests used, counts, and the scenarios it affects.

## Files

- `rnaas_rev1.md`, `Paliwal_TIC375939912_RNAAS_rev1.docx` (built with `tools/mkdocx.py`, the same builder as the 82088 R1; 1 image), `figure1_combined_TIC375939912.png`, `response.txt` (draft letter, for use if asked).
