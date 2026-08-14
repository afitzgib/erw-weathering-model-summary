NOTE ON FIGURES: I've added some figures to the repo. Note that there is still a bug visible at
50 degrees of latitute (data coarseness). Overall these figures are meant to catch these bugs
and show the differences between different model runs. In real life, they will be clipped to existing 
cropland (real deployment possibility), but are not clipped now to show the extent of spatial 
variation in weathering potential.

# Gross CDR comparison set — what changed, and what it did

Four maps, each isolating exactly one additional parameter/assumption on top
of the last, all sharing the same color scale (log, 0.05–12 tCO₂/ha/yr) so
color means the same absolute value in every map — the pattern shifting
between maps is the actual signal, not a rescaled axis. All built from one
representative basalt composition (the real median Ca/Mg across 96,202
kriged source samples, not a textbook assumption), fed through the
composition-weighted kinetics blend and real per-pixel stoichiometric
ceiling built earlier this session.

| File | What's added vs. the previous map | Global mean (tCO₂/ha/yr) | Global median |
|---|---|---|---|
| `gross_cdr_v1_year1.png` | Baseline: dissolution fraction × real ceiling, Year-1 single application. No lifecycle, CCE, leaching efficiency, or system loss. | 1.52 | 0.81 |
| `gross_cdr_v2_10yr_cohort.png` | + 10-year cohort-averaged reapplication (matches the live routing pipeline's convention instead of a Year-1 snapshot) | 4.83 | 3.40 |
| `gross_cdr_v3a_10yr_cohort_leffAnnual.png` | + leaching efficiency, live (stale) annual formula | 1.45 | 0.65 |
| `gross_cdr_v3b_10yr_cohort_leffMonthly.png` | + leaching efficiency, corrected monthly-aggregated formula | 1.74 | 0.89 |

**Reading the table:** cohort-averaging alone more than triples the global
mean (real accumulation effect, not a bug — see the Cascade doc §2.2b/§7 on
why this is the correct convention). Leaching efficiency then cuts that back
down substantially — even the *monthly* (more defensible) version still
lands below the pre-leaching-efficiency number, confirming leaching
efficiency is a real, large, currently-inactive penalty in the live model
(see the Cascade doc §3). The monthly version reads noticeably higher than
the stale annual one (1.74 vs. 1.45 mean) — consistent with the annual
formula under-crediting seasonal climates, exactly as documented.

**Two data bugs were found and fixed while building this set, and both are
baked into every map here** (previously they were baked into *every* map
this project has ever produced from this pipeline):

1. **CHIRPS only covers 50°S–50°N.** The moisture-modifier inputs (λ, α)
   mixed CHIRPS precipitation (latitude-limited) with CRU TS wet-days
   (global), which broke exactly at that boundary — visible as a hard
   latitudinal line in earlier maps. Fixed by falling back to CRU TS's own
   precipitation beyond CHIRPS' real range. **This is not cosmetic:**
   cropland above 50°N (Canada, Russia, Ukraine, Northern Europe) previously
   had dissolution crushed to near-zero (regional mean dissolution fraction
   0.0009) because feeding `alpha=0` into the moisture-modifier formula
   forces the modifier itself to exactly zero — worse than "no penalty," it
   zeroed out dissolution entirely for a huge share of global grain-belt
   cropland. With the fix, that same region's mean rises to 0.046 — roughly
   50×. The cutoff line is visibly gone from all four maps above.
2. **CHIRPS' files use −9999 as a missing-data sentinel but never declare it
   in their own metadata**, so averaging 30 years of files without excluding
   it corrupted ~4.3% of pixels globally (scattered, not boundary-specific —
   confirmed pixel-identical against the original live raster, so this
   predates and is independent of bug #1). Masked out post-hoc rather than
   re-derived from source, since the affected pixels are physically
   identifiable (precipitation-derived quantities can't be negative).

Full technical detail on both bugs and the fix: `scripts/run_fixed_chirps_pipeline.py`.

**A third, smaller, still-open issue — disclosed here rather than fixed —
concerns spatial resolution, not values.** CHIRPS' native resolution is
0.05°; CRU TS (used above 50°N/50°S per the fix above) is only 0.5°, a 10×
coarser grid. That means moisture-modifier inputs above 50°N/S are real,
climatologically valid, but spatially blockier than the rest of the map —
visible as a texture change, not a value discontinuity, right at the 50°
line. Quantifying both pieces of this directly rather than leaving it
qualitative:

| | ≤50° (CHIRPS-native, 0.05°) | >50° (CRU-TS-derived, 0.5°) |
|---|---|---|
| Cropland cells | 992,734 | 574,575 (37% of global cropland) |
| Gross CDR mean, v1 (Year-1) | 1.958 | 0.752 |
| Gross CDR mean, v2 (10-yr cohort) | 6.003 | 2.803 |
| Gross CDR mean, v3a (+leff annual) | 1.868 | 0.715 |
| Gross CDR mean, v3b (+leff monthly) | 2.205 | 0.922 |
| Mean storm depth (alpha), local sample | 3.96 mm | 3.92 mm |
| **Local std-dev of storm depth** (same sample, texture only) | **0.256 mm** | **0.072 mm** |

The lower CDR mean above 50° is real and expected — cooler high-latitude
climates genuinely dissolve rock slower (Arrhenius temperature scaling) —
**not** an artifact of the resolution mismatch. The resolution issue shows up
cleanly in the last row instead: mean storm depth is essentially identical
on both sides of the line (3.96 vs. 3.92 mm), but the *local* variability
(texture) is ~3.5× lower above 50°, because a single 0.5° CRU-TS pixel gets
repeated across many finer output pixels. This is a genuine data-availability
limitation (CHIRPS simply doesn't cover >50°N/S), not a bug — and it's
actively being closed out by pulling ERA5-Land (~0.1° native, genuinely
global) as a replacement for CRU-TS above 50°; that work is in progress.

**What this set deliberately does not yet include**, per the Cascade doc's
§7 open questions — worth raising alongside these figures, not silently
fixing first: CCE, system loss, and lifecycle emissions are all excluded by
design (this is "gross," not "net," CDR); and the moisture modifier (λ×α)
still has its own uncorrected, opposite-direction seasonal aggregation bias
relative to the now-fixed leaching efficiency (see Cascade doc §3/§7) — so
even the v3b (monthly leaching efficiency) map should be read as "one side
of a two-sided fix," not the final word on moisture-driven seasonality.
