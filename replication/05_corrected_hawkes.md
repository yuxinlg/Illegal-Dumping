# Reporting-Corrected Cox-Hawkes Fits by Planning District

Documentation for `[05_corrected_hawkes_by_district.py](05_corrected_hawkes_by_district.py)`,
the supporting model module `[cox_hawkes_offset.py](cox_hawkes_offset.py)`, and the
street-light-only variant driver `[05b_corrected_hawkes_streetlight.py](05b_corrected_hawkes_streetlight.py)`.

Two runs, each in its own folder under `output/` (figures mirror the layout
under `output/figures/`):

| run | NT instrument | decomposition scope | Cox-Hawkes fits | outputs |
| --- | --- | --- | --- | --- |
| main (`05`) | 04 spec A: 5 curated categories | citywide (reuses 04's existing outputs) | all 18 districts | `output/corrected_hawkes/` |
| street-light (`05b`) | **`Street Light Outage` only** — Alley Light Outage, Traffic Signal Emergency, Stop Sign Repair, Dead Animal in Street removed | **single citywide fit** (no district layer; poles as NT exposure) — see below | **Central** (Center City) and **University Southwest** (University City) only | `output/corrected_hawkes_streetlight/` |

> **Note (2026-07-13).** The 05b run was refactored and is now documented in
> its own note, `05b_streetlight_citywide.md`. Its decomposition is no longer
> a per-district refit of 04: it is ONE citywide two-equation fit over all
> CBGs (single β, no `corrected_level`/ā anchor), with the street-pole
> inventory as the NT-equation exposure and PSIP robustness checks. For that
> run the two-layer offset of §2 collapses to just the within-CBG layer,
> `offset_i = β·ŵ_R,i`. §2's per-district derivation continues to describe the
> main (`05`) run only.

---

## 1. Context & intuition

The per-district Cox-Hawkes fits (`03a_analysis_by_division.py`) model the
intensity of *observed* 311 dumping reports, which conflates two things:

> **observed report intensity = (reporting propensity) × (true dumping intensity)**

$$\lambda_{obs}(s,t) = p(s)\cdot\lambda_{true}(s,t)$$

`04_reporting_decomposition.py` already *estimates* the reporting-propensity
side per CBG (fields $w_R$, loadings $\beta_d$) from non-trash 311 counts.
This script implements 04's §7 "downstream hook": it feeds that estimate into
the Cox-Hawkes model as a **known thinning term** — on the log scale, a fixed
additive **offset** in the background intensity — so that everything else the
model fits ($f_{xy}$, covariate weights, $a_0$) describes the
**reporting-corrected true-dumping process**.

**How this differs from the old** `reporting_rate` **covariate.** Before, the
reporting signal entered as one regressor among ~20, with a *fitted* weight
free to be shrunk, confounded, or split against correlated covariates — the
fitted surface still mixed the two components. Here the correction is
*structural*: the offset's coefficient is pinned at 1 (each unit of log
propensity removes exactly one unit of log observed intensity), and
`reporting_rate` is removed from the covariate list so the signal cannot
enter twice.

Three specifications are fit for every district:


| spec           | offset term                                                     | purpose                                                                                                                                                                                                                     |
| -------------- | --------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `baseline`     | none                                                            | comparison fit (same covariates, minus `reporting_rate`)                                                                                                                                                                    |
| `offset_fixed` | $1\cdot\text{offset}(s)$                                        | **main spec**: true thinning                                                                                                                                                                                                |
| `offset_free`  | $\gamma\cdot\text{offset}(s)$, $\gamma \sim \mathcal N(1, 0.5)$ | robustness: if the point process itself estimates $\gamma \approx 1$, the correction is validated by independent data (the *spatiotemporal pattern* of dumping reports, which 04 never saw — 04 only used CBG count totals) |


---



## 2. The offset: computed per district, anchored citywide



### Equation

$$\text{offset}*i = \log p_i = \beta_d\big(\hat w*{R,i} + \hat a_{NT,d} - \bar a\big)$$


| Parameter       | Meaning                                                                                                                                 | Source (04 spec A outputs)                                                |
| --------------- | --------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------- |
| $\hat w_{R,i}$  | CBG $i$'s log reporting-propensity deviation **from its district mean** (04 centres it at 0 within each district); posterior mean       | `reporting_decomp_cbg.geojson`, column `w_R_mean_a`                       |
| $\hat a_{NT,d}$ | district $d$'s log NT-report rate per unit exposure (its "propensity level"); posterior mean                                            | `reporting_decomp_district_summary.csv`, `a_NT_mean` (spec `a`)           |
| $\bar a$        | citywide log NT rate $\log(\sum_i NT_i / \sum_i E_i)$ — one fixed constant, the common anchor                                           | recomputed from `reporting_decomp_cbg.geojson` (`NT_curated`, `exposure`) |
| $\beta_d$       | district $d$'s shared-factor loading: how strongly general 311 propensity carries into *dumping* reports; posterior mean                | `reporting_decomp_district_summary.csv`, `beta_mean` (spec `a`)           |
| $\log p_i$      | log reporting-thinning: the part of CBG $i$'s log observed-report intensity attributable to reporting behavior, on a **citywide scale** | computed here; saved as `log_p_report`                                    |

*(Street-light run: same equation, but every ingredient comes from the
`spec = "sl"` refit in `output/corrected_hawkes_streetlight/` —
`w_R_mean_sl` from `reporting_decomp_cbg_sl.geojson`, `beta_mean` /
`a_NT_mean` from `reporting_decomp_district_summary_sl.csv`, and $\bar a$
from street-light counts only.)*




### Why two layers ("per district, but city scale")

04's decomposition was fit separately per district, absorbing each district's
average propensity into its intercept $a_{NT,d}$ — so $\hat w_{R,i}$ alone is
only meaningful *within* a district. A CBG's full log propensity on one
citywide ruler is the sum of:

```
(district layer)          (within-district layer)
â_NT,d − ā          +          ŵ_R,i
```

both scaled by that district's $\beta_d$. Skipping the $\bar a$ centring was
the known 04 gotcha: $\beta_d$ varies by district, so $\beta_d$ times the raw
*level* of $a_{NT,d}$ (≈ −4 to −8) would make the offset ranking mechanical
in $\beta$.

**Role of each layer in a per-district fit.** Within one district's fit the
district layer $\beta_d(\hat a_{NT,d}-\bar a)$ is the same constant for every
CBG — it cannot change the fitted spatial *shape*. The design intent was that
it would shift that district's $a_0$, making the 18 fitted $a_0$'s
reporting-corrected levels on a common scale. **The full run shows this does
not happen in practice** (see §8, finding 2): $a_0$'s tight $\mathcal N(0,
0.5)$ prior resists large constant shifts and the spatial/temporal GPs have
a free level, so the constant is spread across $a_0$, $f_{xy}$ and $f_t$
rather than landing in $a_0$. Cross-district *levels* should therefore still
be read from 04's `corrected_level_mean`; the Cox-Hawkes contribution is the
corrected within-district **shape** (driven by the CBG-varying layer
$\beta_d\hat w_{R,i}$), the corrected covariate weights and excitation
shares, and the $\gamma$ test.

**Worked example.** District A chatty ($\hat a_{NT,A}-\bar a = +0.7$,
$\beta_A = 1.2$); district B quiet ($-0.7$, $\beta_B = 0.6$). A CBG in A with
$\hat w_R = +0.3$ gets offset $1.2(0.3+0.7) = +1.20$: the model treats
$e^{1.2}\approx 3.3\times$ of its raw intensity as reporting inflation and
discounts it. A CBG in B with $\hat w_R = 0$ gets $0.6(0-0.7) = -0.42$: its
raw counts *understate* dumping and the background is credited upward.

**Data notes (current run).** $\bar a = -5.22$; all 1,338 CBGs matched a 04
estimate (0 filled). Offset range $[-7.85, +3.15]$, 98% within
$[-3.2, +2.2]$; the extreme negative tail is the known exposure-artifact
CBGs (airport / stadium / industrial special-use blocks in Lower South and
Lower Southwest: huge SafeGraph foot traffic, ~zero residents and reports).
They are kept as-is (not clipped) — they contain (almost) no events, so they
influence only the intensity integral, and 04 already flags their propensity
levels as shaky. $\hat w_R$ enters as a plug-in posterior mean; offset
uncertainty is not propagated (deliberate v1 simplification).

---



## 3. Model specification

Identical to `03_analysis.py`'s `Cox_Hawkes_Shared` setup (same priors,
`Temporal_Power_Law` trigger, windows 0.8 / 0.025, SVI with lr 0.01,
20k steps, years 2021–2025) except for the background intensity:

$$\log \mu(s,t) = a_0 + f_t(t) + f_a(t) + f_{xy}(s) + X(s)^\top w + \gamma\text{offset}(s)$$

with total intensity $\lambda(s,t) = \mu(s,t) + \sum_{t_j < t} \alpha
k_t(t-t_j)k_s(s-s_j)$ (the Hawkes self-excitation part, unchanged).


| Term                     | Meaning                                                                                                                                                                           |
| ------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| $a_0$                    | scalar background intercept (prior $\mathcal N(0, 0.5)$). Under `offset_fixed` it becomes the district's reporting-corrected log background level                                 |
| $f_t(t)$                 | long-run temporal GP (trend), VAE-approximated                                                                                                                                    |
| $f_a(t)$                 | seasonal (annual-cycle) GP, VAE-approximated                                                                                                                                      |
| $f_{xy}(s)$              | spatial GP on the 25×25 grid — the residual spatial surface; under `offset_fixed` it is a **corrected** dumping surface                                                           |
| $X(s)^\top w$            | CBG covariates (z-standardized) × fitted weights, prior on $w$ per bstpp default. `reporting_rate` is **excluded** in all three specs                                             |
| $\gamma\text{offset}(s)$ | the reporting-thinning term, **raw scale** (never standardized). $\gamma \equiv 1$ (`offset_fixed`) or $\gamma \sim \mathcal N(1, 0.5)$ (`offset_free`), sampled as `offset_coef` |
| $\alpha$                 | self-excitation reproduction rate (prior Beta(2, 4))                                                                                                                              |


**Interpreting $\gamma$.** $\gamma = 1$: the decomposition's thinning is
exactly consistent with the spatiotemporal point pattern. $\gamma < 1$: the
point process wants only part of the correction; $\gamma > 1$: more than all
of it. Because 04 was fit on count *totals* and the Cox-Hawkes sees the
*pattern*, $\gamma \approx 1$ would be a non-circular validation. Two
caveats when reading $\hat\gamma$: (i) it is identified only from the
CBG-varying layer (the district-constant layer trades off freely with
$a_0$/GP levels regardless of $\gamma$); (ii) it is estimated *conditional
on* ~19 covariates and a flexible spatial GP — any part of $\hat w_R$
correlated with those is already absorbed, which biases $\hat\gamma$ toward
0 relative to the structural thinning claim. $\hat\gamma$ is best read as
"how much *additional* within-district reporting signal the point pattern
demands," a lower bound on the correction, not a referendum on 04.

**v1 caveat — background-only thinning.** Observation thinning applied
rigorously would also affect the Hawkes term (unobserved events cannot
trigger; observed offspring are thinned too). Here the trigger is untouched,
so $\alpha$ partially absorbs reporting differences in offspring
observation; interpret changes in the excitation share qualitatively.
Trigger thinning is the natural v2.

---

## 4. How the two models are embedded for a single district (walkthrough)

Using **Central** (Center City) as the example — the same steps apply to any
district. The two models never share parameters or a likelihood; they are
chained through **one number per CBG**, computed from the first model and
supplied to the second as fixed data.

**Step 1 — count model (04's decomposition), citywide.**
NT and T reports are counted per CBG for the whole city (street-light run:
NT = `Street Light Outage` only). For *each* of the 18 districts, the joint
two-equation Poisson model is fit on that district's CBGs, giving:
$\hat w_{R,i}$ for every CBG $i$ in Central (its log propensity relative to
Central's average), Central's propensity level $\hat a_{NT,d}$ and loading
$\beta_d$, and the citywide anchor $\bar a$. All 18 districts are fit even
when only two will get a Cox-Hawkes fit, because $\bar a$ and the district
levels are citywide quantities.

**Step 2 — the bridge: one offset column.**
$\text{offset}_i = \beta_d(\hat w_{R,i} + \hat a_{NT,d} - \bar a)$ is
computed for every CBG in the city and merged into the covariate
GeoDataFrame (`cov_cbg` + column `log_p_report`). This column is the *only*
thing passed from model 1 to model 2 — a plug-in posterior-mean estimate of
$\log p(s)$, treated as known data downstream.

**Step 3 — point-process model (bstpp Cox-Hawkes), one district.**
Central's polygon becomes the study area; dumping *events* (not counts) inside
it, 2021–2025, are the data. When `Cox_Hawkes_Shared_Offset` is constructed
with the citywide covariate gdf, bstpp already builds two mappings:
`cov_ind` (each event → the CBG it falls in, by spatial join) and `int_df`
(each 25×25 background-grid cell → the CBG pieces it intersects, with areas).
The offset column is read **raw** (it bypasses the z-standardization applied
to covariates) in `cov_ind` order, and is added to the covariate term
$b_0 = X w$ so it enters the likelihood at both places $b_0$ does:

* every event in CBG $i$ contributes
  $\log\mu = a_0 + f_t + f_a + f_{xy} + X_i^\top w + \text{offset}_i$;
* the compensator (expected count) integrates
  $e^{\,b_0 + \text{offset} + f_{xy}}$ over the grid-cell × CBG
  intersections, weighted by their areas.

Nothing else in the Cox-Hawkes changes — same priors, trigger kernel, SVI.

**Step 4 — what the single-district fit does with it.**
Within Central, the CBG-varying part $\beta_d\hat w_{R,i}$ acts as a spatial
discount/credit: events in chatty CBGs are partly "explained" by the offset,
so $f_{xy}$ and the covariate weights no longer need to — the fitted surface
becomes the reporting-corrected one. The district-constant part
$\beta_d(\hat a_{NT,d}-\bar a)$ cannot affect the within-district shape and
gets absorbed into $a_0$/GP levels (§2). In the `offset_free` spec the same
column enters as $\gamma\cdot\text{offset}$ with $\gamma$ estimated — the
single-district test of how much of the decomposition's correction the point
pattern actually demands.

---



## 5. Code changes inventory (nothing existing was touched)


| File                                                                 | Status        | Contents                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| -------------------------------------------------------------------- | ------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `cox_hawkes_offset.py`                                               | **new**       | `spatiotemporal_hawkes_model_shared_offset` — copy of bstpp's model function (~/BSTPP @ 192a599) with the offset added at **both** required points: the event log-intensity (`b_0[cov_ind]`) and the spatial integral (`b_0[int_df.cov_ind]`), via `b_0 = b_0 + γ·offset` so both stay consistent. `Cox_Hawkes_Shared_Offset` — subclass taking `spatial_offset=` (column name; read **raw** from the stored covariate gdf in `cov_ind` order, bypassing z-standardization) and `offset_coef_prior=`; overrides `get_params()` so `offset_coef` is captured in `model.samples`. |
| `05_corrected_hawkes_by_district.py`                                 | **new**       | offset construction, per-district 3-spec batch, outputs & figures (this doc §6–7); clones 03's `setup_and_fit_model` as `setup_and_fit_model_offset`. Writes into per-run folders via `RUN_NAME`/`RUN_OUT`/`FIG_OUT` globals, and its `DECOMP_*` globals let a driver point it at a different decomposition run                                                                                                                                                                                                                                                                  |
| `05b_corrected_hawkes_streetlight.py`                                | **new**       | street-light-only driver: loads 04 as a module and overrides `NT_CURATED = ["Street Light Outage"]` / `NT_SPECS = {"sl": "NT_curated"}` **in memory** (04's functions read them at call time — no file edit), refits the decomposition citywide, saves it under `output/corrected_hawkes_streetlight/`, then loads 05 as a module, points its `DECOMP_*`/`RUN_*` globals at the refit and sets `UNIT_SUBSET = ["Central", "University Southwest"]`                                                                                                                             |
| `05_corrected_hawkes.md`                                             | **new**       | this file                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| `~/BSTPP/*`, `03_analysis.py`, `03a_analysis_by_division.py`, `04_*` | **unchanged** | 03/03a functions are reused by loading them as modules; the bstpp model function is copied, not edited (trade-off: the copy does not follow future ~/BSTPP updates — re-diff before reusing on a newer bstpp)                                                                                                                                                                                                                                                                                                                                                                   |


Why no bstpp edit was needed: `Point_Process_Model.__init__` stores the full
covariate GeoDataFrame as `self.spatial_cov` (only `cov_names` columns are
standardized into `args`), and `run_svi` fits whatever `self.model` holds,
which `Cox_Hawkes_Shared` assigns *after* `super().__init__()` — so a
subclass can swap in the offset-aware model function and extra args
post-construction.

---



## 6. Data outputs

Each run gets its own folder: `output/corrected_hawkes/` (main) and
`output/corrected_hawkes_streetlight/` (street-light variant). File names
below are relative to the run folder. The street-light folder additionally
contains the refit decomposition itself: `reporting_decomp_cbg_sl.geojson`
(same schema as 04's CBG file, spec suffix `_sl`) and
`reporting_decomp_district_summary_sl.csv` (one row per district,
`spec = "sl"`).

### `offset_cbg.geojson` — one row per CBG


| Column               | Definition                                                                                |
| -------------------- | ----------------------------------------------------------------------------------------- |
| `GEOID`, `dist_name` | 2023 CBG id; planning district (from 04's assignment)                                     |
| `log_p_report`       | $\beta_d(\hat w_{R,i} + \hat a_{NT,d} - \bar a)$ — the offset exactly as fed to the model |




### `fit_summary.csv` / `.geojson` — one row per district × spec


| Column                                                 | Definition                                                                                                                             |
| ------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------- |
| `dist_name`, `spec`, `n_events`                        | district; `baseline` / `offset_fixed` / `offset_free`; events in fit                                                                   |
| `prop_excitation_mean`                                 | posterior mean of $I_{excite}/I_{total}$ — share of intensity from self-excitation                                                     |
| `alpha_mean`, `alpha_p_pos`                            | Hawkes reproduction rate: posterior mean; $P(\alpha > 0)$                                                                              |
| `a_0_mean`, `a_0_lo`, `a_0_hi`                         | background intercept: posterior mean and 94% CI. Under `offset_fixed`: the **reporting-corrected, citywide-comparable** district level |
| `w_mean_{cov}`                                         | posterior-mean weight per covariate (standardized scale)                                                                               |
| `{trigger}_mean`                                       | posterior means of trigger-kernel parameters (picked up generically)                                                                   |
| `offset_coef_mean`, `offset_coef_lo`, `offset_coef_hi` | $\hat\gamma$ and 94% CI — `offset_free` rows only                                                                                      |
| `corrected_level_mean`                                 | 04's decomposition-side corrected level, merged in for the validation scatter                                                          |
| `status`, `error`, `runtime_sec`                       | fit bookkeeping (03a pattern)                                                                                                          |




### `cbg_fxy.geojson` — one row per CBG × spec


| Column      | Definition                                                                                                                                                 |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `f_xy_norm` | fitted spatial GP surface aggregated to CBGs, min-max normalized **within district** (03a's `extract_cbg_fxy`); `spec` column distinguishes the three fits |


---



## 7. Figures (`output/figures/{run}/`)

`ch_map_offset_cbg.png` **— the offset itself.**
$\text{color}_i = \log p_i$. Blue (negative): quiet CBGs whose raw counts
understate dumping — the model credits their background upward. Red:
propensity-inflated CBGs, discounted. This is the input, mapped for sanity.

`ch_map_f_xy_norm_baseline.png` **/** `ch_map_f_xy_norm_offset_fixed.png`**.**
$\text{color}*i = \tilde f*{xy,i} \in [0,1]$ within district. The before /
after spatial surfaces (compare within districts only, per 03a's caveat).

`ch_map_f_xy_norm_delta.png` **— where the correction moves the surface.**
$\text{color}*i = \tilde f^{corr}*{xy,i} - \tilde f^{base}_{xy,i} \in [-1,1]$.
Positive: CBG becomes a *stronger* within-district hotspot once reporting
propensity is stripped out (dumping was masked by quiet phones); negative:
raw hotspot partly explained away by chatty reporters. Note both surfaces
are separately min-max normalized within district, so this is a shift in
*relative* within-district position, not in absolute intensity.

`ch_map_a_0_offset_fixed.png` **— fitted intercepts (diagnostic only).**
$\text{color}*d = \hat a*{0,d}$ under `offset_fixed`. Designed as a
citywide-comparable corrected level; the run shows $a_0$ does **not**
capture the offset's constant layer (§8, finding 2), so read this only as a
fit diagnostic — cross-district levels come from 04's figure 4.

`ch_scatter_a0_vs_corrected_level.png` **— cross-method check (negative result).**
$(x_d, y_d) = (\text{04's correctedlevel}*d, \hat a*{0,d})$. Intended
validation scatter; observed $r \approx -0.05$ — flat. Diagnosis: the
district-constant part of the offset does not land in $a_0$
($\operatorname{corr}(\Delta a_0, -\text{district layer}) = -0.53$, wrong
sign): the tight $a_0$ prior and the GPs' free level spread the constant
across components. A level-preserving summary (e.g. mean fitted log
background over the district) would be needed to do this comparison
properly — noted as future work, not attempted in v1.

`ch_scatter_excitation_shift.png` **— does correction change contagion?**
$(x_d, y_d) = (\text{stat}^{base}_d, \text{stat}^{corr}_d)$ for
`alpha_mean` and `prop_excitation_mean`, dashed $y=x$. Points off the
diagonal show where apparent self-excitation was partly reporting structure
(spatially clustered propensity mimics triggering) — or vice versa.

`ch_dumbbell_cov_weights.png` **— covariate weights before/after.**
For each covariate, posterior-mean weight averaged over districts, baseline
vs `offset_fixed`. Weights that shrink under correction were partly proxying
reporting propensity; weights that grow were masked by it.

`ch_forest_offset_coef.png` **— headline validation.**
$\hat\gamma_d$ with 94% CI per district, dashed line at $\gamma = 1$.
Intervals covering 1 ⇒ the point pattern independently endorses 04's
thinning; systematic deviation flags districts where the correction over- or
under-shoots (read together with 04's own $\beta$ forest, figure 9, and the
two $\hat\gamma$ caveats in §3). *Current run:* $\hat\gamma$ ranges −0.30 to
0.69, no CI covers 1 — the pattern demands a **partial** correction, and its
geography mirrors 04's $\beta$ story: $\hat\gamma$ is smallest (even
negative) exactly in the heavy-dumping districts (Upper North −0.30, West
−0.19, Lower Southwest −0.11) and largest in the lower-dumping Northeast /
river districts (Lower Northeast 0.69, River Wards 0.64, West Park 0.57).

---



## 8. Verification & results (runs of 2026-07-09, years 2021–2025)

### Main run (spec-A offset, 18 districts)

**Checks (all passed).**

1. **Equivalence check** (Lower Far Northeast): fitting `offset_fixed` with
  an all-zero offset reproduced the baseline posteriors **exactly** (worst
   |diff| = 0.0e+00 across a_0, α, all covariate weights, trigger
   parameters, loglik) under bstpp's fixed SVI seed — both injection points
   (event term and spatial integral) are consistent.
2. **Offset hand-check**: `log_p_report` recomputed manually for 3 CBGs from
  raw 04 outputs — exact match. All 1,338 CBGs matched a 04 estimate
   (0 filled). $\bar a = -5.2202$.
3. **Smoke test** (one district, 3 specs): all fit; `offset_coef` present in
  samples only for `offset_free`.
4. **Full run**: 54/54 fits successful (18 districts × 3 specs), 37 min
  total on CPU. `git status` confirms `~/BSTPP`, `03_analysis.py`, `03a`
   untouched.

**Findings.**

1. **The point pattern endorses a partial correction, with 04's geography.**
  $\hat\gamma$: mean of district means 0.26, range −0.30 to 0.69; no 94% CI
   covers 1 (2 cover 0). The ordering reproduces 04's β story from
   independent information: lowest/negative γ in the heavy-dumping districts
   (Upper North, West, Lower Southwest — where visible dumping forces its
   way into 311 regardless of propensity), highest in lower-dumping
   districts (Lower Northeast 0.69, River Wards 0.64, West Park 0.57).
   Given the conditioning caveat (§3 — GP + covariates absorb shared signal,
   pushing γ toward 0), read this as: *at minimum* roughly half the 04
   correction survives in the low-dumping districts, and in heavy-dumping
   districts within-district raw counts are already nearly
   propensity-neutral. It does **not** support applying less of the
   correction for cross-district comparisons, which the γ test cannot see.
2. **a_0 is not a usable cross-district corrected level** (negative result):
  the offset's district-constant layer does not land in $a_0$
   ($\operatorname{corr}(\Delta a_0, -\text{layer}) = -0.53$; a_0 vs 04
   corrected_level r = −0.05). Cause: tight $a_0 \sim \mathcal N(0,0.5)$
   prior + free GP levels. Cross-district levels remain 04's job.
3. **Correction shifts contagion estimates down**: α mean −0.11, excitation
  share −0.06 (baseline → offset_fixed) — part of the apparent
   self-excitation was spatially clustered reporting propensity.
4. **Surface changes are concentrated where expected** (delta map): largest
  negative shifts in the special-use extreme-offset CBGs (Lower South /
   Lower Southwest); positive shifts scattered through quiet-reporting
   blocks in North/West Philadelphia.

### Street-light run (`05b`: NT = Street Light Outage only; Cox-Hawkes for Central + University Southwest)

**Decomposition refit** (citywide, 18/18 districts, ~76k street-light
reports): $\bar a = -5.69$; $\beta_d$ range [0.38, 1.11] — Central (1.02),
University Southwest (1.11) and Lower North (1.03) now carry the highest
loadings. Offset range $[-7.50, +2.70]$, all 1,338 CBGs matched.

**Instrument robustness.** The street-light-only propensity field is nearly
identical to the 5-category spec-A field:
$\operatorname{corr}(\hat w_R^{sl}, \hat w_R^{a}) = 0.956$ citywide
(Central 0.936, University Southwest 0.965) — dropping the four other
categories barely changes what the instrument measures, extending 04's
figure-8 robustness down to a single category.

**Cox-Hawkes fits** (6/6 successful; baselines reproduce the main run's
exactly — same seed, offset-independent):

| district | $\hat\gamma$ street-light [94% CI] | $\hat\gamma$ spec A [94% CI] | α (base → fixed) | excit. share |
| --- | --- | --- | --- | --- |
| Central | **0.63** [0.59, 0.66] | 0.45 [0.41, 0.48] | 0.20 → 0.20 | ≈ 0 (both) |
| University Southwest | **−0.03** [−0.08, 0.02] | 0.07 [0.03, 0.12] | 0.20 → 0.20 | ≈ 0 (both) |

Reading: in Central the point pattern endorses substantially *more* of the
street-light-based correction than of the spec-A one (0.63 vs 0.45) — the
purer instrument survives conditioning on the GP/covariates better. In
University Southwest $\hat\gamma \approx 0$ under both instruments: within
that district, whatever propensity variation $w_R$ captures is already fully
absorbed by the covariates and spatial GP, so the offset adds no
*within-district* information there (the cross-district thinning is a
separate question the γ test cannot see, per §3's caveats). Both districts
have essentially no self-excitation (α ≈ 0.2 with ≈ 0 excitation share) in
baseline and corrected fits alike, so the correction leaves the contagion
conclusion unchanged here.

