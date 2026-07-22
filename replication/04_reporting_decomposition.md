# Separating 311 Reporting Propensity from the True Illegal-Dumping Rate

Documentation for [`04_reporting_decomposition.py`](04_reporting_decomposition.py).

---

## 1. Context

The main analysis (`01`–`03`) models 311 "Illegal Dumping" reports with a
Cox-Hawkes spatiotemporal point process. The fundamental measurement problem
is that

> **observed dumping reports = (true dumping) × (willingness to report 311)**.

A neighborhood can show many reports because it has a lot of dumping *or*
because its residents readily call 311 — the two are observationally
equivalent from trash reports alone. Before this script, the only correction
in the pipeline was a `reporting_rate` covariate (non-dumping 311 ÷ SafeGraph
stops, built in `02_covariates.py` §7b) entering the Cox-Hawkes background
intensity as one regressor among many; the fitted intensity therefore still
conflates the two components.

This script implements an explicit **decomposition**: it estimates a latent
reporting-propensity surface $W_R$ and a latent true-dumping surface $W_T$
separately, per planning district (the same analysis geography as
`03a_analysis_by_division.py`), on CBG-level counts.

---

## 2. Intuition

**The trick: a second measurement that shares one unknown but not the other.**

Non-trash (NT) 311 reports — street light outages, broken traffic signals,
dead animals — are filed by the same residents through the same channel, so
they reflect the same *willingness to report* ($W_R$), but they contain **no
information about dumping**. Diagrammatically (the whiteboard sketch):

```
        W_T  (true dumping)
          │  arrows only to T
          ▼
   ┌──── T ────┐          in every spatial unit:
   │           │          T  = trash (Illegal Dumping) reports
   └─── NT ────┘          NT = non-trash reports
          ▲
          │  arrows to BOTH NT and T
        W_R  (reporting propensity)
```

- NT counts act as a **pure thermometer for reporting propensity**: a CBG
  that files many street-light complaints per unit of foot traffic simply
  has engaged 311 users.
- T counts are driven by **both** thermometers: propensity ($W_R$) *and*
  actual trash on the ground ($W_T$).
- Whatever variation in T is **left over** after crediting $W_R$ is
  attributed to $W_T$ — the reporting-corrected dumping rate.

**Worked example.** Two CBGs with identical foot traffic. CBG 1 files 100
street-light reports and 50 dumping reports; CBG 2 files 20 and 40. Raw
counts say CBG 1 is the hotspot (50 > 40). But CBG 1's NT rate shows its
residents report ~5× more readily; per unit of reporting propensity CBG 2
actually has ~4× the dumping (40/20 vs 50/100). The model formalizes this
ratio logic with proper count noise, partial pooling, and posterior
uncertainty.

**Key identifying assumption (state it in the paper).** Cross-CBG variation
in the NT categories reflects *reporting behavior*, not variation in the
true rate of those problems. This is why the main specification uses
infrastructure-failure categories — street lights burn out at roughly the
same rate everywhere — rather than potholes or abandoned autos, whose true
incidence tracks neighborhood conditions. The all-non-dumping specification
is the robustness check on exactly this assumption (Section 5, figure 8).

---

## 3. Technical specification

### Model

Fit **independently for each of the 18 planning districts** (mirroring 03a's
planning-district mode). Within district $d$, for each CBG $i$ whose
representative point falls in $d$, with counts pooled over 2021–2025:

$$NT_i \sim \text{Poisson}\big(E_i \cdot e^{\,a_{NT} + w_{R,i}}\big) \qquad \text{(eq. 1: NT pins down } W_R)$$

$$T_i \sim \text{Poisson}\big(E_i \cdot e^{\,a_T + \beta\, w_{R,i} + w_{T,i}}\big) \qquad \text{(eq. 2: leftover variation} = W_T)$$

$$w_{R,i} \sim \mathcal N(0, \sigma_R^2), \quad w_{T,i} \sim \mathcal N(0, \sigma_T^2), \quad \beta \sim \text{LogNormal}(0, 0.5), \quad \sigma_R, \sigma_T \sim \text{HalfNormal}(1)$$

Piece by piece:

| Term | Role |
|---|---|
| $E_i$ | **Exposure offset**, coefficient fixed at 1 (enters as $\log E_i$ inside the exponent). Converts counts to rates: 10× the foot traffic mechanically yields ~10× the reports at identical propensity/dumping. $E_i$ = SafeGraph avg monthly stop counts (`alloc_avg_s_cnt` in `cov_cbg.geojson`); population fallback if the column is absent. |
| $a_{NT}, a_T$ | **District intercepts.** Absorb the district-average rate of each report type, so $w_{R,i}, w_{T,i}$ are relative (log) deviations from the district mean — which is why the within-district 0–1 normalized mosaic is the natural display. |
| $\beta$ | **Shared-factor loading**: how strongly reporting propensity carries from NT into T. $\beta = 1$: a CBG filing 2× the street-light reports also files 2× the dumping reports at equal true dumping. Estimated, not pinned, so carry-over can be partial ($\beta<1$) or amplified ($\beta>1$). |
| $w_{R,i}, w_{T,i}$ | **Latent CBG fields** (log scale): $e^{w_{T,i}}$ = CBG $i$'s true-dumping rate relative to its district mean, holding propensity fixed. Normal priors give partial pooling: tiny-count CBGs shrink toward the district mean instead of producing wild ratios. |

**Identification.** $w_R$ is determined by eq. 1 alone (no other latent term
appears there); the two equations are fit *jointly*, so uncertainty in $w_R$
propagates honestly into $w_T$. District-average propensity is absorbed into
the intercepts — acceptable because per-district CBG fields are only
interpreted within-district (03a's own caveat for its $f_{xy}$ mosaic); the
cross-district comparison is done by a dedicated scalar, below.

**Cross-district scalar.** Let
$\bar a = \log\!\big(\sum_{\text{all }i} NT_i / \sum_{\text{all }i} E_i\big)$
be the citywide NT rate (a fixed constant). Under the identifying assumption,
a district's NT-rate deviation $(a_{NT,d} - \bar a)$ is pure propensity, so
the log reporting-thinning is $\log p_d = \beta_d (a_{NT,d} - \bar a)$ and

$$\text{corrected\_level}_d \;=\; a_{T,d} \;-\; \beta_d\,\big(a_{NT,d} - \bar a\big)$$

computed **per posterior sample** is the log reporting-corrected dumping rate
per unit exposure, comparable across districts.
*The $\bar a$ centring is essential:* $\beta_d$ varies by district, so without
it the scalar would be $\beta_d$ times the *level* of $a_{NT,d}$ (≈ −4 to −8)
and the ranking would be driven mechanically by $\beta$, not by dumping.

**Estimation.** numpyro NUTS (non-centred parameterisation,
`target_accept=0.9`), 1,000 warmup + 1,000 samples, one chain per
district × spec; ~3 s per fit on CPU. 18 districts × 2 NT specs = 36 fits.

**Why CBG-level data inside district-level fits?** 18 district counts alone
give two equations ≈ 18 observations — almost nothing to separate two latent
fields plus $\beta$. The within-district contrast between high-NT and low-NT
CBGs is what identifies $\beta$. Same philosophy as 03a: district-level
analysis consuming CBG-resolution inputs.

**Deliberate v1 simplifications.** Poisson likelihood (upgrade to negative
binomial if the posterior-predictive check shows overdispersion — it
currently does not); independent $w$'s rather than an ICAR spatial prior
(the district layer already provides coarse spatial structure).

**Synthetic-data check.** On simulated data with known truth
($\beta = 0.8$, $n = 70$ CBGs), the model recovers $\hat\beta = 0.81$,
$\operatorname{corr}(\hat w_R, w_R^{\text{true}}) = 0.96$,
$\operatorname{corr}(\hat w_T, w_T^{\text{true}}) = 0.87$.

### Data

| Ingredient | Definition | Source |
|---|---|---|
| **T** | "Illegal Dumping" 311 reports per 2023 CBG, 2021–2025 pooled (spatial join of report points to TIGER CBGs) | `data/311/{year}/public_cases_fc.shp` |
| **NT spec A** (main, curated) | `Street Light Outage`, `Alley Light Outage`, `Traffic Signal Emergency`, `Stop Sign Repair`, `Dead Animal in Street` — infrastructure failures with plausibly uniform true incidence (~120k reports) | same |
| **NT spec B** (robustness) | All non-dumping 311 **excluding** `Information Request` (phone-based, not a geolocated incident) and trash-adjacent categories overlapping the outcome (`Rubbish/Recyclable Material Collection`, `Sanitation Violation`, `Dumpster Violation`, `Sanitation / Dumpster Violation`) (~880k reports) | same |
| **Exposure $E_i$** | `alloc_avg_s_cnt` — SafeGraph avg monthly stop counts, area-weighted to 2023 CBGs | `output/cov_cbg.geojson` (from `02_covariates.py`) |
| **District assignment** | CBG representative point within district polygon (`dist_name`, 18 districts) | `data/Planning_Districts.geojson` |
| Validation (optional) | litter survey; PPR cleanup operations — checks skip gracefully when files are absent | `data/CSVFiles/litter_surveys_raw_20250730.csv`; `data/PPR_BaseFiles_IllegalDumpingModel.gdb` |

**Exposure caveat.** Foot traffic and *resident population* diverge in a few
districts (Lower South = stadium/airport: huge traffic, few residents), which
makes their propensity level $a_{NT}$ — and hence their corrected level —
shakier. A population-exposure re-run is a cheap third robustness check.

### Notation used below

$S = 1000$ posterior draws; a hat denotes the posterior mean,
$\hat x = \tfrac1S \sum_s x^{(s)}$. "Spec A/B" superscripts denote which NT
definition was used in the fit.

---

## 4. Data outputs (`replication/output/`)

### `reporting_decomp_cbg.geojson` — one row per CBG

| Column(s) | Equation / definition |
|---|---|
| `T_count`, `NT_curated`, `NT_all` | Raw pooled 2021–2025 report counts per CBG |
| `exposure` | $E_i$ |
| `dist_name` | Planning district of the CBG's representative point |
| `w_R_mean_{a,b}`, `w_R_sd_{a,b}` | $\hat w_{R,i}$ and posterior sd, per spec |
| `w_T_mean_{a,b}`, `w_T_sd_{a,b}` | $\hat w_{T,i}$ and posterior sd, per spec |
| `w_R_norm_{a,b}`, `w_T_norm_{a,b}` | Within-district min-max normalisation, e.g. $\tilde w_{T,i} = \dfrac{\hat w_{T,i} - \min_{j \in d}\hat w_{T,j}}{\max_{j\in d}\hat w_{T,j} - \min_{j\in d}\hat w_{T,j}}$ |
| `NT_fitted_{a,b}`, `T_fitted_{a,b}` | Posterior-mean fitted counts (see figure 10 equations) |
| `reporting_rate` | Existing covariate from `02_covariates.py` §7b, kept for validation |

### `reporting_decomp_district_summary.csv` / `.geojson` — one row per district × spec

| Column | Equation / definition |
|---|---|
| `n_cbg`, `T_total`, `NT_total` | Fit inputs |
| `raw_level` | $\log\big(\sum_{i\in d} T_i / \sum_{i\in d} E_i\big)$ — uncorrected benchmark |
| `a_NT_mean`, `a_T_mean` | $\hat a_{NT,d}$, $\hat a_{T,d}$ |
| `beta_mean`, `beta_lo`, `beta_hi` | $\hat\beta_d$ and 94% CI $[Q_{0.03}, Q_{0.97}]$ |
| `sigma_R_mean`, `sigma_T_mean` | Posterior means of the field scales |
| `corrected_level_mean`, `corrected_level_sd` | Posterior mean/sd of $a_{T,d} - \beta_d (a_{NT,d} - \bar a)$ — **the headline cross-district quantity** |
| `status`, `error`, `runtime_sec` | Fit bookkeeping (03a pattern) |

---

## 5. Figures (`output/figures/reporting_decomposition/`): equation + interpretation

### Headline CBG mosaics

**1. `rd_map_w_T_norm_cbg.png` — corrected dumping $W_T$, normalized within district.**

$$\text{color}_i \;=\; \tilde w_{T,i} \;=\; \frac{\hat w_{T,i} - \min_{j \in d}\hat w_{T,j}}{\max_{j \in d}\hat w_{T,j} - \min_{j \in d}\hat w_{T,j}} \in [0,1]$$

*Interpretation.* Within its own district, how much true dumping does each
CBG have after stripping out residents' eagerness to call 311? (0 = district
minimum, 1 = district maximum — do **not** compare colors across a district
boundary; figures 3–4 do that job.) Known dumping corridors (North
Philly/Kensington, industrial Southwest edges) survive the correction —
evidence they are not artifacts of chatty neighborhoods.

**2. `rd_map_w_R_norm_cbg.png` — reporting propensity $W_R$, same normalization.**

$$\text{color}_i \;=\; \tilde w_{R,i} \quad \text{(formula as above with } R \text{ for } T)$$

*Interpretation.* Which CBGs punch above their district in *willingness to
report anything*. Compare with figure 1: the surfaces differ visibly, direct
visual evidence the model isn't splitting one signal in half. Bright here +
dark in fig 1 ⇒ raw counts **overstate** dumping; dark here + bright in
fig 1 ⇒ raw counts **understate** it (blocks a raw hotspot analysis misses).

### District choropleths (main spec A)

**3. `rd_map_raw_level.png` — raw dumping level.**

$$\text{color}_d \;=\; \text{raw}_d \;=\; \log\!\left(\frac{\sum_{i\in d} T_i}{\sum_{i\in d} E_i}\right)$$

*Interpretation.* The "before" picture: reports per unit of foot traffic, no
correction. Bright band across North / Lower North / River Wards / West /
Upper North; far Northeast and Lower South darkest.

**4. `rd_map_corrected_level_mean.png` — reporting-corrected dumping level.**

$$\text{color}_d \;=\; \frac1S \sum_{s=1}^{S}\Big[a_{T,d}^{(s)} - \beta_d^{(s)}\big(a_{NT,d}^{(s)} - \bar a\big)\Big]$$

*Interpretation.* The "after" picture and the cross-district headline. Broad
structure survives (North / Lower North still peak; far Northeast still the
floor) but the middle of the ranking reshuffles — see figure 7. This is the
only district map designed to be comparable across districts under the
model's assumptions.

**5. `rd_map_a_NT_mean.png` — district propensity level.**

$$\text{color}_d \;=\; \hat a_{NT,d}, \qquad e^{\hat a_{NT,d}} \approx \text{NT reports per stop, district-average CBG}$$

*Interpretation.* District-scale reporting culture. Highest: Lower Southwest.
Lowest by far: Lower South — partly an exposure artifact (stadium/airport
foot traffic with few residents), which is why its corrected rank (16th) is
more believable than its suspiciously clean raw rank (18th).

**6. `rd_map_beta_mean.png` — shared-factor loading.**

$$\text{color}_d \;=\; \hat\beta_d \qquad \big(\text{a CBG with } e^{w_R} = 2 \text{ files } 2^{\beta} \times \text{the dumping reports at equal true dumping}\big)$$

*Interpretation.* Spatially structured, not noise. $\beta \approx 1.1$–$1.25$
in mostly better-off districts (Upper Northwest, Central, University
Southwest, Lower Southwest): dumping reports scale one-for-one with general
311 engagement, so raw counts there are heavily propensity-contaminated.
$\beta \approx 0.5$–$0.7$ in the heavy-dumping districts (Upper North, River
Wards, West, North): where dumping is severe and visible, even
low-propensity blocks report it — the signal partially forces its way
through, so raw counts are less biased there.

### Diagnostics, robustness, validation

**7. `rd_rank_raw_vs_corrected.png` — does the correction matter?**

$$(x_d, y_d) \;=\; \Big(\text{rank}_{\downarrow}(\text{raw}_d),\; \text{rank}_{\downarrow}(\text{corrected}_d)\Big), \qquad \text{dashed: } y = x$$

*Interpretation.* Spearman ≈ 0.90: the correction refines rather than
overturns (a total reshuffle would be suspicious). Movers are the policy
story — current run: Upper Northwest **falls 6** places (raw 7th →
corrected 13th; engaged reporters exaggerate its problem) while Lower
Southwest **rises 5** (13th → 8th) and Central rises 3: districts whose
dumping was masked by quiet phones. A cleanup allocation based on raw 311
would over-serve the former and under-serve the latter.

**8. `rd_spec_robustness.png` — is $W_T$ sensitive to the NT definition?**

$$(x_i, y_i) \;=\; \big(\hat w^{A}_{T,i},\; \hat w^{B}_{T,i}\big), \qquad r = \operatorname{corr}\big(\hat w^{A}_{T}, \hat w^{B}_{T}\big), \qquad \text{dashed: } y = x$$

*Interpretation.* Two very different NT definitions (5 curated categories vs
~30 categories, ~120k vs ~880k reports) yield nearly the same corrected
surface (current run: $r = 0.842$, $n = 1{,}338$). This answers the referee
question "aren't your results driven by which categories you called
non-trash?"

**9. `rd_beta_forest.png` — $\beta$ posteriors, both specs, all districts.**

$$\text{point} = \hat\beta_d, \qquad \text{bar} = \big[Q_{0.03}(\beta_d),\; Q_{0.97}(\beta_d)\big] \;(94\%\ \text{CI}), \qquad \text{dashed: } \beta = 1$$

*Interpretation.* (i) The two specs agree within uncertainty almost
everywhere (echoes fig 8). (ii) Districts split into a clearly-below-1 group
and an at-or-above-1 group — the geography of fig 6 is statistically real,
not posterior-mean noise. (iii) Spec B intervals are tighter simply because
it uses ~7× more NT reports.

**10. `rd_ppc_fit.png` — posterior-predictive fit (spec A).**

$$\hat T_i = \frac1S\sum_s E_i\, e^{\,a_T^{(s)} + \beta^{(s)} w_{R,i}^{(s)} + w_{T,i}^{(s)}}, \qquad \hat{NT}_i = \frac1S\sum_s E_i\, e^{\,a_{NT}^{(s)} + w_{R,i}^{(s)}}$$

plotted as $(T_i + 0.5,\ \hat T_i + 0.5)$ and $(NT_i + 0.5,\ \hat{NT}_i + 0.5)$
on log-log axes; dashed $y = x$.

*Interpretation.* Points hug the diagonal across three orders of magnitude —
no systematic Poisson misfit. The visible bottom-left deviation (0–2
observed reports fitted above the line) is partial pooling doing its job: a
CBG with one report over five years shouldn't be declared a pristine
outlier. Severe shrinkage at *high* counts would argue for the
negative-binomial upgrade; none is present.

**11. `rd_validation_reporting_rate.png` — sanity check against the existing covariate.**

$$(x_i, y_i) \;=\; \big(\text{reporting\_rate}_i,\; \hat w_{R,i}\big), \qquad r_s = \operatorname{corr}\big(\text{rank}(x), \text{rank}(y)\big)$$

where `reporting_rate` is the 02_covariates §7b measure (monthly non-dumping
311 ÷ allocated stops, averaged over months), log x-axis.

*Interpretation.* Reassuringly positive (current run: Spearman 0.64,
$n = 1{,}338$) — the latent $w_R$ recovers the same signal as the hand-built
ratio — but deliberately not 1.0: the covariate includes incidence-driven
categories that spec A excludes; $w_R$ is shrunk/pooled while the raw ratio
is noisy for small CBGs; and $w_R$ is measured relative to each district's
mean, which mechanically lowers a citywide correlation.

*(Two further ground-truth checks — Spearman correlations of $w_T$ vs
litter-survey scores and vs PPR cleanup counts, compared against the raw T
rate — run automatically when `data/CSVFiles/litter_surveys_raw_20250730.csv`
and `data/PPR_BaseFiles_IllegalDumpingModel.gdb` are present; they are
skipped gracefully otherwise.)*

---

## 6. What the figures collectively say (current run, 2021–2025 data)

The decomposition works and is stable: the corrected dumping surface keeps
the known hotspots (figs 1, 3, 4), is insensitive to the NT definition
(figs 8, 9: $r = 0.84$), fits the data (fig 10), and its propensity component
agrees with the existing reporting-rate measure (fig 11: $r_s = 0.64$). The
substantively new findings are the $\beta$ geography (figs 6, 9: raw counts
are most propensity-contaminated in well-off districts, least in
heavy-dumping ones) and the ranking corrections (fig 7: Upper Northwest
overstated by raw counts; Lower Southwest and Central understated).

## 7. Downstream hook (not implemented)

$\beta \cdot w_{R,i}$ is the estimated log reporting-thinning $\log p(s)$; it
can enter the per-district Cox-Hawkes fits (`03a_analysis_by_division.py`)
as a **fixed offset** in the background intensity, making each district's
$f_s$ a reporting-corrected dumping surface.
