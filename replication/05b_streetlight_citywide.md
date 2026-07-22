# Citywide Street-Light Decomposition with Pole-Inventory Exposure

Documentation for `05b_corrected_hawkes_streetlight.py`.
This note records what changed in this move, why, the full model, and the
validation evidence. It supersedes the street-light-variant description in
`05_corrected_hawkes.md` §1 (the "05b" column of its run table): the sl run is
no longer a per-district refit of 04 — it is a single citywide fit with a new
exposure for the NT equation.

---

## 1. Context

Observed 311 illegal-dumping counts conflate two things: how much dumping
truly happens, and how willing residents are to call 311. The decomposition
(04) separates them with a second report stream that shares the *willingness*
component but carries no dumping signal: street-light-outage reports are filed
by the same residents through the same channel, so cross-CBG variation in them
pins down the latent reporting propensity $w_R$; whatever variation is left in
dumping reports after crediting $w_R$ is the reporting-corrected dumping field
$w_T$.

Two weaknesses in the previous implementation motivated this move:

1. **Layered pooling.** The model was fit independently in each of the 18
   planning districts (own $\beta_d$, $a_{NT,d}$, $a_{T,d}$, $\sigma$'s), and
   CBG results were only made citywide-comparable through an anchored scalar
   $a_{T,d} - \beta_d(a_{NT,d} - \bar a)$ plus a two-layer offset — with a
   known gotcha (the $\bar a$ centring) and a documented negative result (the
   district-constant offset layer does not land in the Cox-Hawkes $a_0$).
2. **Wrong denominator for outages.** Both equations used SafeGraph foot
   traffic as exposure, so $w_R$ was identified from "outage reports per
   foot-traffic unit". But the population at risk of an outage is *lights that
   can break*, not visits. A CBG with few poles and lots of foot traffic would
   look like a low-propensity CBG mechanically.

**This move:** one citywide fit over all Philadelphia CBGs (no district layer
anywhere in the decomposition), with the pole inventory
(`data/Street_Poles.geojson`, 203,036 poles) as the NT-equation exposure.

## 2. Intuition

- With a single fit, $w_{R,i}$ and $w_{T,i}$ are deviations from the
  *citywide* level: every CBG is directly comparable to every other, with no
  within-district-only caveat, no anchor scalar, no centring gotcha. The
  cross-district level variation formerly absorbed by $a_{NT,d}$ now lives
  inside $w_R$ itself (so its range widens relative to the old fit — that is
  signal, not a bug).
- Putting poles in the NT denominator makes the instrument read as it should:
  "given how many lights this CBG has, how many outage reports does it
  generate?" Under the identifying assumption below, that ratio is pure
  reporting propensity.
- The Cox-Hawkes offset collapses from the two-layer
  $\beta_d(\hat w_{R,i} + \hat a_{NT,d} - \bar a)$ to simply
  $\beta\,\hat w_{R,i}$: with one citywide $\beta$ and one citywide intercept,
  there is no district-constant layer left to carry.

### 2.1 Mechanism: one equation, two unknowns

For each CBG we observe one number — dumping reports — but believe it is the
product of two things we cannot see separately:

$$\text{observed dumping reports} \;=\; (\text{true dumping}) \times (\text{willingness to report})$$

With only this equation, a CBG with few reports is ambiguous: is it clean, or
is it a place where nobody bothers to call 311? Two unknowns per CBG, one
measurement — unidentifiable, no matter the model.

The fix is a second measurement that shares exactly **one** of the unknowns.
Under identifying assumption (ii) — lights truly break at the same per-pole
rate everywhere —

$$\text{outage reports} \;=\; \text{poles} \times (\text{citywide breakage constant}) \times (\text{willingness to report})$$

contains **no CBG-specific "true rate" unknown**: breakage is a constant,
poles are data, and the only CBG-varying term left is willingness. This
equation is therefore solvable on its own — a CBG's outage reports per pole,
relative to the city's, *is* its reporting propensity $w_R$. Plugged back into
the dumping equation, only one unknown remains, and the leftover is $w_T$.
The whole design is: *use the stream with one unknown to solve for the shared
unknown, then subtract it out of the stream with two.*

In log form (what the code fits, §4), the anatomy is visible directly:

$$
\begin{aligned}
\log \mathbb{E}[NT_i] &= \log P_i + a_{NT} + w_{R,i}
  && \leftarrow \text{one unknown per CBG: solvable} \\
\log \mathbb{E}[T_i]  &= \log E_i + a_T + \beta\, w_{R,i} + w_{T,i}
  && \leftarrow w_R \text{ known from above; leftover} = w_T
\end{aligned}
$$

**Why $\beta$?** Willingness to report a broken streetlight and willingness to
report dumping are related but not identical behaviors; $\beta$ is the
exchange rate between them, fit from how strongly outage-per-pole rates and
dumping rates co-move across the 1,335 CBGs. This run found $\beta = 0.58$: a
CBG filing $2\times$ the average outage reports tends to file about
$2^{0.58} \approx 1.5\times$ the dumping reports *at equal true dumping*.

### 2.2 A worked example (this run's fitted numbers)

Citywide constants from the fit: outage rate
$e^{\hat a_{NT}} = e^{-1.08} \approx 0.34$ reports per pole over 2021–25;
dumping rate $e^{\hat a_T} \approx 0.0052$ per foot-traffic unit;
$\beta = 0.58$. Take two hypothetical CBGs that look **identical** in the
dumping data — same foot traffic (10,000 stops), same 52 dumping reports, both
exactly at the citywide average rate — differing only in outages:

| | poles | outage reports | outages/pole vs city | $w_R$ | $e^{\beta w_R}$ | predicted dumping | observed | $w_T$ |
|---|---|---|---|---|---|---|---|---|
| CBG A ("chatty") | 100 | 60 | $1.76\times$ | $+0.57$ | $1.39$ | $\approx 72$ | $52$ | $-0.33$ |
| CBG B ("quiet") | 100 | 17 | $0.50\times$ | $-0.69$ | $0.67$ | $\approx 35$ | $52$ | $+0.40$ |

The arithmetic, spelled out for CBG A: outages per pole $60/100 = 0.60$
against the citywide $0.34$ gives the ratio $1.76$, so
$w_R = \ln 1.76 = +0.57$; the propensity multiplier entering the dumping
equation is $e^{\beta w_R} = e^{0.58 \times 0.57} = 1.39$; propensity alone
therefore predicts $52 \times 1.39 \approx 72$ dumping reports, and the
shortfall gives $w_T = \ln(52/72) = -0.33$. (CBG B mirrors it:
$w_R = \ln 0.50 = -0.69$, prediction $\approx 35$,
$w_T = \ln(52/35) = +0.40$.)

Same raw dumping counts, opposite conclusions. A's residents call 311 about
everything (76% above average on lights); given that chattiness we would
expect ~72 dumping reports and see only 52 — A is **truly cleaner** than it
looks. B's residents rarely call (half the average); given that silence, 52
reports is a lot — B is **truly dirtier**. Raw counts rank A and B equal; the
corrected field $w_T$ ranks B well above A. This reversal is the entire
purpose of the model.

### 2.3 Why a joint Bayesian fit instead of that arithmetic

The table above is the deterministic skeleton, and for large-count CBGs the
model does essentially that. The Poisson + priors machinery earns its keep at
the edges:

- **Small counts.** A CBG with 20 poles and 2 outage reports has a wildly
  noisy ratio; the prior $w_R \sim \mathcal N(0, \sigma_R^2)$ shrinks it
  toward the citywide average instead of taking $2/20$ literally (the bend in
  `sl_ppc_fit.png`).
- **Simultaneity.** $\beta$, the citywide constants, and all ~2,670 latent
  values are estimated together, so uncertainty in one propagates into the
  others (per-CBG `w_R_sd`, `w_T_sd` in the output).
- **No division by zero.** Everything is counts and offsets; zero-report CBGs
  are handled naturally.

One caveat the trick cannot escape: it corrects *relative* geography, not the
absolute amount of unreported dumping — a citywide uniform under-reporting
factor is invisible to any within-city comparison.

## 3. Identifying assumptions

(i) **Reporting behavior is uniform within a CBG** — the same latent
    propensity drives both outage reports and dumping reports from a CBG's
    residents.

(ii) **The true physical rate of streetlight breakage per pole is uniform
    across the city.** Under (ii), variation in reported outages *per pole*
    reflects reporting propensity, not real breakage — which is what lets
    outage reports proxy the reporting rate we divide out of dumping reports.

**The PSIP threat to (ii).** The Philadelphia Street Lighting Improvement
Program replaced old bulbs with LEDs, which fail far less. The pole file shows
the conversion landed mostly **inside the study window's tail**: of ~117k
dated relights (`light_date`), only ~21k predate 2024 vs ~96k in 2024–2026,
and the per-CBG LED share spans roughly 45%→75% (p10→p90). So during 2024–25
true breakage fell *non-uniformly* across CBGs — exactly the kind of real-rate
variation assumption (ii) rules out. Two robustness levers address it (§6):
a pre-conversion time slice (NT counts from 2021–2023 only, when the stock was
still mostly old bulbs) and an LED-share sensitivity check (does $w_R$
correlate with conversion geography?).

Design decisions fixed in advance: the denominator is **all poles regardless
of owner** (residents report outages without knowing whether Streets or PECO
owns the pole, and we have no 311 routing data); `pole_date` is **not** a
temporal filter (nearly all poles predate the window); the denominator is pole
count, not luminaire count (`nlumin` is missing for all PECO poles and
contains zeros/junk).

## 4. Model

One fit over all CBGs $i$ with $E_i > 0$ and $P_i > 0$ (counts pooled over
2021–2025):

$$
\begin{aligned}
NT_i &\sim \mathrm{Poisson}\!\left( P_i \, e^{\,a_{NT} + w_{R,i}} \right)
  && \text{(outage equation)} \\
T_i  &\sim \mathrm{Poisson}\!\left( E_i \, e^{\,a_T + \beta w_{R,i} + w_{T,i}} \right)
  && \text{(dumping equation)}
\end{aligned}
$$

- $P_i$ = street poles in CBG $i$ (exposure, coefficient fixed at 1)
- $E_i$ = SafeGraph average monthly stop counts (`alloc_avg_s_cnt`), unchanged
  from 04 — the dumping equation still measures dumping per unit of activity
- $NT_i$ = Street Light Outage 311 reports (single category; the 5-category
  spec-A instrument lives in 04's main run)
- $T_i$ = Illegal Dumping 311 reports

Priors (identical to 04; non-centred parameterisation for NUTS):

| parameter | prior | role / interpretation |
|---|---|---|
| $a_{NT}$ | $\mathcal N(0, 5^2)$ | citywide log outage-report rate per pole |
| $a_T$ | $\mathcal N(0, 5^2)$ | citywide log dumping-report rate per unit foot traffic |
| $\beta$ | $\mathrm{LogNormal}(0, 0.5)$ | loading of reporting propensity into dumping reports (one citywide value; formerly 18 $\beta_d$) |
| $\sigma_R$, $\sigma_T$ | $\mathrm{HalfNormal}(1)$ | spread of the two latent fields |
| $w_{R,i} = \sigma_R z_{R,i}$, $\; z_{R,i} \sim \mathcal N(0,1)$ | — | CBG $i$'s log reporting-propensity deviation from the citywide level |
| $w_{T,i} = \sigma_T z_{T,i}$, $\; z_{T,i} \sim \mathcal N(0,1)$ | — | CBG $i$'s log reporting-corrected dumping deviation — the target quantity |

Retired relative to the per-district spec: $\beta_d$, $a_{NT,d}$, $a_{T,d}$,
$\sigma_{\cdot,d}$ (18 of each → 1), the anchor
$\bar a = \log(\sum NT / \sum E)$, and the deterministic
$\text{corrected\_level}_d = a_{T,d} - \beta_d(a_{NT,d} - \bar a)$. With one
fit, cross-CBG comparability is built in and the corrected citywide level is
just $a_T$.

**Cox-Hawkes offset.** The reporting-thinning offset entering the background
log-intensity of the per-district Cox-Hawkes fits (05's machinery, unchanged)
becomes

$$\mathrm{offset}_i = \hat\beta \, \hat w_{R,i} \qquad \text{(plug-in posterior means)}$$

CBGs without a decomposition estimate get 0 (the citywide mean). The
district-constant layer of the old offset is gone by construction; the
`corrected_level` column consumed by 05's $a_0$-comparison figure is retired
(returned as NaN, which makes that figure skip itself — it was a documented
negative result anyway).

**Fitting.** NUTS (numpyro), `target_accept = 0.9`, 1000 warmup + 1000
samples, 1 chain, seed 0 — 04's settings, one run over ~1,325 CBGs
(~2,650 latent $z$'s) instead of 18 small runs. Point-process side: bstpp
Cox-Hawkes per district for Central and University Southwest, three specs
(baseline / offset_fixed / offset_free), exactly as in 05.

## 5. The pole inventory (EDA facts)

`data/Street_Poles.geojson`, 203,036 poles, EPSG:4326 → 26918.

- **Coverage**: 202,386 poles fall inside a 2023 CBG (650 outside the layer).
  Poles per CBG: median 131 (p5 54, p95 312). Only 3 CBGs have zero poles
  (excluded from the fit); 9 have <10 (kept — the Normal prior shrinks them).
- **Owner**: Streets 169,971 (83.7%), PECO 23,463 (11.6%), blank 6,878, plus
  13 minor owners (PennDOT, Private, CCD, …). PECO poles carry no attribute
  data: 98.7% `tap_id = 0`, 99.9% `psip_status = "N/A"`, `bulb_type UNKNOWN`,
  `nlumin`/`lum_size`/`height` missing. Spatially PECO is a small minority
  nearly everywhere — median CBG share 9%, p90 22%, only 4 CBGs above 50% —
  so PECO heterogeneity cannot drive cross-CBG variation. All poles stay in
  the denominator per the design decision.
- **Attributes**: `type` is dominated by WP (103k) / SNP (24k) / AEL (18k);
  `nlumin`/`lum_size`/`height` are 33.8% missing; `light_date` 42.2% missing;
  `block`/`plate` 91% missing. `nlumin` has 6,913 zeros and junk values
  (3614, 320) — reinforcing the pole-count (not luminaire-count) denominator.
- **PSIP**: `psip_status` COMPLETED 123,513 (60.8%), N/A 76,991 (mostly PECO +
  minor owners), PENDING/OTHER ~1.3k each. `bulb_type` LED 123,726 (60.9%).
  Dated relights: 2023 ≈ 11k, 2024 ≈ 65k, 2025 ≈ 30k — see §3.

Outputs: `poles_cbg.geojson` (per-CBG counts and LED/PECO/pre-2024 shares),
`pole_eda_*.csv` tables, and figures `pole_map_per_cbg.png`,
`pole_map_led_share.png`, `pole_hist_per_cbg.png`, `pole_light_date_hist.png`.

## 6. Validation and robustness results

Run of 2026-07-13. All correlations are in `validation_summary.csv`; scatters
under `output/figures/corrected_hawkes_streetlight/val_*.png`.

**Model fit.** 1,335 / 1,338 CBGs (excluded: 3 zero-pole CBGs, all in the
9892xx special-use tract; none lacked exposure). T = 116,751 dumping reports,
NT = 75,956 outage reports. NUTS: ~10 s, **0 divergences**.
$\hat\beta = \mathbf{0.579}$ $[0.470,\, 0.686]$ (94% CI) — comfortably inside
the old per-district range and below 1, echoing 04's finding that propensity
carries into dumping reports with a loading well under proportional.
$\hat a_{NT} = -1.077$ ($\approx 0.34$ outage reports per pole over the
window), $\hat a_T = -5.264$, $\hat\sigma_R = 0.620$, $\hat\sigma_T = 1.488$.
PPC in `sl_ppc_fit.png`; citywide maps `sl_map_w_R.png`, `sl_map_w_T.png`.

| check | what it tells us | pearson | spearman | n |
|---|---|---|---|---|
| $w_R$ citywide vs archived district-wise fit | how much the refactor changed the field (old fit was within-district centred, foot-traffic exposure) | 0.68 | 0.63 | 1,335 |
| $w_R$ vs log outages-per-pole | mechanical sanity check — should be ≈ 1 | 0.99 | 1.00 | 1,324 |
| $w_R$ full window vs 2021–23 NT window | PSIP time-slice robustness — high corr means assumption (ii) survives the LED conversion | **0.95** | **0.95** | 1,335 |
| $w_R$ vs LED share | PSIP sensitivity | 0.28 | 0.22 | 1,335 |
| $w_R$ vs pre-2024-relit share | PSIP sensitivity | −0.15 | −0.19 | 1,335 |
| log poles vs log foot traffic | how different the two denominators are | 0.46 | 0.51 | 1,335 |

**Decision rule** (pre-registered in the script): if
$|\mathrm{spearman}(w_R, \text{LED share})| > 0.3$, the 2021–23 NT window
becomes the headline spec. Result: $0.22 \le 0.3$, and the pre-conversion time
slice reproduces the full-window field at $r = 0.95$ (its
$\hat\beta = 0.588 \approx 0.579$) — **the full-window NT stands as the
headline spec; assumption (ii) is robust to the PSIP conversion.** The
moderate 0.63 correlation with the old district-wise fit is expected, not
alarming: the citywide field now carries the cross-district level variation
the old fit absorbed into its 18 intercepts, and the denominator changed
(poles vs foot traffic — themselves only 0.51 correlated).

**Cox-Hawkes** (Central + University Southwest × baseline / offset_fixed /
offset_free; 6/6 fits successful; offset range fed to the model
$[-1.26, +1.12]$, vs $[-7.8, +3.1]$ under the old two-layer offset):

| district | spec | $a_0$ | $\alpha$ | excitation share | $\gamma$ (94% CI) |
|---|---|---|---|---|---|
| Central | baseline | −1.28 | 0.20 | ~0 | — |
| Central | offset_fixed | −1.32 | 0.20 | ~0 | — |
| Central | offset_free | −1.32 | 0.20 | ~0 | **1.22 [1.15, 1.29]** |
| University Southwest | baseline | −1.50 | 0.20 | ~0 | — |
| University Southwest | offset_fixed | −1.51 | 0.20 | ~0 | — |
| University Southwest | offset_free | −1.50 | 0.20 | ~0 | −0.11 [−0.27, 0.05] |

The headline is Central's $\gamma$: under the old district-layer,
foot-traffic-based offset it was 0.63 (already the highest after conditioning
on the GP and covariates, which absorb shared signal and push $\gamma$ down);
with the citywide, pole-denominated offset the point process now estimates
$\gamma$ **at or above 1** — the spatiotemporal pattern of dumping reports,
data 04 never saw, independently endorses the correction as full thinning.
University Southwest stays ≈ 0 as before (University City's dumping pattern
carries little within-district reporting signal; its offset field is
comparatively flat). Both districts remain essentially free of
self-excitation, consistent with the earlier runs. The
`ch_scatter_a0_vs_corrected_level` figure is retired with `corrected_level`
itself (it documented a negative result).

## 7. Figure guide

All under `output/figures/corrected_hawkes_streetlight/`, in pipeline order.

### Pole inventory EDA (Step 0)

**`pole_light_date_hist.png` — PSIP relighting dates vs the study window.**
Bar chart of `light_date` years, 2021–2025 shaded. The visual statement of the
assumption-(ii) threat: relights are negligible until 2023, then ~11k (2023),
~65k (2024), ~30k (2025) — the LED conversion happened *inside* the window's
tail, which is why the robustness checks in §6 exist.

**`pole_map_per_cbg.png` / `pole_hist_per_cbg.png` — the NT denominator.**
Poles per CBG as choropleth + histogram. Coverage is dense and even (median
131/CBG, no holes) — what makes poles a usable exposure. The 3 gray zero-pole
CBGs (airport/navy-yard tract) are the excluded ones.

**`pole_map_led_share.png` — the geography of the PSIP threat.** Share of each
CBG's poles with `bulb_type = LED`. If this map mirrored the $w_R$ map,
"propensity" would partly be "conversion timing". It is patchy — low-LED
pockets in parts of North Philly, Grays Ferry, near the airport — but does not
mirror $w_R$ (spearman 0.22 < 0.3 threshold).

### Citywide decomposition (Step 1)

**`sl_map_w_R.png` — reporting propensity, one citywide scale.** Log deviation
from the citywide outage-per-pole rate ($+1 \approx 2.7\times$ the citywide
reporting rate, $-1 \approx 0.37\times$). Unlike the old within-district maps,
colors are comparable anywhere on the map. High propensity through Center
City, the South Philly rowhouse core, scattered Northeast pockets;
systematically depressed ($-1.5$ to $-2$) across the Lower Southwest and the
industrial riverfronts.

**`sl_map_w_T.png` — the target: reporting-corrected dumping, citywide
comparable.** After stripping willingness-to-report, the elevated band runs
across North Philadelphia and Kensington, with pockets in West and Southwest
Philly; the far Northeast is uniformly low. Put this next to the raw
dumping-rate map to show what the correction changes.

**`sl_ppc_fit.png` — posterior-predictive check.** Observed vs fitted (log-log),
one panel per equation; both hug the 45° line over 3+ orders of magnitude. The
bend at the bottom-left (fitted > observed for CBGs with 0–3 reports) is
partial pooling shrinking tiny-count CBGs toward the citywide mean — intended,
not pathological (0 divergences).

### Proxy validation & PSIP robustness

**`val_scatter_wR_vs_outage_per_pole.png` — mechanical sanity check.**
$w_R$ vs raw log(outages/pole); spearman ≈ 1.0 by construction ($w_R$ *is* the
shrunk, intercept-removed version of this ratio). Departures from the line are
the low-count CBGs being shrunk.

**`val_scatter_wR_full_vs_pre24.png` — the key robustness figure.**
$w_R$ with NT from 2021–2023 only (pole stock still mostly old bulbs) vs $w_R$
from the full window: a tight 45° band, spearman = 0.945, no curvature or
fanning from LED-converted CBGs. Even though true breakage fell non-uniformly
in 2024–25, it barely moves the propensity field — assumption (ii) survives.

**`val_scatter_wR_vs_led_share.png` — no-refit sensitivity.** Diffuse cloud,
spearman 0.22. Had this been strong, "high propensity" would really have meant
"not yet converted". It isn't → full-window spec stands (§6 decision rule).

**`val_scatter_wR_city_vs_district.png` — what the refactor did.** Archived
district-wise $w_R$ (foot-traffic exposure, within-district centred; note its
extreme $-7..+5$ range from the old exposure artifacts) vs the new citywide
estimate (compressed to about $\pm 2$). $r = 0.63$: clearly related,
deliberately not identical — the citywide fit re-injects the cross-district
level differences the old fit hid in its 18 intercepts, and the denominator
changed (poles vs foot traffic, themselves only 0.51 correlated).

### Cox-Hawkes step (Central + University Southwest)

**`ch_map_offset_cbg.png` — the offset fed to the point process.** Red =
above-average propensity (intensity thinned *down*), blue = suppressed
reporting (scaled *up*); range $[-1.26, +1.12]$ vs the old offset's
$[-7.8, +3.1]$. **Known cosmetic caveat:** the title still prints the old
two-layer formula $\beta_d(\hat w_R + \hat a_{NT,d} - \bar a)$ — it is
hardcoded in 05 (untouched); the data plotted is the new $\beta\,\hat w_R$.

**`ch_map_f_xy_norm_baseline.png` / `ch_map_f_xy_norm_offset_fixed.png`.**
The fitted GP spatial surface $f_{xy}$ aggregated to CBGs, min-max normalized
within each district, before/after correction — where *within each district*
the background intensity concentrates (03a's caveat applies: compare within a
district only).

**`ch_map_f_xy_norm_delta.png` — corrected − baseline.** Red = a CBG becomes a
relatively stronger hotspot once reporting is corrected; blue = it was partly
a reporting artifact. The core of Center City turns blue (high-propensity
residents inflated its apparent intensity) while Central's northern and
southern fringes turn red; University Southwest shifts more mildly the same
way. Magnitudes are small ($\pm 0.1$ on a 0–1 scale): the correction reshapes
hotspots at the margin, it does not redraw them.

**`ch_map_a_0_offset_fixed.png`.** Background intercept $a_0$ under the main
spec (Central −1.32, USW −1.51). Record-keeping with two districts; per the §8
negative result in `05_corrected_hawkes.md`, do not read $a_0$ as a
cross-district corrected level.

**`ch_scatter_excitation_shift.png`.** $\alpha$ and excitation share, baseline
vs corrected, one point per district. Both districts sit on the diagonal at
excitation ≈ 0: no meaningful self-excitation, unchanged by the correction
(unlike the 18-district spec-A run, where correction lowered excitation).

**`ch_dumbbell_cov_weights.png`.** Posterior-mean covariate weights, baseline
vs offset_fixed, averaged over the two districts. Covariates whose dots move
were partly proxying for reporting propensity in the baseline.

**`ch_forest_offset_coef.png` — the validation headliner.** offset_free lets
the data scale the offset by $\gamma$; $\gamma = 1$ (dashed line) = the
correction is exactly right. Central: $\gamma = 1.22$ $[1.15, 1.29]$ — at or
above full thinning, up from 0.63 under the old offset: the spatiotemporal
dumping pattern (data the decomposition never saw) independently endorses the
pole-based correction. University Southwest: $\gamma = -0.11$ $[-0.27, 0.05]$
≈ 0 — its point pattern carries little within-district reporting signal (its
offset field is comparatively flat), same as before.

## 8. Files

| file | status | role |
|---|---|---|
| `05b_corrected_hawkes_streetlight.py` | rewritten | pole EDA + join, citywide two-exposure fit, PSIP robustness, offset override, Cox-Hawkes driver |
| `04_reporting_decomposition.py` | untouched | count/exposure builders reused as modules; its own per-district spec-A run is unaffected |
| `05_corrected_hawkes_by_district.py` | untouched | Cox-Hawkes batch reused; its `build_offset` is replaced in memory by the citywide version |
| `cox_hawkes_offset.py`, `~/BSTPP`, 03/03a | untouched | hard constraint |
| `output/corrected_hawkes_streetlight/reporting_decomp_cbg_sl.geojson` | overwritten | citywide per-CBG `w_R`/`w_T` (+ `w_R_mean_pre24`, pole columns) |
| `…/reporting_decomp_city_summary_sl.csv` | new | two rows: full-window and pre-24 citywide scalars |
| `…/reporting_decomp_*_districtfit.*` | archived | the pre-refactor district-wise outputs, kept for the comparison in §6 |
