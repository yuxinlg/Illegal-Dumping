"""
cox_hawkes_offset.py
====================
Reporting-corrected Cox-Hawkes model: `Cox_Hawkes_Shared` plus a fixed (or
free-coefficient) per-CBG offset in the background log-intensity.

    log mu(s,t) = a_0 + f_t + f_a + f_xy + X.w + gamma * offset(s)

where offset(s) is the log reporting-thinning log p(s) estimated by
05a_streetlight_decomposition.py (formerly 04's per-district decomposition,
now archived), and gamma is either pinned at 1 (a true
thinning/filter, the main spec) or given a prior (robustness spec: gamma ~ 1
validates the correction from the point-process side).

This module deliberately does NOT modify the bstpp package (editable install
at ~/BSTPP). It works because of two properties of bstpp's wiring:

  * `Point_Process_Model.__init__` stores the full covariate GeoDataFrame as
    `self.spatial_cov`; only the `cov_names` columns are z-standardized into
    `args['spatial_cov']`. A raw offset column therefore stays accessible
    untransformed, in `cov_ind` row order.
  * `run_svi` fits whatever `self.model` holds, and `Cox_Hawkes_Shared`
    assigns its model function after `super().__init__()`, so a subclass can
    swap in its own model function and extra args post-construction.

`spatiotemporal_hawkes_model_shared_offset` below is a copy of
`bstpp.cox_hawkes_shared.spatiotemporal_hawkes_model_shared`
(~/BSTPP @ commit 192a599, 2026-02-03) with the additions marked `# OFFSET`.
If ~/BSTPP is updated, this copy does not follow automatically — re-diff
against the source before reusing it on a newer bstpp.

Scope caveat (v1): the thinning is applied to the *background* intensity
only. The Hawkes trigger term is untouched, so alpha partially absorbs
reporting differences in offspring observation; interpret the excitation
share accordingly.
"""

import numpy as np
import jax
import jax.numpy as jnp
import numpyro
import numpyro.distributions as dist

from bstpp.cox_hawkes_shared import Cox_Hawkes_Shared, hawkes_intensity_sum
from bstpp.vae_functions import (
    vae_decoder_temporal,
    vae_decoder_seasonal,
    vae_decoder_spatial,
)


def spatiotemporal_hawkes_model_shared_offset(args):
    """
    Copy of bstpp.cox_hawkes_shared.spatiotemporal_hawkes_model_shared with a
    per-CBG background offset (coefficient fixed at 1 unless
    args['offset_coef_prior'] is a numpyro distribution). Changes are marked
    with `# OFFSET`.
    """
    hawkes    = args["hawkes_design"]
    t_events  = args["t_events"]
    xy_events = args["xy_events"]
    N         = t_events.shape[0]

    # intercept for background rate log(mu(t,s))
    a_0 = numpyro.sample("a_0", args['priors']['a_0'])

    ## ---- Temporal Gaussian Process ---------------------------------------------------

    # gaussian vector to feed into VAE
    z_temporal = numpyro.sample("z_temporal", dist.Normal(
        jnp.zeros(args["z_dim_temporal"]), jnp.ones(args["z_dim_temporal"])
    ))
    decoder_nn_temporal = vae_decoder_temporal(
        args["hidden_dim_temporal"],
        args["n_t"]
    )
    decoder_params = args["decoder_params_temporal"]

    # approximate Gaussian Process with VAE
    v_t = numpyro.deterministic("v_t", decoder_nn_temporal[1](decoder_params, z_temporal))
    f_t = numpyro.deterministic("f_t", v_t[0:args["n_t"]])
    rate_t = numpyro.deterministic("rate_t",jnp.exp(f_t + a_0))

    # calculate temporal integral over LGCP
    Itot_t = numpyro.deterministic("Itot_t", jnp.sum(rate_t)/args["n_t"]*args["T"])

    # temporal part of log(mu(t,s))
    f_t_events = f_t[args["indices_t"]]

    ## ---- Seasonal Gaussian Process ---------------------------------------------------

    # gaussian vector to feed into VAE
    z_seasonal = numpyro.sample("z_seasonal", dist.Normal(
        jnp.zeros(args["z_dim_seasonal"]), jnp.ones(args["z_dim_seasonal"])
    ))
    decoder_nn_seasonal = vae_decoder_seasonal(
        args["hidden_dim1_seasonal"],
        args["hidden_dim2_seasonal"],
        args["n_s"]
    )
    decoder_params = args["decoder_params_seasonal"]

    # approximate Gaussian Process with VAE
    v_a = numpyro.deterministic("v_a", decoder_nn_seasonal[1](decoder_params, z_seasonal))
    f_a = numpyro.deterministic("f_a", v_a[0:args["n_s"]])
    rate_a = numpyro.deterministic("rate_a",jnp.exp(f_a))

    # calculate integral over LGCP
    Itot_a = numpyro.deterministic("Itot_a", jnp.sum(rate_a)/args["n_s"]*args["S"])

    # seasonal part of log(mu(t,s))
    f_a_events = f_a[args["indices_a"]]

    ## ---- Spatial Gaussian Process ----------------------------------------------------

    # Generate gaussian vector to feed into VAE
    z_spatial = numpyro.sample("z_spatial", dist.Normal(
        jnp.zeros(args["z_dim_spatial"]), jnp.ones(args["z_dim_spatial"])
    ))
    decoder_nn = vae_decoder_spatial(
        args["hidden_dim1_spatial"],
        args["hidden_dim2_spatial"],
        args["n_xy"]
    )
    decoder_params = args["decoder_params_spatial"]

    # Generate Gaussian Process from VAE
    f_xy = numpyro.deterministic("f_xy",
        jnp.exp(args['sp_var_mu']) * decoder_nn[1](decoder_params, z_spatial)
    )
    f_xy_events = f_xy[args["indices_xy"]]

    ## ---- Spatial Covariates ----------------------------------------------------------

    # calculate spatial intensity
    if 'spatial_cov' in args:
        spatial_cov = args['spatial_cov']
        if spatial_cov.ndim == 1:
            spatial_cov = spatial_cov[:, None]
        args['spatial_cov'] = spatial_cov

        # weights for linear combination
        w = numpyro.sample("w", args['priors']['w'])
        b_0 = numpyro.deterministic("b_0", args['spatial_cov'] @ w)

        # OFFSET: reporting-thinning term log p(s), per CBG, raw scale
        # (bypasses the z-standardization applied to cov_names columns).
        # gamma pinned at 1 -> true thinning; prior given -> free coefficient.
        if 'spatial_offset' in args:                                     # OFFSET
            if args.get('offset_coef_prior') is not None:                # OFFSET
                gamma_off = numpyro.sample("offset_coef",                # OFFSET
                                           args['offset_coef_prior'])    # OFFSET
            else:                                                        # OFFSET
                gamma_off = 1.0                                          # OFFSET
            b_0 = b_0 + gamma_off * args['spatial_offset']               # OFFSET
        # b_0 now feeds BOTH the event log-intensity and the spatial
        # integral below, so the likelihood stays internally consistent.

        f_xy_events = f_xy_events + b_0[args['cov_ind']]
        spatial_integral = jnp.exp(
            b_0[args['int_df']['cov_ind'].values] + f_xy[args['int_df']['comp_grid_id'].values]
        ) @ args['int_df']['area'].values
    else:
      # OFFSET: the offset is defined at covariate (CBG) level; without
      # spatial_cov there is no cov_ind mapping to place it.
      if 'spatial_offset' in args:                                       # OFFSET
          raise ValueError("spatial_offset requires spatial_cov/cov_names.")  # OFFSET
      rate_xy = numpyro.deterministic("rate_xy",jnp.exp(f_xy))
      spatial_integral = jnp.sum(rate_xy[args['spatial_grid_cells']])/args['n_xy']**2

    # spatial part of log(mu(t,s)) (includes GP and spatial covariates)
    Itot_xy = numpyro.deterministic("Itot_xy", spatial_integral)

    # calculate total background integral
    Itot_txy_back = numpyro.deterministic("Itot_txy_back", Itot_t * Itot_a * Itot_xy)

    ## ---- Hawkes Process --------------------------------------------------------------

    # reproduction (self-excitation) rate
    alpha = numpyro.sample("alpha", args['priors']['alpha'])

    # spatial gaussian kernel parameters
    t_pars = args['t_trig'].sample_parameters()
    sp_pars = args['sp_trig'].sample_parameters()

    T,x_min,x_max,y_min,y_max = args['T'],args['x_min'],args['x_max'],args['y_min'],args['y_max']

    l_hawkes_sum = hawkes_intensity_sum(
        alpha, t_pars, sp_pars,
        hawkes["hawkes_i"],
        hawkes["hawkes_dt"],
        hawkes["hawkes_dx"],
        hawkes["hawkes_dy"],
        args["t_trig"],
        args["sp_trig"],
        hawkes["n_events"],
    )
    l_hawkes = numpyro.deterministic('l_hawkes', alpha * l_hawkes_sum)

    # background log intensity
    log_back = a_0 + f_t_events + f_a_events + f_xy_events

    # Hawkes log intensity — ensure non-negativity, clip tiny values
    log_hawkes = jnp.log(jnp.clip(l_hawkes, a_min=1e-12))

    # Stable mixture term
    log_intensity = jax.nn.logsumexp(jnp.stack([log_hawkes, log_back], axis=0), axis=0)

    ell_1 = numpyro.deterministic("ell_1", jnp.sum(log_intensity))

    #### hawkes integral
    temp_part = alpha * args['t_trig'].compute_integral(t_pars, T - t_events)

    sp_limits = jnp.stack((x_max-xy_events[0],xy_events[0]-x_min,
                           y_max-xy_events[1],xy_events[1]-y_min)
                         ).reshape(2,2,-1)

    sp_part = args['sp_trig'].compute_integral(sp_pars,sp_limits)

    Itot_excite = numpyro.deterministic("Itot_excite",jnp.sum(temp_part*sp_part))

    ## ---- Total for Cox-Hawkes --------------------------------------------------------

    # total integral
    Itot_txy = numpyro.deterministic("Itot_txy", Itot_excite + Itot_txy_back)
    loglik = numpyro.deterministic('loglik', ell_1 - Itot_txy)
    numpyro.factor("events", loglik)

    ## ---- Litter Survey Ordinal Logistic Regression -----------------------------------

    if "litter_hawkes_design" in args:
        litter = args["litter_hawkes_design"]
        log_bg_litter = (
            a_0 +
            f_t[args['litter_indices_t']] +
            f_a[args.get('litter_indices_a')] +
            f_xy[args['litter_indices_xy']]
        )
        l_hawkes_litter = hawkes_intensity_sum(
            alpha, t_pars, sp_pars,
            litter["hawkes_i"],
            litter["hawkes_dt"],
            litter["hawkes_dx"],
            litter["hawkes_dy"],
            args["t_trig"],
            args["sp_trig"],
            litter["n_events"],
        )
        l_hawkes_litter = jnp.clip(l_hawkes_litter, 1e-12)

        # combined background and hawkes intensity
        log_total = jax.nn.logsumexp(
            jnp.stack([log_bg_litter, jnp.log(l_hawkes_litter)]),
            axis=0,
        )

        # intercept, scale factor, and cutpoints
        gamma_l0 = numpyro.sample('gamma_l0', dist.Normal(0., 0.5))
        gamma_l1 = numpyro.sample('gamma_l1', dist.Normal(1., 0.5))
        cutpoints = jnp.cumsum(
            jnp.sort(
                numpyro.sample(
                    'cutpoints_unordered',
                    dist.Normal(jnp.array([-1., 0., 1.]), 0.5)
                )
            )
        )

        eta = jnp.clip(gamma_l0 + gamma_l1 * log_total, -20, 20)
        numpyro.sample(
            'litter_obs',
            dist.OrderedLogistic(eta, cutpoints),
            obs = args['litter_score'],
        )

    ## ---- Dumping Report Linear Regression --------------------------------------------

    if "dumping_report_hawkes_design" in args:
        dumping = args["dumping_report_hawkes_design"]
        log_bg_dumping = (
            a_0 +
            f_t[args['dumping_report_indices_t']] +
            f_a[args.get('dumping_report_indices_a')] +
            f_xy[args['dumping_report_indices_xy']]
        )
        l_hawkes_dumping = hawkes_intensity_sum(
            alpha, t_pars, sp_pars,
            dumping["hawkes_i"],
            dumping["hawkes_dt"],
            dumping["hawkes_dx"],
            dumping["hawkes_dy"],
            args["t_trig"],
            args["sp_trig"],
            dumping["n_events"],
        )
        l_hawkes_dumping = jnp.clip(l_hawkes_dumping, 1e-12)

        # combined background and hawkes intensity
        log_total = jax.nn.logsumexp(
            jnp.stack([log_bg_dumping, jnp.log(l_hawkes_dumping)]),
            axis=0,
        )

        # intercept, scale, and variance
        gamma_d0 = numpyro.sample('gamma_d0', dist.Normal(0., 0.5))
        gamma_d1 = numpyro.sample('gamma_d1', dist.Normal(1., 0.5))
        sigma_dumping = numpyro.sample("sigma_dumping", dist.HalfNormal(1.0))

        nu = numpyro.sample("nu_dumping", dist.Exponential(1.0))

        eta = jnp.clip(gamma_d0 + gamma_d1 * log_total, -20, 20)

        numpyro.sample(
            'dumping_report_obs',
            dist.StudentT(nu, eta, sigma_dumping),
            obs=args['dumping_report_log_cleanup_cost'],
        )


class Cox_Hawkes_Shared_Offset(Cox_Hawkes_Shared):
    """
    Cox_Hawkes_Shared with a per-CBG background offset.

    Extra keyword-only parameters (everything else passes through unchanged):

    spatial_offset : str | None
        Column name in the covariate GeoDataFrame holding the log
        reporting-thinning log p(s) per CBG. Read raw — never standardized.
        None reproduces Cox_Hawkes_Shared exactly (the offset-aware model
        function is still installed, but takes the identical code path).
    offset_coef_prior : numpyro distribution | None
        None (default) pins the offset coefficient at 1 (true thinning).
        A distribution (e.g. dist.Normal(1.0, 0.5)) makes it a sampled
        parameter, exposed in `model.samples["offset_coef"]`.
    """

    def __init__(self, *args, spatial_offset=None, offset_coef_prior=None, **kwargs):
        # captured here so they never reach bstpp's __init__, which would
        # otherwise swallow unknown kwargs into the priors dict
        super().__init__(*args, **kwargs)
        if spatial_offset is not None:
            if spatial_offset not in self.spatial_cov.columns:
                raise ValueError(
                    f"offset column {spatial_offset!r} not in spatial_cov"
                )
            # cov_ind row order — b_0 is indexed by cov_ind, so the offset
            # vector must follow the same ordering
            off_vals = (
                self.spatial_cov.sort_values("cov_ind")[spatial_offset]
                .to_numpy(dtype=float)
            )
            if np.isnan(off_vals).any():
                raise ValueError(
                    f"offset column {spatial_offset!r} contains NaN — fill "
                    "missing CBGs before fitting"
                )
            self.args["spatial_offset"] = jnp.asarray(off_vals)
            self.args["offset_coef_prior"] = offset_coef_prior
        self.model = spatiotemporal_hawkes_model_shared_offset

    def get_params(self):
        # run_svi builds its sample-site list from get_params(), so the free
        # coefficient must be registered here to appear in model.samples
        pars = super().get_params()
        if self.args.get("offset_coef_prior") is not None:
            pars["offset_coef"] = 1
        # bstpp's Hawkes_Model.get_params omits the b_0 deterministic (unlike
        # the Point_Process_Model base), so SVI never returns samples['b_0']
        # and plot_spatial(include_cov=True) raises KeyError. Register it
        # here. Note the deterministic is X@w only — the offset term is added
        # after registration, so plotted surfaces show the covariate field.
        if "spatial_cov" in self.args:
            pars["b_0"] = 0
        return pars
