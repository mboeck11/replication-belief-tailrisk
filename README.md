# Replication files for *Belief Shocks and Implications of Expectations about Growth-at-Risk*.

This page contains the replication files for Boeck, M. and Pfarrhofer, M. (2025) *Belief Shocks and Implications of Expectations about Growth-at-Risk*, Journal of Applied Econometrics, forthcoming. Codes are written in R.

Based on the codes provided, all figures of the article can be replicated. Please run the files in the given order in `/scripts`, as each file will create output stored in `/results`. Figures will be stored in `/figures`. All functions needed for estimation are stored in `/scripts/functions`. Set your working directory to `/`. For Figure 3, you have to run the script twice using different samples: the parameter `sample` has to be set once to `original` and once to `extended`. The time-varying parameter quantile regression (TVP-QR) is based on 2.000 draws, where we discard the first 1.000 draws as burnins. The VAR estimation is based on 15.000 MCMC draws, where we discard the first 5.000 as burnins. Please do note that the estimation can take some time (2-3h). Results are obtained using a Mac computer (mac OS Sequoia 15.1) and R 4.4.2.

- Figures:
  + Figure 1. Real GDP growth, estimated quantiles, and SPF nowcast distribution. (gdp_quantiles_rt_rw.pdf)
  + Figure 2. NEs computed with actual realizations of real GDP growth and the indicated selected GaR probabilities. (nowcast_errors_rt_rw.pdf)
  + Figure 3. Impulse response functions to a non-belief and belief shock. (beliefshocks_linear_original.pdf and beliefshocks_linear_extended.pdf)
  + Figure 4. Impulse response functions to belief shocks using GaR-NEs. (beliefshocks_tails_p=4_vs_ekm_rt_rw.pdf)
  + Figure A1. Scatterplot of NEs computed with actual realizations of real GDP and the indicated selected GaR probabilities and density estimates. (nowcast_scatter_v2_rt_rw.pdf)

**Abstract** This paper revisits the question of how shocks to expectations of market participants can cause business cycle fluctuations. We use a vector autoregression to estimate dynamic causal effects of belief shocks which are extracted from nowcast errors about output growth. In a first step, we replicate and corroborate the findings of Enders, Kleemann, and Müller (2021). The second step computes nowcast errors about growth-at-risk at various quantiles. This involves both recovering the quantiles of the nowcast distribution of output growth from the Survey of Professional Forecasters; and, since the true quantiles of output growth are unobserved, estimating them with quantile regressions. We document a lack of distinct patterns in response to shocks arising from nowcasts misjudging macroeconomic risk. Although the differences are statistically insignificant, belief shocks about downside risk seem to produce somewhat sharper business cycle fluctuations.

**Links:** [(Latest Version Sept 2024)](https://mboeck11.github.io/papers/BP2024JAE.pdf)
