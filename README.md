# 2sCOPE vs. BMW: instrument-free control-function methods for endogeneity

Replication code and data for the article *Evaluating Two Alternative Instrument-Free Control-Function Methods for Addressing Endogeneity*.

The article compares two instrument-free control-function methods for linear regression with an endogenous regressor:

- **2sCOPE**, the two-stage copula endogeneity correction of Yang, Qian and Xie (2025), *Journal of Marketing Research*, 62(4), 601–623;
- **BMW**, the nonparametric control function of Breitung, Mayer and Wied (2024), *The Econometrics Journal*, 27(3), 362–383.

The comparison has two parts: Monte Carlo simulations under the assumptions of each method, and five real-data applications with instrumental-variable (IV) benchmarks.

## Repository structure

```
Monte-Carlo simulations/
    simulations_under_2sCOPE_hypothesis.do    data generated under the 2sCOPE assumptions
    simulations_under_BMW_hypothesis.do       data generated under the BMW assumptions
    corr_z_u_versus_rho.dta                   values of rho and corresponding corr(z,u) (Table 1)
    simulations_results_segun_copula_1_new_change_sample_size.dta   output, 2sCOPE assumptions
    simulations_results_segun_copula_2_new_change_sample_size.dta   output, BMW assumptions
Empirical_experiments/
    Original databases/      the five datasets as CSV files
    Cleaned databases/       the same datasets in Stata format, with harmonized variable names
    Estimations_with_bootstrap.do
    Results/Resultados_Estimaciones_bootstrap.xlsx
```

## Requirements

Stata 14 or later. Only built-in commands are used (`regress`, `bootstrap`, `lincom`, `sfrancia`, `export excel`).

## 1. Monte Carlo simulations (Section 4, Tables 2–9, Figures 1–9)

Run each do-file from the `Monte-Carlo simulations` folder (it reads `corr_z_u_versus_rho.dta` and saves its output there).

The model is y = βx + γz + u, with β = 0.7, γ = 0.4, x ~ N(0, 9), and seven degrees of endogeneity, corr(z,u) = 0, 0.1, …, 0.6 (Table 1). For each of the sample sizes N = 250, 500, 1,000, 2,000 and 4,000, there are 1,000 replications per degree of endogeneity. The seed (20261109) is reset for each sample size.

| Do-file | Data-generating process | Output |
|---|---|---|
| `simulations_under_2sCOPE_hypothesis.do` | Gaussian copula between z and u (Section 4.1) | `simulations_results_segun_copula_1_new_change_sample_size.dta` |
| `simulations_under_BMW_hypothesis.do` | u depends on a nonlinear transformation of the first-stage error (Section 4.2) | `simulations_results_segun_copula_2_new_change_sample_size.dta` |

Each output file has one row per replication. The variables are:

| Variable | Content |
|---|---|
| `rho`, `corr_z_u_teorico`, `corr_z_u` | endogeneity parameter, target and realized corr(z,u) |
| `beta_sin_correccion`, `gamma_sin_correccion` | OLS estimates |
| `beta_con_cupula1`, `gamma_con_cupula1`, `rho_con_cupula1` | 2sCOPE estimates |
| `beta_con_cupula2`, `gamma_con_cupula2`, `rho_con_cupula2` | BMW estimates |
| `N` | sample size |

The bias and the mean squared error reported in Tables 2–9 (N = 500) and in Figures 1–9 are the averages over replications of (estimate − true value) and (estimate − true value)², by method, degree of endogeneity and sample size.

## 2. Real-data applications (Section 5, Tables 10–16)

Run `Estimations_with_bootstrap.do` from the `Empirical_experiments` folder. It reads every `.dta` file in `Cleaned databases` and writes `Results/Resultados_Estimaciones_bootstrap.xlsx`, with one sheet per dataset.

All datasets follow the same convention: `y` is the dependent variable, `z` the endogenous regressor, `X1`, `X2`, … the exogenous regressors, `iv` or `iv1`, `iv2`, … the instruments (when available), and `cluster` the bootstrap cluster (when relevant).

| File | Application | Table | N | Instrument | Source |
|---|---|---|---|---|---|
| `MEPS.dta` | Health insurance and medical expenditure | 10 | 10,089 | SSI ratio | Medical Expenditure Panel Survey, as distributed by Katchova (2013), Econometrics Academy |
| `toothpaste_store111.dta`, `toothpaste_store89.dta` | Toothpaste demand | 11, 12 | 371 per store | price at the other store | Dominick's Finer Foods, Kilts Center for Marketing, University of Chicago Booth School of Business |
| `BLP.dta` | Automobile demand | 13 | 2,217 | the eight BLP instruments | Berry, Levinsohn and Pakes (1995), as distributed with PyBLP (Conlon and Gortmaker, 2020) |
| `diapers.dta` | Diaper demand | 14 | 261 | none | `diapers` dataset of the R package Rcope (Obrzut, Qian and Xie, 2026) |
| `airfare.dta` | Air travel demand | 15 | 4,596 | market concentration (`bmktshr`) | `AIRFARE` dataset of Wooldridge (2020) |

For each dataset, the do-file estimates OLS, 2sCOPE, BMW and IV (by control function) in the same bootstrap replications (1,000 replications, seed 20261110, resampling clusters when `cluster` exists). This gives standard errors for each estimator and for the differences between estimators. The results file has three specifications:

- `lineal`: the linear model (Tables 10–15);
- `con z2`: the model with z² added to the outcome equation, with the control function built for z only (Table 16);
- `cambio z2`: the change in the effect of z when z² is added, with a bootstrap standard error that accounts for the covariance between both specifications (Table 16).

## Data licenses

The datasets are redistributed here for replication purposes only. Their use is subject to the terms of the original providers listed above.

## Citation

[Reference of the article, to be completed after publication.]
