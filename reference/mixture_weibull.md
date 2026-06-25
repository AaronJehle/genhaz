# Mixture Weibull survival, cumulative hazard, and hazard functions

A two-component mixture Weibull model with a general hazard (GH)
parametrisation. These functions are primarily intended for generating
synthetic survival data for simulation studies.

The survival function is \$\$S(t \mid X) = \bigl\[p \exp(-\lambda_1
t^{\gamma_1} e^{X \beta_1 \gamma_1}) + (1-p)\exp(-\lambda_2 t^{\gamma_2}
e^{X \beta_1 \gamma_2})\bigr\]^{\exp(\beta_2 X)}\$\$
