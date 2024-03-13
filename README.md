# Bayesian Fatigue Modeling
This repository contains code to fit a Bayesian model of stress vs. lifetime fatigue data. The data come in pairs $\{S,N\}$ representing the stress and lifetime. The data generating process is modeled as 2-stage hierarchical model. For a given stress value $S$, we assume there is some probability of a the material reaching 10^7 stress cycles without fracturing, defined as: 
$$\alpha(S) := p(log(N)=7 | S) = \exp{1/(b_0 + b_1 S)}$$ 
where above, $\alpha(S)$ is modeled as a logistic function of $S$.

For the observed value $N$, with probability $\alpha(S)$,
$$
\frac{-1}{N} = \beta_0 + \beta_1 S + \epsilon
$$
(where $\epsilon \sim N(0,\sigma^2)$), or $N=10^7$ otherwise.
