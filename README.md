# Bayesian Fatigue Modeling
This repository contains code to fit a Bayesian model of stress vs. lifetime fatigue data. The data come in pairs $(S_i,N_i)$ representing the stress and lifetime, respectively. The data generating process is modeled as 2-stage hierarchical model. For a given stress value $S$, we assume there is some probability of a the material reaching $10^7$ stress cycles without fracturing, defined as: 
$$\alpha(S) := p(\log(N)=7 | S) = \frac{1}{1 + e^{-(b_0 + b_1 S)}}$$ 
where above, $\alpha(S)$ is modeled as a logistic function of $S$.

Then, for the observed value $N$, with probability $\alpha(S)$,
$$\frac{-1}{N} = \beta_0 + \beta_1 S + \epsilon$$
(where $\epsilon$ is assumed to be a mean zero Gaussian random variable). With probability $1 - \alpha(S)$, $N=10^7$.
