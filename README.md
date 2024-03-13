# Bayesian Fatigue Modeling
This repository contains code to fit a Bayesian model of stress vs. lifetime fatigue data. The data are assumed to come in pairs $(S_i,N_i)$ representing the stress and lifetime, respectively. The data generating process is modeled as 2-stage hierarchical model. For a given stress value $S$, we assume there is some probability of a the material reaching $10^7$ stress cycles without fracturing, defined as: 
$$\alpha(S) := p(\log(N)=7 | S) = \frac{1}{1 + e^{-(b_0 + b_1 S)}}$$ 
where above, $\alpha(S)$ is modeled as a logistic function of $S$.

Then, for the observed value $N$, with probability $\alpha(S)$,
$$\frac{-1}{N} = \beta_0 + \beta_1 S + \epsilon$$
(where $\epsilon$ is assumed to be a mean zero Gaussian random variable). With probability $1 - \alpha(S)$, $N=10^7$.

Using a Bayesian approach has several advantages over maximum likelihood in this case. Of particular note are:
 * uncertainty estiamtes based on the posterior distribution appear much more stable and reliable than, e.g., inverting the Fisher information matrix (which is likely a poor approximation to the parameters' behavior, even if it were easily computable)
 * since the dataset size is small, the regularization effect of the Bayesian model is quite helpful, especially considering the complexity of the model
 * there is no appeal to asymptotics for uncertainties of parameters estimates; we get the exact posterior distribution (up to Monte Carlo error)
 * all of the usual advantages, such as (1) rich information on the parameters in the form of probability distributions, (2) we are able to gently guide the model to better estimates with prior scientific knowledge, and (3) there is no need to come up with any sort of formula or estimate for each of the parameters, these are simply derived from the posterior samples
