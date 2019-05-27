# Classify
A set of tools related to computing classification accuracy between two Gaussian distributions or two sets of observations. This provides some functionalities that are not provided by the [classify function in MATLAB's Statistics and Machine Learning Toolbox](https://www.mathworks.com/help/stats/classify.html). Things you can do with this toolbox:

* Given the parameters of two Gaussian distributions, and optionally their prior probabilities,
  * compute the optimal classification accuracy,
  * return the coefficients of the quadratic equation of the boundary, and a set of points on it,
  * compute the discriminability/sensitivity index d' between the distributions, both exactly and approximately.
  * produce a plot

* Input two sets of observations instead of Gaussian parameters. Gaussian parameters are then estimated and the above results are computed. In addition, a second quadratic boundary is returned that is optimized to better separate the observations, than just the boundary between the fitted Gaussians. This is more useful when the data are less Gaussian.

* Supply custom quadratic coefficients of the boundary.

## Quick start
Download the repository, and check out `demo.m` for a set of examples.

## Syntax
* Classify between two Gaussians with specified means `mu_a` and `mu_b`, and covariance matrices `v_a` and `v_b`:

  `results=classify([mu_a,v_a],[mu_b,v_b])`  

* Specify prior probability of distribution a (assumed equal by default):

  `results=classify([mu_a,v_a],[mu_b,v_b],'p_a',.7)`  

* Input observations instead of parameters:

  `results=classify(obs_a,obs_b,'type','obs')`
  
* Supply custom boundary coefficients:
  
  `results=classify(obs_a,obs_b,'type','obs','custom_bd_coeffs',custom_bd_coeffs)`

dPrime_gauss will be wrt custom boundary when provided. dPrime_aprx doesn't care about any boundaries. If a flipped custom boundary is provided, accuracy may be <50%.

In the case of observation input, if you specify p_a, this will be interpreted as the prior probability of the fitted Gaussian a, and will be reflected in the gaussian outputs requested, i.e. d' gauss and acc gauss etc. But for the purely data outputs, category probabilities are implied simply by the counts of observations from each distribution, so they will not reflect any p_a input.
