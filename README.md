# Classify
A set of tools related to computing classification accuracy between two Gaussian distributions or two sets of observations. This provides some functionalities that are not provided by the [classify function in MATLAB's Statistics and Machine Learning Toolbox](https://www.mathworks.com/help/stats/classify.html).

Author: Abhranil Das, Wilson Geisler Lab, Center for Perceptual Systems, The University of Texas at Austin.

Things you can do with this toolbox:

* Given the parameters of two Gaussian distributions, and optionally their prior probabilities,
  * compute the optimal classification accuracy (up to 3D),
  * return the coefficients of the quadratic equation of the boundary, and a set of points on it,
  * compute the discriminability/sensitivity index d' between the distributions, both exactly (up to 3D) and approximately (any dimensions).
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

## Inputs
For Gaussian parameter inputs,
* `mu_a`, `mu_b`: column vectors of means (scalars in 1D)
* `v_a`, `v_b`: variance-covariance matrices (scalar variances in 1D)

 

* `p_a`: prior probability of Gaussian a. Unless specified, equal priors are assumed by default.

For observation inputs,
* `obs_a`, `obs_b`: matrices of observations. Rows are observations, columns are variables.
<br><br/>
* `type`: input type. Specify as `obs` if inputs are observations instead of Gaussian parameters. 
* `custom_bd_coeffs`: coefficients of a custom quadratic boundary.


In the case of observation input, if you specify p_a, this will be interpreted as the prior probability of the fitted Gaussian a, and will be reflected in the gaussian outputs requested, i.e. d' gauss and acc gauss etc. But for the purely data outputs, category probabilities are implied simply by the counts of observations from each distribution, so they will not reflect any p_a input.

* `bPlot`: specify 0 if you don't want a plot. 1 by default.

## Outputs
The output is a `struct` with fields:

### `d_gauss_aprx`
Approximate d' between the distributions. This takes the covariance matrix of each distribution to be equal to the average, and computes the Mahalanobis distance between the distributions using this average covariance matrix. This approximation is useful when the distributions have too little overlap to compute d' exactly, or in >3 dimensions.
<br><br/>
### `bd_coeffs_gauss_opt`
Coefficients `a2`, `a1` and `a0` of the optimal quadratic boundary, written in matrix form:

<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{x}'&space;\boldsymbol{a_2}&space;\boldsymbol{x}&space;&plus;&space;\boldsymbol{a_1}'&space;\boldsymbol{x}&space;&plus;&space;a_0&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{x}'&space;\boldsymbol{a_2}&space;\boldsymbol{x}&space;&plus;&space;\boldsymbol{a_1}'&space;\boldsymbol{x}&space;&plus;&space;a_0&space;=&space;0" title="\boldsymbol{x}' \boldsymbol{a_2} \boldsymbol{x} + \boldsymbol{a_1}' \boldsymbol{x} + a_0 = 0" /></a>

If d is the number of dimensions, `a2` is a d x d matrix (scalar for 1D), `a1` is a length-d column vector (scalar for 1D), and `a0` is a scalar.
<br><br/>
### `bd_pts_gauss_opt`
A set of points on the optimal boundary between the Gaussians.
<br><br/>
### `acc_gauss`,`acc_gauss_a`, `acc_gauss_b`
Classification accuracy for the Gaussians, and for the individual Gaussians a and b.
<br><br/>
### `d_gauss`
Discriminability/sensitivity index d' between the Gaussian distributions. This always assumes equal priors and the optimal boundary.
<br><br/>
### `bd_coeffs_obs_opt`
Returned only for observation inputs. Matrix-form coefficients `a2`, `a1` and `a0` of the quadratic boundary that optimally separates the observations. We start with the optimal boundary between Gaussians fitted to the data, then optimize the quadratic coefficients to maximize the classification accuracy of the data.
<br><br/>
### `bd_pts_gauss_opt`
Returned only for observation inputs. A set of points on the above optimized boundary between the observations.
<br><br/>
### `acc_obs`,`acc_obs_a`, `acc_obs_b`
Returned only for observation inputs. Classification accuracy for the observation samples, and for the individual samples a and b.
<br><br/>

 ![alt text](https://github.com/abhranildas/classify/blob/master/summary_image.jpg)
