# Classify
A set of tools related to computing classification accuracy between two Gaussian distributions or two sets of observations. This provides some functionalities that are not provided by the [classify function in MATLAB's Statistics and Machine Learning Toolbox](https://www.mathworks.com/help/stats/classify.html).

## Credits
Abhranil Das, R Calen Walshe and Wilson Geisler, Center for Perceptual Systems, The University of Texas at Austin.

If you use this code, please cite: [A new method to compute classification error, Journal of Vision.](https://jov.arvojournals.org/article.aspx?articleid=2750251)

## Things you can do with this toolbox

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

  `results=bayes_classify([mu_a,v_a],[mu_b,v_b])`  

* Specify prior probability of distribution a (assumed equal by default):

  `results=bayes_classify([mu_a,v_a],[mu_b,v_b],'p_a',.7)`  

* Input observations instead of parameters:
  `results=bayes_classify(obs_a,obs_b,'type','obs')`
  
* Supply custom boundary coefficients:
  
  `results=bayes_classify(obs_a,obs_b,'type','obs','custom_bd_coeffs',custom_bd_coeffs)`

## Inputs
For Gaussian parameter inputs,
* `mu_a`, `mu_b`: column vectors of means (scalars in 1D)
* `v_a`, `v_b`: variance-covariance matrices (scalar variances in 1D)

For observation inputs,
* `obs_a`, `obs_b`: matrices of observations. Rows are observations, columns are variables.
<br><br/>
* `p_a`: prior probability of distribution a. Default is 0.5. For observation inputs, prior probabilities are by default the relative sample sizes. If `p_a` is specified in this case, it is taken to be the prior of the Gaussian fitted to the data, and affects only the Gaussian-fit-based outputs, not the observation-based outputs.
* `type`: input type. Specify as `obs` if inputs are observations instead of Gaussian parameters. 
* `custom_bd_coeffs`: coefficients of a custom quadratic boundary.
* `bPlot`: specify 0 if you don't want a plot. 1 by default.

## Outputs
The output is a `struct` with fields:

### `err_gauss`
Array of three numbers: classification error for the Gaussians, for Gaussian a, and for Gaussian b. For observation inputs, these correspond to Gaussians fitted to the data.
<br><br/>
### `log_err_gauss_max`
The log (base 10) of the Chernoff upper bound for classification error of the Gaussians, see ch. 2 of 'Pattern Classification' by Duda, Hart and Stork. This is expressed as the log to be able to represent very small errors. For observation inputs, this corresponds to Gaussians fitted to the data.
<br><br/>
### `d_gauss`
Discriminability/sensitivity index d' between the Gaussian distributions. This always assumes equal priors and the optimal boundary. For observation inputs, this corresponds to Gaussians fitted to the data.
<br><br/>
### `d_gauss_aprx`
Approximate d' between the distributions. This takes the covariance matrix of each distribution to be equal to the average, and computes the Mahalanobis distance between the distributions using this average covariance matrix. This approximation is useful when the distributions have too little overlap to compute d' exactly, or in >3 dimensions. For observation inputs, this corresponds to Gaussians fitted to the data.
<br><br/>
### `d_gauss_min`
The lower bound for d' corresponding to the equal-prior Chernoff lower bound for Gaussian classification accuracy. For observation inputs, this corresponds to Gaussians fitted to the data.
<br><br/>
### `bd_coeffs_gauss_opt`
Coefficients `a2`, `a1` and `a0` of the optimal quadratic boundary, written in matrix form:

<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{x}'&space;\boldsymbol{a_2}&space;\boldsymbol{x}&space;&plus;&space;\boldsymbol{a_1}'&space;\boldsymbol{x}&space;&plus;&space;a_0&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{x}'&space;\boldsymbol{a_2}&space;\boldsymbol{x}&space;&plus;&space;\boldsymbol{a_1}'&space;\boldsymbol{x}&space;&plus;&space;a_0&space;=&space;0" title="\boldsymbol{x}' \boldsymbol{a_2} \boldsymbol{x} + \boldsymbol{a_1}' \boldsymbol{x} + a_0 = 0" /></a>

If d is the number of dimensions, this is an equation of the length-d vector variable x. `a2` is a d x d matrix (scalar for 1D), `a1` is a length-d column vector (scalar for 1D), and `a0` is a scalar.

For observation inputs, this corresponds to Gaussians fitted to the data.
<br><br/>
### `bd_pts_gauss_opt`
A set of points on the optimal boundary between the Gaussians. For observation inputs, this corresponds to Gaussians fitted to the data.
<br><br/>
### `err_obs`
Returned only for observation inputs. Array of three numbers: classification error for the observation samples, for sample a, and for sample b.
<br><br/>
### `d_obs`
Returned only for observation inputs, and if the two samples are of equal size. d' between the observed distributions based on the classification accuracy of the data.
### `bd_coeffs_obs_opt`
Returned only for observation inputs. Matrix-form coefficients `a2`, `a1` and `a0` of the quadratic boundary that optimally separates the observations. We start with the optimal boundary between Gaussians fitted to the data, then optimize the quadratic coefficients to maximize the classification accuracy of the data.
<br><br/>
### `bd_pts_obs_opt`
Returned only for observation inputs. A set of points on the above optimized boundary between the observations.
<br><br/>

 ![Summary image](https://github.com/abhranildas/classify/blob/master/summary_image.png)
