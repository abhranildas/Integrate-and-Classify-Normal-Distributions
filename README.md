# Integrate and Classify Normal Distributions [![View Integrate and classify normal distributions on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/82410-integrate-and-classify-normal-distributions)
A set of tools related to computing classification accuracy between two Gaussian distributions or two sets of observations. This provides some functionalities that are not provided by the [classify function in MATLAB's Statistics and Machine Learning Toolbox](https://www.mathworks.com/help/stats/classify.html).


## Credits
Abhranil Das, R Calen Walshe and Wilson Geisler, Center for Perceptual Systems, The University of Texas at Austin.

If you use this code, please cite: [A new method to compute classification error, Journal of Vision.](https://jov.arvojournals.org/article.aspx?articleid=2750251)

## Things you can do with this toolbox

* Given the parameters of two Gaussian distributions, and optionally their prior probabilities,
  * compute the optimal classification accuracy (up to 3D),
  * return the coefficients of the quadratic equation of the boundary, and a set of points on it,
  * compute the discriminability/sensitivity index d' between the distributions, both exactly (up to 3D) and approximately (any dimensions).
  * produce a plot of the distributions and the decision boundary.

* Instead of Gaussian parameters, input sets of observations from the two distributions, whether or not Gaussian. Gaussian parameters are then estimated and the above results are computed. In addition, a second quadratic boundary is returned that is optimized to better separate the observations, than just the boundary between the fitted Gaussians. This is more useful when the data are less Gaussian.

* Perform classification involving outcome rewards and penalties.

* Supply a custom boundary.

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
* `custom_bd_coeffs`: coefficients `a2`, `a1` and `a0` of a custom quadratic boundary, written in matrix form: <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{x}'&space;\boldsymbol{a_2}&space;\boldsymbol{x}&space;&plus;&space;\boldsymbol{a_1}'&space;\boldsymbol{x}&space;&plus;&space;a_0&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{x}'&space;\boldsymbol{a_2}&space;\boldsymbol{x}&space;&plus;&space;\boldsymbol{a_1}'&space;\boldsymbol{x}&space;&plus;&space;a_0&space;=&space;0" title="\boldsymbol{x}' \boldsymbol{a_2} \boldsymbol{x} + \boldsymbol{a_1}' \boldsymbol{x} + a_0 = 0" /></a>. If d is the number of dimensions, this is an equation of the length-d vector variable x. `a2` is a d x d matrix (scalar for 1D), `a1` is a length-d column vector (scalar for 1D), and `a0` is a scalar.
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
Returned in the optimal case (equal priors, no custom outcome values or custom boundary). Discriminability/sensitivity index d' between the Gaussian distributions. For observation inputs, this corresponds to Gaussians fitted to the data.
<br><br/>
### `d_gauss_aprx`
Returned in the optimal case. Approximate d' between the distributions. This takes the covariance matrix of each distribution to be equal to the average, and computes the Mahalanobis distance between the distributions using this average covariance matrix. This approximation is useful when the distributions have too little overlap to compute d' exactly, or in >3 dimensions. For observation inputs, this corresponds to Gaussians fitted to the data.
<br><br/>
### `d_gauss_min`
Returned in the optimal case. The lower bound for d' corresponding to the Chernoff lower bound for Gaussian classification accuracy. For observation inputs, this corresponds to Gaussians fitted to the data.
<br><br/>
### `bd_coeffs`
Coefficients `a2`, `a1` and `a0` of the quadratic boundary between the Gaussians, written in matrix form:

<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{x}'&space;\boldsymbol{a_2}&space;\boldsymbol{x}&space;&plus;&space;\boldsymbol{a_1}'&space;\boldsymbol{x}&space;&plus;&space;a_0&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{x}'&space;\boldsymbol{a_2}&space;\boldsymbol{x}&space;&plus;&space;\boldsymbol{a_1}'&space;\boldsymbol{x}&space;&plus;&space;a_0&space;=&space;0" title="\boldsymbol{x}' \boldsymbol{a_2} \boldsymbol{x} + \boldsymbol{a_1}' \boldsymbol{x} + a_0 = 0" /></a>

If d is the number of dimensions, this is an equation of the length-d vector variable x. `a2` is a d x d matrix (scalar for 1D), `a1` is a length-d column vector (scalar for 1D), and `a0` is a scalar.

For observation inputs, this corresponds to the boundary between Gaussians fitted to the data. If custom boundary coefficients are supplied, these are the custom coefficients.
<br><br/>
### `bd_pts`
A set of points on the boundary between the distributions. For observation inputs, it is the boundary between Gaussians fitted to the data. If custom boundary coefficients are supplied, these are points on the custom boundary.
<br><br/>
### `outcome_counts_obs`
Returned when inputs are observations. 2x2 matrix of counts of the four classification outcomes. Entry (1,1): true category was *a* and was classified as *a*, entry (1,2): true category was *a* but was classified as *b*, etc.
<br><br/>
### `err_obs`
Returned only for observation inputs. Array of three numbers: classification error for the observation samples, for sample a, and for sample b.
<br><br/>
### `d_obs`
Returned when inputs are observations and in the optimal case (if the two samples are of equal size, and no custom boundary is specified). d' between the observed samples based on the optimal classification accuracy.
<br><br/>
### `bd_coeffs_obs`
Returned only for observation inputs. Matrix-form coefficients `a2`, `a1` and `a0` of the quadratic boundary that best separates the observations. We start with the calculated boundary between Gaussians fitted to the data, then optimize the quadratic coefficients to maximize the classification accuracy of the data (or the expected outcome value, if outcome values are supplied).
<br><br/>
### `bd_pts_obs`
Returned only for observation inputs. A set of points on the above boundary between the observations.
<br><br/>
### `outcome_vals_gauss`
Returned when outcome values are supplied in the input. A 2x2 matrix of the expected values of the four classification outcomes, given the masses of the two Gaussians in the two decision regions. Entry (1,1): true category was *a* and was classified as *a*, entry (1,2): true category was *a* but was classified as *b*, etc.
<br><br/>
### `ex_val_gauss`
Returned when outcome values are supplied in the input. Overall expected classification outcome value given the masses of the two Gaussians in the two decision regions.
<br><br/>
### `outcome_vals_obs`
Returned when inputs are observations and outcome values are supplied. A 2x2 matrix of total values from the four outcomes. Entry (1,1): true category was *a* and was classified as *a*, entry (1,2): true category was *a* but was classified as *b*, etc.
<br><br/>
### `ex_val_obs`
Returned when inputs are observations and outcome values are supplied. Overall expected classification outcome value based on the observations and the observation-based boundary.

 ![Summary image](https://github.com/abhranildas/classify/blob/master/summary_image.png)
