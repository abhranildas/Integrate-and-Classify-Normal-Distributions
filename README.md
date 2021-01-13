# Integrate and Classify Normal Distributions [![View Integrate and classify normal distributions on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/82410-integrate-and-classify-normal-distributions)
Matlab toolbox to
 * integrate normal (Gaussian) distributions in any dimensions with any parameters, in any domain
 * compute pdf/cdf/inverse cdf of any function of a normal vector
 * compute measures of classification performance among any number of multivariate normal classes, like classification error, error matrix, and discriminability d'.

## Author
Abhranil Das, Center for Perceptual Systems, The University of Texas at Austin.
Bugs/comments/questions/suggestions to abhranil.das@utexas.edu.

If you use this code, please cite the accompanying paper: [A method to integrate and classify normal distributions.](https://arxiv.org/abs/2012.14331)

## Installation
Within Matlab's Home tab, select Add-Ons > Get Add-Ons > Search for 'Integrate and Classify Normal Distributions' and install.

## Quick Start
After installation, begin with the Getting Started live script with interactive examples:
 * run the command:
open(strcat(fileparts(which('integrate_normal')),filesep,'doc',filesep,'GettingStarted.mlx'))
 * or, after installation, click 'Open Folder' (or at any time, go to Matlab Home tab > Add-Ons > Manage Add-Ons > this toolbox > Open Folder). Then open the toolbox folder/doc/GettingStarted.mlx.

## Documentation
For function help, type:

    doc integrate_normal
    doc classify_normals
    doc classify_normals_multi
    doc norm_fun_cdf
    doc norm_fun_pdf
    doc norm_fun_inv
