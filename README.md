# compute_error_rate


dPrime_gauss will be wrt custom boundary when provided. dPrime_aprx doesn't care about any boundaries. If a flipped custom boundary is provided, accuracy may be <50% and dPrime_gauss will be negative. This is not a bug.

In the case of observation input, if you specify p_a, this will be interpreted as the prior probability of the fitted Gaussian a, and will be reflected in the gaussian outputs requested, i.e. d' gauss and acc gauss etc. But for the purely data outputs, category probabilities are implied simply by the counts of observations from each distribution, so they will not reflect any p_a input.
