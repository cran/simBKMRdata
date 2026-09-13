# simBKMRdata 0.2.2

* Added a comprehensive `testthat` edition 3 test suite covering all exported
  functions, including PIP-threshold calculation, transformations, moment and
  Gamma-parameter estimation, multivariate Gamma generation, and grouped
  Gaussian/Gamma simulation workflows.
* Moved package tests to the standard `tests/testthat/` location used by
  `R CMD check`.
* Clarified the documented Gaussian simulation behavior: when
  `sampCorr_mat` is supplied to `MASS::mvrnorm()`, it is used directly as
  `Sigma`. Thus, a correlation matrix produces standardized Gaussian variables
  with unit marginal variances. The Gaussian simulation implementation itself
  is unchanged in this release.
* Added the package website/source repository and bug-report URLs to
  `DESCRIPTION`.
* Updated package references in `DESCRIPTION`.

# simBKMRdata 0.2.1

* Corrected estimates in `calculate_pip_threshold()`.
* Collected feedback from remaining authors and made corresponding changes in
  the vignettes; updated package author information.
* Prepared vignette materials for journal submission.

# simBKMRdata 0.1.1

* Fixed issues identified during CRAN review.

# simBKMRdata 0.1.0

* Initial release.
