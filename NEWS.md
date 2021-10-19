# scutrboot 0.2.2

NEW FEATURES

* None

SIGNIFICANT USER-VISIBLE CHANGES

* `testTwoSample` respects the `minCellsPerGene` argument.

BUG FIXES

* Previously, the `minCellsPerGene` argument was effectively being ignored.

# scutrboot 0.2.1

NEW FEATURES

* Adjusts estimated p-value in two-sample test based on number of non-NA producing
  bootstraps.

SIGNIFICANT USER-VISIBLE CHANGES

* `testTwoSample` now returns number of valid bootstraps. Typically, these will 
  match the user-specified `nBootstraps`, but when bootstraps yield samples with 
  no expression in the tested gene, they become invalid. This occurs when a gene
  is expressed in too few cells. The p-values are adjusted to the number of 
  valid bootstraps. Users may want to exclude genes with low `bootstraps` values.
  Alternatively, the `minCellsPerGene` argument can be increased to help avoid 
  testing such genes.

BUG FIXES

* Removes `NA` bootstraps when estimating p-values for `testTwoSample`.
* P-value estimate includes *as extreme* bootstraps.

# scutrboot 0.2.0

NEW FEATURES

* Convert to new two-sample test.
* Generalizes usage index. User can specify index of choice (e.g., SUI, LUI, IPA)
* Adds Wasserstein Distance statistic.

SIGNIFICANT USER-VISIBLE CHANGES

* Two-sample test is now `testTwoSample`.
* Usage index tests use `statistic="UI"` with a user-provided `featureIndex`.
* Default test is multi-isoform via `statistic="WD"` for *Wasserstein Distance*,
  also known as the *Earth Mover's Distance*.

BUG FIXES

* None
