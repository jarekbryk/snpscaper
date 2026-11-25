# Plot cross-entropy to help choose K

Shows the cross-entropy of every admixture run for each `K` (grey
points) and the best run for each `K` (line). Lower cross-entropy means
the model predicts left-out genotypes better. A sensible choice of `K`
is where the line stops falling steeply or reaches its lowest point. In
small genomic regions the curve is often flat or noisy, which itself
tells you that the region carries little information about population
structure.

## Usage

``` r
plot_cross_entropy(adm)
```

## Arguments

- adm:

  The result of
  [`run_admixture()`](https://jarekbryk.github.io/snpscaper/reference/run_admixture.md).

## Value

A ggplot object.

## Examples

``` r
# \donttest{
lct <- filter_snps(example_lct, min_maf = 0.05, prune_ld = TRUE)
adm <- run_admixture(lct, k = 2:5, reps = 3)
#> Running admixture for K = 2, 3, 4, and 5 (3 runs each, 12 SNPs)...
plot_cross_entropy(adm)

# }
```
