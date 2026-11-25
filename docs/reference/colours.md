# Colours for populations and super-populations

Named vectors of colours used by all snpscaper plots. Use them in your
own plots, for example with
`ggplot2::scale_colour_manual(values = superpop_colours)`, or pass them
to factoextra:
`factoextra::fviz_pca_ind(pca, habillage = "super_pop", palette = superpop_colours)`.

## Usage

``` r
superpop_colours

pop_colours
```

## Format

Named character vectors of hex colour codes.

## Details

- `superpop_colours`: one colour for each of the five super-populations.

- `pop_colours`: one colour for each of the 26 populations, in shades of
  their super-population colour.

## Examples

``` r
superpop_colours
#>       AFR       AMR       EAS       EUR       SAS 
#> "#E69F00" "#D55E00" "#009E73" "#0072B2" "#CC79A7" 
head(pop_colours)
#>       ACB       ASW       ESN       GWD       LWK       MSL 
#> "#996900" "#B27B00" "#CC8D00" "#E69F00" "#EAB02E" "#EEC25D" 
```
