# Map locations of the 1000 Genomes populations

Approximate coordinates used to draw maps. Populations sampled away from
their region of ancestral origin (CEU, GIH, ITU, STU, MXL) are placed at
that region, which is more informative on a map of allele frequencies;
this is recorded in `placed_at`. Populations with recent ancestry from
several continents (ASW, ACB, and the admixed American populations other
than MXL) are placed where they were sampled.

## Usage

``` r
pop_locations
```

## Format

A tibble with 26 rows and 4 columns:

- pop:

  Population code.

- latitude, longitude:

  Approximate location, in decimal degrees.

- placed_at:

  Whether the point marks the sampling location or the ancestral origin.

## Source

Built by `data-raw/pop_locations.R`.
