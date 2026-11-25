# Summarise a snpscape object

Counts how many people from each population are in the data, and shows
each processing step with the number of people and SNPs remaining.

## Usage

``` r
# S3 method for class 'snpscape'
summary(object, ...)
```

## Arguments

- object:

  A `snpscape` object.

- ...:

  Not used.

## Value

A list with `people` (a tibble of counts per population) and `steps`
(the processing log), invisibly. Both are also printed.
