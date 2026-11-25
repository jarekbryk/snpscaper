# Colours are fixed so that every plot, handout and student figure uses the
# same colour for the same population. Super-population colours come from the
# Okabe-Ito palette, which stays distinguishable for people with the common
# forms of colour blindness. Each population gets a shade of its
# super-population colour, so related populations look related.

#' Colours for populations and super-populations
#'
#' Named vectors of colours used by all snpscaper plots. Use them in your own
#' plots, for example with `ggplot2::scale_colour_manual(values =
#' superpop_colours)`, or pass them to factoextra:
#' `factoextra::fviz_pca_ind(pca, habillage = "super_pop", palette =
#' superpop_colours)`.
#'
#' * `superpop_colours`: one colour for each of the five super-populations.
#' * `pop_colours`: one colour for each of the 26 populations, in shades of
#'   their super-population colour.
#'
#' @format Named character vectors of hex colour codes.
#' @examples
#' superpop_colours
#' head(pop_colours)
#' @name colours
NULL

#' @rdname colours
#' @export
superpop_colours <- c(
  AFR = "#E69F00",
  AMR = "#D55E00",
  EAS = "#009E73",
  EUR = "#0072B2",
  SAS = "#CC79A7"
)

# Populations grouped by super-population, in the order used in plots.
pops_by_superpop <- list(
  AFR = c("ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"),
  AMR = c("CLM", "MXL", "PEL", "PUR"),
  EAS = c("CDX", "CHB", "CHS", "JPT", "KHV"),
  EUR = c("CEU", "FIN", "GBR", "IBS", "TSI"),
  SAS = c("BEB", "GIH", "ITU", "PJL", "STU")
)

# Make n shades of a colour, from darker to lighter.
shades_of <- function(colour, n) {
  dark <- grDevices::colorRampPalette(c("black", colour))(10)[7]
  light <- grDevices::colorRampPalette(c(colour, "white"))(10)[6]
  grDevices::colorRampPalette(c(dark, colour, light))(n)
}

#' @rdname colours
#' @export
pop_colours <- unlist(
  purrr::imap(
    pops_by_superpop,
    function(pops, superpop) {
      stats::setNames(shades_of(superpop_colours[[superpop]], length(pops)), pops)
    }
  ),
  use.names = FALSE
) |>
  stats::setNames(unlist(pops_by_superpop, use.names = FALSE))

#' Colour and fill scales for populations
#'
#' ggplot2 scales that apply the fixed snpscaper colours ([superpop_colours],
#' [pop_colours]), so the same population always has the same colour.
#'
#' @param ... Passed on to [ggplot2::scale_colour_manual()] or
#'   [ggplot2::scale_fill_manual()], for example `name` or `breaks`.
#'
#' @return A ggplot2 scale.
#' @export
#' @examples
#' library(ggplot2)
#' df <- data.frame(super_pop = names(superpop_colours), n = c(661, 347, 504, 503, 489))
#' ggplot(df, aes(super_pop, n, fill = super_pop)) +
#'   geom_col() +
#'   scale_fill_superpop() +
#'   theme_jb()
scale_colour_superpop <- function(...) {
  ggplot2::scale_colour_manual(values = superpop_colours, ...)
}

#' @rdname scale_colour_superpop
#' @export
scale_fill_superpop <- function(...) {
  ggplot2::scale_fill_manual(values = superpop_colours, ...)
}

#' @rdname scale_colour_superpop
#' @export
scale_colour_pop <- function(...) {
  ggplot2::scale_colour_manual(values = pop_colours, ...)
}

#' @rdname scale_colour_superpop
#' @export
scale_fill_pop <- function(...) {
  ggplot2::scale_fill_manual(values = pop_colours, ...)
}
