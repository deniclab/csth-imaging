# load dependencies
require(ggplot2)
require(proto)

GeomBoxplotJitterOutlier <- ggproto(GeomBoxplotJitterOutlier, ggplot2::GeomBoxplot, {
  draw <- function (., data, ..., outlier.colour = "black", outlier.shape = 16, 
                    outlier.size = 2, outlier.jitter=0) {
    defaults <- with(data, data.frame(x = x, xmin = xmin, xmax = xmax, 
                                      colour = colour, size = size, linetype = 1, group = 1, 
                                      alpha = 1, fill = alpha(fill, alpha), stringsAsFactors = FALSE))
    defaults2 <- defaults[c(1, 1), ]
    if (!is.null(data$outliers) && length(data$outliers[[1]] >= 
                                          1)) {
      pp<-position_jitter(width=outlier.jitter,height=0)
      p<-pp$adjust(data.frame(x=data$x[rep(1, length(data$outliers[[1]]))], y=data$outliers[[1]]),.scale)
      outliers_grob <- GeomPoint$draw(data.frame(x=p$x, y = p$y, colour = I(outlier.colour), 
                                                 shape = outlier.shape, alpha = 1, size = outlier.size, 
                                                 fill = NA), ...)
    }
    else {
      outliers_grob <- NULL
    }
    with(data, ggname(.$my_name(),grobTree(outliers_grob, GeomPath$draw(data.frame(y = c(upper, ymax),
                                                                                   defaults2), ...),
                                           GeomPath$draw(data.frame(y = c(lower, ymin), defaults2), ...),
                                           GeomRect$draw(data.frame(ymax = upper, ymin = lower, defaults), ...),
                                           GeomRect$draw(data.frame(ymax = middle, ymin = middle, defaults), ...))))
  }
  
  objname <- "boxplot_jitter_outlier"
  desc <- "Box and whiskers plot with jittered outlier"
  guide_geom <- function(.) "boxplot_jitter_outlier"
  
})
geom_boxplot_jitter_outlier <- GeomBoxplotJitterOutlier$build_accessor()