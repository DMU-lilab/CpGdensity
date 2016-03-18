library(plotly)

PlotHist2D <- function(plot.data) {
  p <- plot_ly(
    x = plot.data$score, 
    y = plot.data$density, 
    colorbar = list(
      x = 0.999, 
      y = 0.4, 
      len = 0.6, 
      thickness = 10, 
      ticks = "outside", 
      title = "", 
      titleside = "right"
    ), 
    colorscale = "Rainbow", 
    name = "title",
    type = "histogram2d", 
    xbins = list(
      end = 1.0, 
      size = 0.01, 
      start = 0
    ), 
    ybins = list(
      end = 1.0, 
      size = 0.01, 
      start = 0
    ), 
    nbinsy = 100,
    zsmooth = "best"
  )
  return (layout(p,
  xaxis = list(title = ""),
  yaxis = list(title = "")))
}
PlotHist <- function(plot.data) {
  return(plot_ly(y = df.cg$density, type="histogram"))
}
LayoutPlot <- function(plot, xdomain, ydomain, title = '') {
  return(layout(plot, 
                title = tile, 
                xaxis = list(domain = c(0, 0.6)), 
                yaxis = list(domain = c(0, 0.6))))
}

# all cgs

# plot.hist2d.cg.all <- PlotHist2D(df.cg)
# plot.hist2d.cg.island <- PlotHist2D(df.cg.island)
# plot.hist2d.cg.medium <- PlotHist2D(df.cg.medium)
# plot.hist2d.cg.high <- PlotHist2D(df.cg.high)

# plot.hist.cg.all.density <- PlotHist(df.cg)
# plot.hist.cg.all.score <- PlotHist(df.cg)
# plot.hist.cg.island.density <- PlotHist(df.cg.island)
# plot.hist.cg.medium.density <- PlotHist(df.cg.medium)
# plot.hist.cg.high.density <- PlotHist(df.cg.high)

# layout.hist2d.cg.all <- LayoutPlot(plot.hist2d.cg.all, c(0, 0.6), c(0, 0.6))
# layout.hist.cg.all.density <- LayoutPlot(plot.hist.cg.all.density, c(0.7, 1), c(0.7, 1))

PlotCGData <- function(plot.data, slice.data = NULL) {
  # construct the plots
  
  plot.hist2d <- plot_ly(
    x = plot.data$score, 
    y = plot.data$density, 
    colorbar = list(
      x = 0.99, 
      y = 0.3, 
      len = 0.7, 
      thickness = 20, 
      ticks = "outside", 
      title = "", 
      titleside = "right"
    ), 
    colorscale = "Rainbow", 
    name = "Density Map", 
    type = "histogram2d", 
    nbinsx = 100, 
    nbinsy = 100, 
    zsmooth = "best"
  )
  plot.hist.density <- plot_ly(y = plot.data$density, type="histogram", name = "CpG Density", histnorm = "probability density")
  plot.hist.score <- plot_ly(x = plot.data$score, type = "histogram", name = "Methylation Level", histnorm = "probability density")
  plot.line.pdf <- plot_ly(x = slice.data$dist, y = slice.data$pdf, 
                           type = "scatter", mode = "lines+markers", name = "Distribution")
  plot.line.cdf <- plot_ly(x = slice.data$dist, y = slice.data$cdf, 
                           type = "scatter", mode = "lines+markers", name = "Accumulation")
  
  # layout the plots
  
  layout.hist2d <- layout(plot.hist2d, 
                          title="", 
                          xaxis = list(domain = c(0, 0.7),
                                       title = 'Methylation Level'), 
                          yaxis = list(domain = c(0, 0.7),
                                       title = 'CpG Density'))
  layout.hist.density <- layout(plot.hist.density, 
                                xaxis = list(domain = c(0.75, 1),
                                             title = ''), 
                                yaxis = list(domain = c(0, 0.7), 
                                             title = '',
                                             showticklabels = FALSE))
  layout.hist.score <- layout(plot.hist.score,
                              xaxis = list(domain = c(0, 0.7),
                                           title = '',
                                           showticklabels = FALSE), 
                              yaxis = list(domain = c(0.75, 1),
                                           title = ''))
  if(!is.null(slice.data)) {
    layout.line.pdf <- layout(plot.line.pdf,
                              xaxis = list(domain = c(0.35, 0.65),
                                           title = ''), 
                              yaxis = list(domain = c(0.35, 0.65),
                                           title = '',
                                           showticklabels = FALSE))
    layout.line.cdf <- layout(plot.line.cdf,
                              xaxis = list(domain = c(0.75, 1),
                                           title = ''), 
                              yaxis = list(domain = c(0.75, 1),
                                           title = '',
                                           showticklabels = FALSE))
  }
  
  # do the plot
  if(!is.null(slice.data)) {
    subplot(layout.hist2d, layout.hist.density, layout.hist.score, layout.line.pdf, layout.line.cdf)
  } else
  {
    subplot(layout.hist2d, layout.hist.density, layout.hist.score)
  }
}

p.line <- add_trace(p, x = c(0, 1), y = c(1, 0), type = "scatter")
plot.line <- layout(p.line, title="abc", xaxis = list(domain = c(0, 0.6)), yaxis = list(domain = c(0, 0.6)))

plot.hist2d <- layout(p, title="abc", xaxis = list(domain = c(0, 0.6)), yaxis = list(domain = c(0, 0.6)))

plot.hist.density <- layout(plot_ly(y = df.cg$density, type="histogram"), xaxis = list(domain = c(0.7, 1)), yaxis = list(domain = c(0, 0.6)))
plot.hist.score <- layout(plot_ly(x = df.cg$score, type = "histogram"), xaxis = list(domain = c(0, 0.6)), yaxis = list(domain = c(0.7, 1)))

p.line.step <- plot_ly(x = v.band, y = v.step, type = "scatter", mode = "lines+markers")
plot.steps <- layout(p.line.step, xaxis = list(domain = c(0.7, 1)), yaxis = list(domain = c(0.7, 0.85)))

p.line.expand <- plot_ly(x = v.band, y = v.count, type = "scatter", mode = "lines+markers")
plot.expand <- layout(p.line.expand, xaxis = list(domain = c(0.7, 1)), yaxis = list(domain = c(0.85, 1)))

plot.all <- subplot(plot.hist2d, plot.hist.density, plot.hist.score, plot.steps, plot.expand)
plot.all

#pdf <- plotly_IMAGE(plot.all, format = "pdf", out_file="cgdensity.pdf")
