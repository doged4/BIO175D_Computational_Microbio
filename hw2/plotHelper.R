# Bio185d
# Joe Wirth, September 2022

library(docstring)
savePlotAsPng <- function(plot.in, out.file){
    #' savePlotAsPng
    #'
    #' Saves a plot as a png file
    #'
    #' plot.in:  the plot to be saved
    #'
    #' out.file: the filename for saving
    #'
    #' does not return


    # open file, plot, then close file
    png(out.file)
    plot(plot.in)
    dev.off()
}