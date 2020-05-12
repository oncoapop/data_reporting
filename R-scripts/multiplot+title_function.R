# Load R libraries
library(ggplot2)

# Outputs
Expt="NMD_Expr_vs_Highest_Isoform"
pdfname=paste(Expt,"pdf",sep=".")
graphtitle=gsub("_"," ",Expt)
makepdf = TRUE

#############################################
# Multiple plot function
# Modified by Damian Yap to add title row
# Jan 2017
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
# Added a title using one row on the bottom of all graphs
# title needs 'graphtitle' variable
# might need to change the heights ratio
graphspace = 29
titlespace = 1
##############################################

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    intlayout <- matrix(seq(from = 1, by = 1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
    layout <- rbind(intlayout,(numPlots+1))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout), heights = c(rep_len(graphspace, nrow(layout)-1), titlespace) ) ))
    grid.text(graphtitle, vp = viewport(layout.pos.row = ceiling(numPlots/cols)+1, layout.pos.col = 1:cols))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# multiplot(p1, p2, p3, cols=3)

