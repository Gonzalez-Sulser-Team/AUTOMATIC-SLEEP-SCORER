
library(ggplot2)
library(reshape2)
library(magrittr)
library(viridisLite)
library(cowplot) # plot_grid(g1,g2,ncol=1)

mgg.plot <- function(...) {
  ggplot(...) +
    scale_x_continuous(expand=expansion(mult=c(.02,.02))) +
    scale_y_continuous(expand=expansion(mult=c(.02,.02))) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}

mgg.col <- function(i=1,k=256) {
  if (i==1) colorRampPalette(c("green","white","red"))(k)
  else if (i==2) colorRampPalette(c("darkblue","cyan","orange","red","darkred"),space="Lab",interpolate="spline")(k)
  else if (i==3) colorRampPalette(c("black","blue","green","brown","red","yellow","white"))(k)
  else if (i==4) magma(k)
}

mgg.grid <- function(g,xb=NULL,yb=NULL,n=5) {
  approx <- function(d,r=c(1,2,5,10)) {
    l <- log10(d)
    e <- floor(l)
    l <- l - e
    k <- which.min(abs(l-log10(r)))
    10^e*r[k]
  }
  s <- layer_scales(g)
  rx <- s$x$range$range
  if (is.null(xb)) xb <- approx((rx[2]-rx[1])/n)
  if (xb<Inf) xb <- xb*(ceiling(rx[1]/xb):floor(rx[2]/xb))
  if (xb==Inf || xb[1]>=xb[length(xb)]) xb <- NULL
  g$scales$scales[[1]]$breaks <- xb
  ry <- s$y$range$range
  if (is.null(yb)) yb <- approx((ry[2]-ry[1])/n)
  if (yb<Inf) yb <- yb*(ceiling(ry[1]/yb):floor(ry[2]/yb))
  if (yb==Inf || yb[1]>=yb[length(yb)]) yb <- NULL
  g$scales$scales[[2]]$breaks <- yb
  g
}

mgg.img <- function( g=NULL, m = if(class(g)[1]!="gg") g, xlim=c(1,ncol(m)), ylim=c(1,nrow(m)), marg=.5, q=c(.05,.95),
                     zlim=quantile(m,q,na.rm=T), col=mgg.col(4), gcol="#80808080", xb=NULL, yb=NULL ) {
  if (class(m)!="matrix") m <- t(as.matrix(m))
  if (length(marg)<2) marg <- c(marg,marg)
  dx <- marg[1]*(xlim[2]-xlim[1]+1)/ncol(m); xlim <- xlim+c(-dx,dx)
  dy <- marg[2]*(ylim[2]-ylim[1]+1)/nrow(m); ylim <- ylim+c(-dy,dy)
  m <- matrix(col[pmin(pmax((m[nrow(m):1,]-zlim[1])/(zlim[2]-zlim[1])*(length(col)-1)+1,1),length(col))],nrow(m),ncol(m))
  if (class(g)[1]!="gg") { ing <- T; g <- mgg.plot() } else ing <- F
  g <- g + expand_limits(x=xlim,y=ylim)
  if (ing | !is.null(xb) | !is.null(yb)) {
    g <- g + theme( panel.ontop = TRUE,
                    panel.background = element_rect(color = NA, fill = NA),
                    panel.grid = element_line(colour=gcol) )
    g <- mgg.grid(g,xb,yb)
  }
  return( g + annotation_raster(m,xlim[1],xlim[2],ylim[1],ylim[2]) )
}

mgg.data <- function(x, y=NULL, vname = if(cx) c("x","y") else "x", cx=T) {
  if (is.null(y)) { y <- x; x <- NULL }
  if (class(y)!="matrix"||is.data.frame(x)) y <- t(as.matrix(y))
  if (is.null(x)) { if (cx || nrow(y)==1) y <- rbind(1:ncol(y),y) }
  else y <- rbind(x,y,deparse.level=0)
  nn <- vname[1:min(length(vname),nrow(y))]
  k <- nrow(y) - length(nn)
  if (k>0) nn <- c(nn[-length(nn)],rownames(provideDimnames(cbind(1:(k+2)),base=vname[length(vname)]))[-1])
  n <- row.names(y)
  if (is.null(n)) n <- nn
  else n[n==""] <- nn[n==""]
  d <- data.frame(y[1,],val=y[2,],var=n[2])
  if (nrow(y)>2) for (i in 3:nrow(y)) d <- rbind(d,data.frame(y[1,],val=y[i,],var=n[i]))
  names(d)[1] <- n[1]
  d
}

mgg.lin <- function(g=NULL, x=NULL, y=NULL, xb=NULL, yb=NULL, vname=c("x","y")) {
  if (is.null(g) && is.null(x) && is.null(y)) stop("No data")
  if (!is.null(g) && class(g)[1]!="gg") {
    if (is.null(x) && is.null(y)) { y <- g } 
    else if (is.null(x)) x <- g
    else if (is.null(y)) { y <- x; x <- g }
    else stop("No plot")
    g <- NULL
  }
  d <- mgg.data(x=x,y=y,vname=vname)
  n <- names(d)
  if (is.null(g))
    return((mgg.plot(data=d) + geom_path(aes_string(x=n[1], y=n[2], color=n[3]))) %>% mgg.grid(xb,yb))
  else
    return(g + geom_path(data=d, aes_string(x=n[1], y=n[2], color=n[3])))
}

mgg.scat <- function(g=NULL, x=NULL, y=NULL, xb=NULL, yb=NULL, vname="x", shape=19) {
  if (is.null(g) && is.null(x) && is.null(y)) stop("No data")
  if (!is.null(g) && class(g)[1]!="gg") {
    if (is.null(x) && is.null(y)) { y <- g } 
    else if (is.null(x)) x <- g
    else if (is.null(y)) { y <- x; x <- g }
    else stop("No plot")
    g <- NULL
  }
  d <- mgg.data(x=x,y=y,vname=vname,cx=F)
  n <- names(d)
  if (is.null(g))
    return((mgg.plot(data=d) + geom_point(aes_string(x=n[1], y=n[2], color=n[3]), shape=shape)) %>% mgg.grid(xb,yb))
  else
    return(g + geom_point(data=d, aes_string(x=n[1], y=n[2], color=n[3]), shape=shape))
}

mgg.save <- function( filename, plot=last_plot(), path=NULL, scale=1,
                      font=NULL, dim=c(3.5,3.5), units=c("in", "cm", "mm"), dpi=300 ) {
  if (!is.null(font)) plot <- plot + theme(text=element_text(size=font))
  ggsave(filename, plot, path=path, scale=scale, width=dim[1], height=dim[2], units=units, dpi=dpi)
}
