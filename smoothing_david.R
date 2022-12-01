

###########
# Method 1
###########

# Idea: assuming periodicity of the underlying process 
# and assuming the data are observed over the period [0,T], 
# add to the original dataset: 
# (i) (0,yi) for all yi observed at time xi = T
# (ii) (T,yi) for all yi observed at time xi = 0
# In other words, duplicate the observations at the endpoints 
# and exchange their time coordinates (0 <--> T)
# Then use regular smoothing splines with cross-validation 
# separately for each gene profile


library(splines)
load("xy_common_genes.RData")

## Arrange the availability data in a time x gene data
o <- order(x$GeneID, x$adj.pt)
ntimes <- nrow(x) / ngenes
xmat <- matrix(x$log2.expr[o], ntimes, ngenes)

## Extend the dataset by duplicating data at boundaries
## of observation interval and switching times 0 and T
ptx <- x$adj.pt[1:ntimes]
uptx <- unique(ptx)
start <- which(ptx == ptx[1])
end <- which(ptx == ptx[ntimes])
x.ext <- xmat[c(end, 1:ntimes, start), ]
pt.ext <- ptx[c(end, 1:ntimes, start)]

## Prepare spline smoothing
xs <- matrix(, length(uptx), ngenes)
wx <- matrix(, length(pt.ext), ngenes) # data weights
sparx <- numeric(ngenes)
intrvl <- diff(ptx)
tol <- min(intrvl[intrvl > 0]) / 2
thresx <- 0.25 # tolerance for zero in data
lbx <- 1.1 # lower bound for smoothing parameter

## Initial fit by CV
for (i in 1:ngenes) {
	wx[,i] <- ifelse(x.ext[,i] < thresx, 1, 2)
	fit <- smooth.spline(x = pt.ext, y = x.ext[,i], w = wx[,i], 
		tol = tol, keep.data = FALSE)
	xs[,i] <- fit$y
	sparx[i] <- fit$spar
	if (i %% 100 == 0) cat(i, " ")
}

## Refit for genes whose profiles are undersmoothed
idx <- which(sparx < lbx)
for (i in idx) {
	fit <- smooth.spline(x = pt.ext, y = x.ext[,i], w = wx[,i], 
		spar = lbx, tol = tol, keep.data = FALSE)
	xs[,i] <- fit$y
}

## Display data + smoothed time courses
set.seed(2022)
idx <- sample(1:ngenes, 4)
dev.new()
pdf("availability.pdf")
par(mfrow = c(2, 2))
# idx <- 1:ngenes
for (i in idx)
{
	plot(ptx, xmat[,i], xlab = "Pseudotime", ylab = "Availability", 
		main = levels(x$GeneID)[i], ylim = c(0, 2.25))
	lines(uptx, xs[,i], col = "red", lwd = 1.5)
	Sys.sleep(.75)
}
dev.off()


####


## Repeat the previous steps for the expression data

o <- order(y$GeneID, y$adj.pt)
ntimes <- nrow(y) / ngenes
ymat <- matrix(y$log2.expr[o], ntimes, ngenes)
pty <- y$adj.pt[1:ntimes]
upty <- unique(pty)
start <- which(pty == pty[1])
end <- which(pty == pty[ntimes])
y.ext <- ymat[c(end, 1:ntimes, start), ]
pt.ext <- pty[c(end, 1:ntimes, start)]
ys <- matrix(, length(upty), ngenes)
wy <- matrix(, length(pt.ext), ngenes)
spary <- numeric(ngenes)
intrvl <- diff(upty)
tol <- min(intrvl[intrvl > 0]) / 2
thresy <- .25
lby <- 0.9
for (i in 1:ngenes) {
	wy[,i] <- ifelse(y.ext[,i] < thresy, 1, 2)
	fit <- smooth.spline(x = pt.ext, y = y.ext[,i], w = wy[,i],
		tol = tol, keep.data = FALSE)
	ys[,i] <- fit$y
	spary[i] <- fit$spar
	if (i %% 100 == 0) cat(i, " ")
}
idx <- which(spary < lby)
for (i in 1:ngenes) {
	fit <- smooth.spline(x = pt.ext, y = y.ext[,i], w = wy[,i], 
		spar = lby, tol = tol, keep.data = FALSE)
	ys[,i] <- fit$y
}
set.seed(2022)
idx <- sample(1:ngenes, 4)
dev.new()
pdf("expression.pdf")
par(mfrow = c(2,2))
# idx <- 1:ngenes
for (i in idx)
{
plot(pty, ymat[,i], xlab = "Pseudotime", ylab = "Expression", 
	main = levels(y$GeneID)[i]) # , ylim = c(0, 2.25))
lines(upty, ys[,i], col = "red", lwd = 1.5)
Sys.sleep(.75)
}
dev.off()




