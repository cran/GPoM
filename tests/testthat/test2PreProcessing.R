#
# TestPreProcessing


# This file aims to show a simple example of time series
# preprocessing before using global modelling
#
# in particular we show the usefulness to resample the
# data - that will improve estimates of the derivatives -
# and the limitations reached, especially when considering
# derivatives of higher degree.

######################
# Example 1 ##########
######################
#devtools::wd(path='./tests/testthat/')

test_that("Single time series analysis", {

# original data set
#
# load
data("Ross76")
# plot
tin <- Ross76[,1]
data <- Ross76[,2:4]
dev.new()
plot(tin, data[,1], type='l')

# undersampled data
#
# subsampling:
Ross76us <- Ross76[seq(1,4000,by=50),]
# plot
tus <- Ross76us[,1]
datus <- Ross76us[,2:4]
lines(tus, datus[,1], type='p', col='red')
# resampling:
# (a spline is used here, however, alternative
# methods may be preferred depending on the dynamical
# behavior considered and on degree of subsampling,
# on the signal quality, etc.)
xout <- seq(min(tus), max(tus), by = 0.01)
rspl <- spline(tus, y = datus[,1], method = "fmm", xout = xout)
# plot resampled
lines(rspl$x, rspl$y, type='l', col='green')

# compute derivatives
#
# from original signal
drv1 <- drvSucc(tin, data[,1], nDeriv=2)
# from undersampled signal
drv2 <- drvSucc(tus, datus[,1], nDeriv=2)
# from resampled signal
drv3 <- drvSucc(rspl$x, rspl$y, nDeriv=2)

# plot resulting output as a function of time
dev.new()
par(mfrow = c(3, 1))
# output variable (smoothed) and its derivatives
plot(drv1$tout, drv1$seriesDeriv[,1], type='l', xlab = 'time', ylab = 'x(t)')
lines(drv2$tout, drv2$seriesDeriv[,1], type='p', cex = 1, col = 'red')
lines(drv3$tout, drv3$seriesDeriv[,1], type='l', col = 'green')
# first derivative
plot(drv1$tout, drv1$seriesDeriv[,2], type='l', xlab = 'time', ylab = 'dx(t)/dt')
lines(drv2$tout, drv2$seriesDeriv[,2], type='p', cex = 1, col = 'red')
lines(drv3$tout, drv3$seriesDeriv[,2], type='l', col = 'green')
# second derivative
plot(drv1$tout, drv1$seriesDeriv[,3], type='l', xlab = 'time', ylab = 'd²x(t)/dt²')
lines(drv2$tout, drv2$seriesDeriv[,3], type='p', cex = 1, col = 'red')
lines(drv3$tout, drv3$seriesDeriv[,3], type='l', col = 'green')

# plot resulting output as phase portraits
dev.new()
par(mfrow = c(2,2), pty = 's')
# output variable (smoothed) and its derivatives
 plot(drv1$seriesDeriv[,1], drv1$seriesDeriv[,2], type='l', xlab = 'x(t)', ylab = 'dx(t)/dt')
lines(drv2$seriesDeriv[,1], drv2$seriesDeriv[,2], type='l', col = 'red')
lines(drv3$seriesDeriv[,1], drv3$seriesDeriv[,2], type='l', col = 'green')
# first derivative
 plot(drv1$seriesDeriv[,1], drv1$seriesDeriv[,3], type='l', xlab = 'x(t)', ylab = 'd²x(t)/dt²')
lines(drv2$seriesDeriv[,1], drv2$seriesDeriv[,3], type='l', col = 'red')
lines(drv3$seriesDeriv[,1], drv3$seriesDeriv[,3], type='l', col = 'green')
# second derivative
 plot(drv1$seriesDeriv[,2], drv1$seriesDeriv[,3], type='l', xlab = 'dx(t)/dt', ylab = 'd²x(t)/dt²')
lines(drv2$seriesDeriv[,2], drv2$seriesDeriv[,3], type='l', col = 'red')
lines(drv3$seriesDeriv[,2], drv3$seriesDeriv[,3], type='l', col = 'green')

})


test_that("Multiple time series analysis", {

# original data set
#
# load
data("Ross76")
# plot
tin <- Ross76[,1]
data <- Ross76[,2:4]
dev.new()
par(mfrow = c(3, 1))
plot(tin, data[,1], type='l', xlab = 'time', ylab = 'x(t)')
plot(tin, data[,2], type='l', xlab = 'time', ylab = 'y(t)')
plot(tin, data[,3], type='l', xlab = 'time', ylab = 'z(t)')
dev.new()
par(mfrow = c(2, 2), pty = 's')
plot(data[,1], data[,2], type='l', xlab = 'x(t)', ylab = 'y(t)')
plot(data[,1], data[,3], type='l', xlab = 'y(t)', ylab = 'z(t)')
plot(data[,2], data[,3], type='l', xlab = 'x(t)', ylab = 'z(t)')

# undersampled data
#
# subsampling:
Ross76us <- Ross76[seq(1,4000,by=100),]
# plot
tus <- Ross76us[,1]
datus <- Ross76us[,2:4]
lines(tus, datus[,1], type='p', col='red')
# resampling:
# (a spline is used here, however, alternative
# methods may be preferred depending on the dynamical
# behavior considered and on degree of subsampling,
# on the signal quality, etc.)
xout <- seq(min(tus), max(tus), by = 0.01)
rsplx <- spline(tus, y = datus[,1], method = "fmm", xout = xout)
rsply <- spline(tus, y = datus[,2], method = "fmm", xout = xout)
rsplz <- spline(tus, y = datus[,3], method = "fmm", xout = xout)
# plot resampled
dev.new()
par(mfrow = c(2, 2), pty = 's')
plot(data[,1], data[,2], type='l', xlab = 'x(t)', ylab = 'y(t)')
lines(rsplx$y, rsply$y, type='l', col='green')
plot(data[,1], data[,3], type='l', xlab = 'y(t)', ylab = 'z(t)')
lines(rsplx$y, rsplz$y, type='l', col='green')
plot(data[,2], data[,3], type='l', xlab = 'x(t)', ylab = 'z(t)')
lines(rsply$y, rsplz$y, type='l', col='green')

# derivatives
# from original signal
drv1x <- drvSucc(tin, data[,1], nDeriv=2)
drv1y <- drvSucc(tin, data[,2], nDeriv=2)
drv1z <- drvSucc(tin, data[,3], nDeriv=2)
# from resampled signal
drv3x <- drvSucc(rsplx$x, rsplx$y, nDeriv=2)
drv3y <- drvSucc(rsply$x, rsply$y, nDeriv=2)
drv3z <- drvSucc(rsplz$x, rsplz$y, nDeriv=2)
# plot phase portraits
# from original signal
# from resampled signal
drv3x <- drvSucc(rsplx$x, rsplx$y, nDeriv=2)
drv3y <- drvSucc(rsply$x, rsply$y, nDeriv=2)
drv3z <- drvSucc(rsplz$x, rsplz$y, nDeriv=2)

})
