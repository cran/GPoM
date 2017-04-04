## ---- eval = TRUE--------------------------------------------------------
    # load data
    data("Ross76")

## ---- eval = TRUE, fig.align='center'------------------------------------
    # time vector
    tin <- Ross76[,1]
    # single time series
    data <- Ross76[,3]
    # plot
    plot(tin, data, type='l',
         xlab = 't', ylab = 'y(t)', main = 'Original time series')

## ---- eval = TRUE, echo = TRUE-------------------------------------------
    # model dimension
    nVar = 3
    # maximum polynomial degree
    dMax = 2
    # generalized Polynomial modelling
    out1 <- gPoMo(data, tin=tin, dMax = dMax, nS = nVar, show = 0,
                  nPmin = 9, nPmax = 11, verbose = FALSE, IstepMax = 3000)

## ---- eval = TRUE--------------------------------------------------------
    which(out1$okMod == 1)

## ---- eval = TRUE, fig.align='center'------------------------------------
      # obtained model #3
      plot(out1$stockoutreg$model3[,1], out1$stockoutreg$model3[,2],
           type='l', col = 'red',
           xlab = 'y(t)', ylab = 'dy(t)/dt', main = 'Phase portrait')
      # original phase portrait
      lines(out1$filtdata[,1], out1$filtdata[,2])

## ---- eval = TRUE--------------------------------------------------------
      # obtained model #3
      visuEq(nVar, dMax, out1$models$model3)

## ---- eval = TRUE--------------------------------------------------------
    # multiple (three) time series
    data <- Ross76[,2:4]

## ---- eval = TRUE, fig.align='center'------------------------------------
      # x(t)
      plot(tin, data[,1], ylim = c(-6.5,12),
           type='l', col = 'blue',
           xlab = 't', ylab = 'x(t)', main = 'Original time series')
      # y(t)
      lines(tin, data[,2], type = 'l', col = 'brown')
      # z(t)
      lines(tin, data[,3], type = 'l', col = 'orange')

## ---- eval = TRUE--------------------------------------------------------
    # generalized Polynomial modelling
    out2 <- outGP <- gPoMo(data, tin = tin, dMax = 2, nS = c(1,1,1),  show = 0,
                           IstepMin = 10, IstepMax = 3000, nPmin = 7, nPmax = 8)

## ---- eval = TRUE, fig.align='center'------------------------------------
      # obtained model #5
      plot(out2$stockoutreg$model5[,1], out2$stockoutreg$model5[,2],
           type='l', col = 'red',
           xlab = 'y(t)', ylab = 'dy(t)/dt', main = 'Phase portrait')
      # original phase portrait
      lines(out2$filtdata[,1], out2$filtdata[,2])

## ---- eval TRUE----------------------------------------------------------
    visuEq(3, 2, outGP$models$model5, approx = 4)

## ---- eval = FALSE-------------------------------------------------------
#      # time vector
#      tin <- Ross76[,1]
#      # multiple (two) time series
#      data <- Ross76[,2:3]

## ---- eval = TRUE, fig.align='center'------------------------------------
      # obtained model #5
      plot(tin, data[,1], ylim = c(-6.5,6.5),
           type='l', col = 'blue',
           xlab = 't', ylab = 'x(t)', main = 'Original time series')
      # original phase portrait
      lines(tin, data[,2], type = 'l', col = 'brown')

## ---- eval = TRUE, fig.align='center'------------------------------------
      # correlation between Rössler-x and Rössler-y
      cor(data[,1], data[,2])

## ---- eval = TRUE, fig.align='center'------------------------------------
      plot(data[,1], data[,2],
           type='p', cex = 0.01, col = 'blue',
           xlab = 'x', ylab = 'y', main = 'Scatter plot (x,y)')

## ---- eval = TRUE--------------------------------------------------------
    # model template:
    EqS <- matrix(1, ncol = 3, nrow = 10)
    EqS[,1] <- c(0,0,0,1,0,0,0,0,0,0)
    EqS[,2] <- c(1,1,0,1,0,1,1,1,1,1)
    EqS[,3] <- c(0,1,0,0,0,0,1,1,0,0)
    visuEq(3, 2, EqS, substit = c('X','Y','Z'))

## ---- eval = TRUE--------------------------------------------------------
    # generalized Polynomial modelling
    out3 <- gPoMo(data, tin=tin, dMax = 2, nS=c(2,1), EqS = EqS,
                  show = 0, verbose = FALSE,
                  IstepMin = 10, IstepMax = 2000, nPmi = 9, nPmax = 11)

## ---- eval = TRUE, fig.align='center'------------------------------------
      # obtained model #2
      plot(out3$stockoutreg$model2[,1], out3$stockoutreg$model2[,2],
           type='l', col = 'red',
           xlab = 'y(t)', ylab = 'dy(t)/dt', main = 'Phase portrait')
      # original phase portrait
      lines(out3$filtdata[,1], out3$filtdata[,2])

## ---- eval = TRUE--------------------------------------------------------
    visuEq(3, 2, out3$models$model2, approx = 4)

## ---- eval = TRUE--------------------------------------------------------
    # load data
    data(sprottK)
    data(rossler)
    # multiple (six) time series
    data <- cbind(rossler,sprottK)[1:400,]

## ---- eval = TRUE, fig.align='center'------------------------------------
      tin = (0:(dim(data)[1] - 1)) * 1/20
      # For the Rössler-1976 system
      # x(t)
      plot(tin, data[,1], ylim = c(-6.5,8.5),
           type='l', col = 'red',
           xlab = 't', ylab = 'x(t)', main = 'Original time series')
      # y(t)
      lines(tin, data[,2], type = 'l', col = 'brown')
      # z(t)
      lines(tin, data[,3], type = 'l', col = 'orange')
      # For the Sprott-K system
      # u(t)
      lines(tin, data[,4], type = 'l', col = 'green')
      # v(t)
      lines(tin, data[,5], type = 'l', col = 'darkgreen')
      # w(t)
      lines(tin, data[,6], type = 'l', col = 'lightgreen')

## ---- eval = TRUE--------------------------------------------------------
    # generalized Polynomial modelling
    out4 <- gPoMo(data, dtFixe = 1/20, dMax = 2, nS = c(1,1,1,1,1,1),
                  show = 0, method = 'rk4',
                  IstepMin = 2, IstepMax = 3,
                  nPmin = 13, nPmax = 13)

## ---- eval = TRUE, fig.show='hold'---------------------------------------
      KL <- out4$models$model347
      v0 <- as.numeric(head(data,1))
      outNumi <- numicano(nVar = 6, dMax = 2, Istep=5000, onestep=1/20, KL=KL, v0=v0)
      # obtained model #347
      plot(outNumi$reconstr[,2], outNumi$reconstr[,3],
           type='l', col = 'red',
           xlab = 'x(t)', ylab = 'y(t)', main = 'Phase portrait')
      # original phase portrait
      lines(out4$filtdata[,1], out4$filtdata[,2])
      # original phase portrait
      plot(outNumi$reconstr[,5], outNumi$reconstr[,6],
           type='l', col = 'green',
           xlab = 'u(t)', ylab = 'v(t)', main = 'Phase portrait')
      # original phase portrait
      lines(out4$filtdata[,4], out4$filtdata[,5])

## ---- eval = TRUE--------------------------------------------------------
    visuEq(6, 2, out4$models$model347, approx = 2, substit = 1)

## ---- eval = TRUE--------------------------------------------------------
    # load data
    data("Ross76")
    # time vector
    tin <- Ross76[1:500,1]
    # single time series
    series <- Ross76[1:500,3]
    # plot
    plot(tin, series, type = 'l', col = 'gray')

## ---- eval = TRUE--------------------------------------------------------
    # some noise is added
    series[1:100] <- series[1:100] + 0.1 * runif(1:100, min = -1, max = 1)
    series[301:320] <- series[301:320] + 0.5 * runif(1:20, min = -1, max = 1)
    plot(tin, series, type = 'l', col = 'black')
    #
    # weighting function
    W <- tin * 0 + 1
    W[1:100] <- 0  # the first fourty values will not be considered
    W[301:320] <- 0  # twenty other values will not be considered either
    lines(tin, W, type = 'l', col = 'brown')

## ---- eval = TRUE--------------------------------------------------------
    reg <- gloMoId(series, dt=1/100, weight = W, nVar=3, dMax=2, show=1)

## ---- eval = TRUE--------------------------------------------------------
    visuEq(3, 2, reg$K, approx = 4)

## ---- eval = TRUE--------------------------------------------------------
    # first weight which value not equal to zero:
    i1 = which(reg$finalWeight == 1)[1]
    v0 <-  reg$init[i1,1:3]
    reconstr <- numicano(nVar=3, dMax=2, Istep=5000, onestep=1/250, PolyTerms=reg$K,
                          v0=v0, method="rk4")

## ---- eval = TRUE--------------------------------------------------------
    plot(reconstr$reconstr[,2], reconstr$reconstr[,3], type='l', lwd = 3,
         main='phase portrait', xlab='time t', ylab = 'x(t)', col='orange')
    # original data:
    lines(reg$init[,1], reg$init[,2], type='l',
          main='phase portrait', xlab='x', ylab = 'dx/dt', col='black')
    # initial condition
    lines(v0[1], v0[2], type = 'p', col = 'red')

