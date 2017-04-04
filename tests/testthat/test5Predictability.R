#
# TestPredictability
context("TestPredictability")
# This file aims to show how to estimate automatically
# the predictability of the models obtained with gPoMo

test_that("First example", {

# load data
data("Ross76")
# time vector
tin <- Ross76[seq(1, 3000, by = 8), 1]
# single time series
data <- Ross76[seq(1, 3000, by = 8), 3]
dev.new()
plot(tin, data)

# global modelling
# results are put in list outputGPoM
outputGPoM <- gPoMo(data[1:300], tin = tin[1:300], dMax = 2, nS=c(3),
                    show = 0, method = 'rk4',
                    nPmax = 12, IstepMin = 400, IstepMax = 401)
#
expect_equal_to_reference(outputGPoM, './Meta/TestPredictability_outGPoM.rds')
visuOutGP(outputGPoM)

#######################
# test predictability #
#######################
outpred <- predictab(outputGPoM, hp = 15, Nech = 30)
expect_equal_to_reference(outpred, './Meta/TestPredictability_outpred.rds')

})
