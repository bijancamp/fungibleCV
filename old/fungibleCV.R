##----------------------------------------------------------------------------------#
## fungibleCV - Simulation on the Cross-Validation Performance of Fungible Weights
## Bijan Camp
## April 6, 2014
##

#################################
######## OUTLINE OF CODE ########
#################################
##  ---------------------------------------------------------------------------------
##  I. Initial variable definitions:    Initialization of basic parameters and
##                                      storage structures for the results of
##                                      the simulation.
##
##  II. 'fungibleCV' function:          Main function that generates population
##                                      correlation matrices with fungibleR,
##                                      afterwards executing internal functions
##                                      defined in III to produce data for structures
##                                      defined in I.
##
##  III. Internal functions:
##
##      A. 'collect' function:          Simulates standard calibration/holdout data,
##                                      along with fungible weights associated with a
##                                      given R2 value. Passes its results to the
##                                      'analyze' function.
##
##      B. 'analyze' function:          Computes cross-validation metrics on collected
##                                      data. Saves the results in arrays defined in
##                                      part I.
##
##      C. 'summarize' function:        Tabulates the results.
##
##
##  Requirements:                       In addition to the libraries included below,
##                                      this simulation requires fungibleR.R for
##                                      producing correlation matrices (Waller, 2013,
##                                      under review) and fungibleWeights.R for
##                                      producing alternate weights that share a common
##                                      R2 value (Waller, 2008).
##
##  References:                         Waller, N.G. (2008). Fungible weights in
##                                           multiple regression. Psychometrica, 73,
##                                           691-703.
##                                      Waller, N.G. (2013, under review). Fungible
##                                           Correlation Matrices: A New Tool for
##                                           Evaluating Penalized Regression Models.
##



###########################
######## LIBRARIES ########
###########################
library(MASS) # provides mvrnorm and more
library(svMisc) # provides nice progress meter



##############################################
######## INITIAL VARIABLE DEFINITIONS ########
##############################################

## R.count:               number of fungible matrices to simulate data from
## predictor.count:       number of predictors for the data samples
## sample.count:          number of pairs of calibration/holdout samples to simulate
##                        from each fungible matrix
## sample.size:           number of data points to simulate in each calibration and
##                        holdout sample
## parameter.names:       descriptors of parameter sets for summarizing results
## parameter.count:       number of parameter sets to run simulation over
## summary.names:         descriptors of statistics computed on the data for
##                        summarizing results
## summary.count:         number of quantities to report in summary
## fungible.weight.count: number of fungible weights to cross-validate in each sample
R.count <- 100
predictor.count <- 4
sample.count <- 1000
sample.size <- 100
parameter.names <- c("Normal 1", "Normal 2", "Normal 3")
parameter.count <- length(parameter.names)
summary.names <- c("R2.pop", "R2.OLS.calibrated", "R2.OLS.holdout", "sd.R2.OLS.calibrated", "sd.R2.OLS.holdout", "R2.alt.calibrated", "R2.alt.holdout", "sd.R2.alt.holdout", "range.R2.alt.holdout", "RMS.R2.alt.calibrated.holdout")
summary.count <- length(summary.names)
fungible.weight.count <- 20

## Parameter space (array) consisting of sets of beta vectors and seed
## matrices. The first dimension of the array represents the index of
## a parameter pair, i.e., parameters[1, , ] represents the first
## parameter pair. The second and third dimensions together make up a
## matrix where the first row of the matrix is a population beta
## vector, and the remaining rows comprise a seed correlation matrix
## for use with fungibleR.
parameters <- array(0, c(parameter.count, predictor.count + 1, predictor.count))
parameters[1, , ] <- rbind(c(.1, .2, .3, .4), diag(4))
parameters[2, , ] <- rbind(c(.1, .2, .5, .6), diag(4))
parameters[3, , ] <- rbind(c(.1, .2, .4, .8), diag(4))

## Arrays containing the results of the simulation. Most of the arrays
## are three-dimensional. The first dimension of these
## (three-dimensional) arrays indexes the parameter sets, where each set
## consists of a beta vector and seed matrix. The second dimension of the
## arrays indexes the sampled data sets simulated from each fungible matrix.
## The third dimension indexes the fungible matrices.
R2.pop <- array(0, c(parameter.count))
R2.OLS.calibrated <- array(0, c(parameter.count, sample.count, R.count))
R2.OLS.holdout <- array(0, c(parameter.count, sample.count, R.count))
R2.alt.calibrated <- array(0, c(parameter.count, sample.count, R.count))
R2.alt.holdout <- array(0, c(parameter.count, sample.count, R.count))
sd.R2.alt.holdout <- array(0, c(parameter.count, sample.count, R.count))
range.R2.alt.holdout <- array(0, c(parameter.count, sample.count, R.count))
RMS.R2.alt.calibrated.holdout <- array(0, c(parameter.count, sample.count, R.count))



##########################################
######## MAIN SIMULATION FUNCTION ########
##########################################

fungibleCV <- function() {
    set.seed(2.718281828)

    ## array that will store each population correlation matrix
    R.fungible <- array(0, c(R.count, predictor.count, predictor.count))

    ## Main loop that runs over all sets of parameters
    for (parameter.counter in 1:parameter.count) {
        ## Retrieve pre-defined population quantities
        Beta <- parameters[parameter.counter, 1, ]
        R.seed <- parameters[parameter.counter, -1, ]

        ## Generate (sample.count) data sets for each fungible matrix
        for (R.counter in 1:R.count) {
            ## population fungible correlation matrix
            R.fungible[R.counter, , ] <- fungibleR(R.seed, Beta, Lp = .00, eps=1e-16)$Rstar
            
            ## population coefficients of determination
            R2.pop[parameter.counter] <<- t(Beta) %*% R.fungible[R.counter, , ] %*% Beta

            ## progress meter
            progress((parameter.counter - 1)*R.count + R.counter, parameter.count*R.count)
            cat("Using parameter set ", parameter.counter, ", fungible matrix ", R.counter, "...\n", sep = "")

            ## Simulate and analyze data from each fungible matrix
            for (sample.counter in 1:sample.count) {
                ## Simulate calibration/holdout sample, and produce R2 values
                ## associated with the OLS/fungible weights in the calibration/
                ## holdout sample
                data <- collect(R.fungible[R.counter, , ], parameter.counter)

                ## Analyze and store the results
                analyze(data, parameter.counter, sample.counter, R.counter)

            }# end of 'sample.count' loop

        }# end of 'R.count' loop
        
    }# end of 'parameter.count' loop

    ## Tabulate results
    summarize()

}# end of 'fungibleCV' function



####################################
######## INTERNAL FUNCTIONS ########
####################################

collect <- function(R, parameter.counter) {
    ## vector of correlations between criterion and predictors
    rxy <- R %*% Beta

    ## matrix of all correlations among criterion and predictors
    SM <- cbind(R, rxy)
    SM <- rbind(SM, c(t(rxy), 1))

    ## population mean for the Gaussian data sets
    mu <- rep(0, ncol(SM))

    ## Generate Gaussian data set
    XY <- mvrnorm(2*sample.size, mu, SM)
    
    ## Split data set into calibration and holdout samples
    random.rows <- sample(nrow(XY), size = sample.size)
    XY.calibrated <- XY[random.rows, ]
    XY.holdout <- XY[-random.rows, ]

    ## basic statistics computed on calibrated data
    SM.hat.calibrated <- cor(XY.calibrated)
    R.hat.calibrated <- SM.hat.calibrated[1:predictor.count, 1:predictor.count]
    rxy.hat.calibrated <- SM.hat.calibrated[1:predictor.count, predictor.count + 1]
    b.OLS.calibrated <- solve(R.hat.calibrated) %*% rxy.hat.calibrated
    R2.OLS.calibrated <- t(b.OLS.calibrated) %*% R.hat.calibrated %*% b.OLS.calibrated

    ## basic statistics computed on holdout data
    SM.hat.holdout <- cor(XY.holdout)
    R.hat.holdout <- SM.hat.holdout[1:predictor.count, 1:predictor.count]
    rxy.hat.holdout <- SM.hat.holdout[1:predictor.count, predictor.count + 1]
    b.OLS.holdout <- solve(R.hat.holdout) %*% rxy.hat.holdout
    R2.OLS.holdout <- t(b.OLS.holdout) %*% R.hat.holdout %*% b.OLS.holdout

    ## amount of decrease in R2 of alternate weights
    theta <- .01

    ## correlation between linear combinations of the beta weights and linear combinations of any set of alternate weights with the same R2
    r.yhata.yhatb.calibrated <- sqrt(1 - theta/R2.OLS.holdout)

    ## Number of fungible weights to generate and extract a smaller subset from
    sets <- 50*fungible.weight.count

    ## Generate fungible regression weights
    a.alt.calibrated.all <- Fungible(R.hat.calibrated, rxy.hat.calibrated, r.yhata.yhatb.calibrated, sets, print = FALSE)$a

    ## Randomly pick a few sets from these alternate weights
    a.alt.calibrated <- a.alt.calibrated.all[sample(nrow(a.alt.calibrated.all), size = fungible.weight.count), ]

    ## Calculate R2 in the calibrated sample for the alternate weight vectors
    ## (only one alternate weight vector needs to be analyzed as they all have
    ## the same R2)
    R2.alt.calibrated <- (t(a.alt.calibrated[1, ]) %*% rxy.hat.calibrated)^2/(t(a.alt.calibrated[1, ]) %*% R.hat.calibrated %*% a.alt.calibrated[1, ])

    ## R2 in the holdout sample for the alternate weight vectors
    R2.alt.holdout <- apply(a.alt.calibrated, 1, function(a) (t(a) %*% rxy.hat.holdout)^2/(t(a) %*% R.hat.holdout %*% a))

    ## output
    list(R2.OLS.calibrated = R2.OLS.calibrated, R2.OLS.holdout = R2.OLS.holdout, R2.alt.calibrated = R2.alt.calibrated, R2.alt.holdout = R2.alt.holdout)
}# End of 'collect' function



analyze <- function(data, parameter.counter, sample.counter, R.counter) {
    ## R2 of alternate and beta weights in calibrated/holdout samples
    R2.OLS.calibrated <- data$R2.OLS.calibrated
    R2.OLS.holdout <- data$R2.OLS.holdout
    R2.alt.calibrated <- data$R2.alt.calibrated
    R2.alt.holdout <- data$R2.alt.holdout
    mean.R2.alt.calibrated <- mean(R2.alt.calibrated)
    mean.R2.alt.holdout <- mean(R2.alt.holdout)
        
    ## Range of R2 of alternate weights in calibrated/holdout samples
    extrema.R2.alt.calibrated <- range(R2.alt.calibrated)
    min.R2.alt.calibrated <- extrema.R2.alt.calibrated[1]
    max.R2.alt.calibrated <- extrema.R2.alt.calibrated[2]
    range.R2.alt.calibrated <- abs(min.R2.alt.calibrated - max.R2.alt.calibrated)

    extrema.R2.alt.holdout <- range(R2.alt.holdout)
    min.R2.alt.holdout <- extrema.R2.alt.holdout[1]
    max.R2.alt.holdout <- extrema.R2.alt.holdout[2]
    range.R2.alt.holdout <- abs(min.R2.alt.holdout - max.R2.alt.holdout)

    ## standard deviation of the alternate weights
    sd.R2.alt.holdout <- sd(R2.alt.holdout)

    ## Root-mean-square error of the alternate weights in the holdout sample
    ## with respect to the alternate weights in the calibrated sample
    RMS.R2.alt.calibrated.holdout <- sqrt(sum((R2.alt.calibrated - R2.alt.holdout)^2)/fungible.weight.count)

    ## output
    R2.OLS.calibrated[parameter.counter, sample.counter, R.counter] <<- R2.OLS.calibrated
    R2.OLS.holdout[parameter.counter, sample.counter, R.counter] <<- R2.OLS.holdout
    R2.alt.calibrated[parameter.counter, sample.counter, R.counter] <<- mean.R2.alt.calibrated
    R2.alt.holdout[parameter.counter, sample.counter, R.counter] <<- mean.R2.alt.holdout
    sd.R2.alt.holdout[parameter.counter, sample.counter, R.counter] <<- sd.R2.alt.holdout
    range.R2.alt.holdout[parameter.counter, sample.counter, R.counter] <<- range.R2.alt.holdout
    RMS.R2.alt.calibrated.holdout[parameter.counter, sample.counter, R.counter] <<- RMS.R2.alt.calibrated.holdout
}# End of 'analyze' function



## Tabulate results for each parameter set
summarize <- function() {
    cat("\n======================\n")
    cat("Summary of Results\n")
    cat("======================\n")

    ## table to store results
    summary.table <- matrix(0, nrow = summary.count, ncol = length(parameter.names))

    ## table entries
    for (j in 1:parameter.count) {
        summary.table[1, j] <- round(R2.pop[j], 2)
        summary.table[2, j] <- round(mean(R2.OLS.calibrated[j, , ]), 2)
        summary.table[3, j] <- round(mean(R2.OLS.holdout[j, , ]), 2)
        summary.table[4, j] <- round(sd(R2.OLS.calibrated[j, , ]), 2)
        summary.table[5, j] <- round(sd(R2.OLS.holdout[j, , ]), 2)
        summary.table[6, j] <- round(mean(R2.alt.calibrated[j, , ]), 2)
        summary.table[7, j] <- round(mean(R2.alt.holdout[j, , ]), 2)
        summary.table[8, j] <- round(mean(sd.R2.alt.holdout[j, , ]), 2)
        summary.table[9, j] <- round(mean(range.R2.alt.holdout[j, , ]), 2)
        summary.table[10, j] <- round(mean(RMS.R2.alt.calibrated.holdout[j, , ]), 2)
    }

    ## labels
    rownames(summary.table) <- summary.names
    colnames(summary.table) <- parameter.names

    ## output
    print(summary.table)

}# End of 'summarize' function
