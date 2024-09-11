crossfit <- function(train, valid, y, type = c("binomial", "continuous"), id = NULL, learners, bound = FALSE) {
    preds <- mlr3superlearner(data = train,
                              target = y,
                              library = learners,
                              outcome_type = match.arg(type),
                              folds = NULL,
                              newdata = valid,
                              group = id)$preds
    lapply(preds, function(x) bound(x))
}

bound <- function(x) {
    perc_above_099 <- mean(x > 0.99)
    perc_below_001 <- mean(x < 0.01)
    
    if (perc_above_099 < 0.05 && perc_below_001 < 0.05) {
        lower_bound <- 0.01
        upper_bound <- 0.99
    } else {
        lower_bound <- 0.001
        upper_bound <- 0.999
    }
    
    x <- pmax(pmin(x, upper_bound), lower_bound)
    
    x
}
