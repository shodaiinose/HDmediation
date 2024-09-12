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
    
    if (perc_above_099 + perc_below_001 < 0.05) {
        p <- 1e-02
    } else {
        p <- 1e-03
    }
    
    pmax(pmin(x, 1 - p), p)
}
