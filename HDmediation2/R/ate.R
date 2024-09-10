ate <- function(data, A, W, Y, cens,
                            family, folds = 1, partial_tmle, bounds,
                            learners_g = "glm",
                            learners_b = "glm",
                            learners_cens = "glm") {
    npsem <- Npsem$new(A = A, W = W, Y = Y, cens = cens)
    folds <- make_folds(data, folds)
    
    bounds <- scale_y(data[[npsem$Y]], family, bounds)
    data[[npsem$Y]] <- Y <- bounds$y
    Y <- ifelse(is.na(Y), -999, Y)
    A <- data[[npsem$A]]
    
    gg <- g(data, npsem, folds, learners_g)
    qq <- mY(data, npsem, family, folds, learners_b)
    
    # probability censoring for ATE, excluding Z and M
    if (!is.null(cens)) {
        prob_obs_ate <- pObs_ate(data, npsem, folds, learners_cens)
    }
    
    # probability of observation, excluding Z and M
    if (!is.null(cens)) {
        obs <- data[[npsem$cens]]
        ipcw_a1_ate <- obs / prob_obs_ate[, "P(delta=1|A=1,W)"]
        ipcw_a0_ate <- obs / prob_obs_ate[, "P(delta=1|A=0,W)"]
    } else {
        ipcw_a1 <- ipcw_a0 <- 1
    }
    
    Y <- ifelse(is.na(data[[npsem$Y]]), -999, data[[npsem$Y]])
    mYa <- data[[npsem$A]]*qq[, "Q(1,W)"] + (1 - data[[npsem$A]])*qq[, "Q(0,W)"]
    
    # H_a for ATE, excluding Z and M
    H_a_ate <- ((data[[npsem$A]] / gg[, "g(1|w)"])*ipcw_a1_ate) - 
        (((1 - data[[npsem$A]]) / gg[, "g(0|w)"])*ipcw_a0_ate)
    
    eif_ate <- (Y - mYa)*H_a_ate + qq[, "Q(1,W)"] - qq[, "Q(0,W)"]
    
    ans <- data.frame(ate = mean(eif_ate))
    
    eif_ate <- eif_ate - ans$ate
    
    ans$var_ate <- var(eif_ate)

    ci_ate <- ans$ate + c(-1, 1) * qnorm(0.975) * sqrt(ans$var_ate / nrow(data))

    ans$ci_ate_low <- ci_ate[1]
    ans$ci_ate_high <- ci_ate[2]
    eif_mat <- matrix(eif_ate + ans$ate)
    colnames(eif_mat) <- c("eif_ate_uncentered")
    results <- list(ans, eif_mat)
    
    results
}
