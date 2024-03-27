ate_est <- function(data, A, W, Z, M, Y, cens = NULL, S = NULL,
                      family = c("binomial", "continuous"), folds = 1,
                      partial_tmle = TRUE, bounds = NULL,
                      learners_g = c("glm"),
                      learners_e = c("glm"),
                      learners_c = c("glm"),
                      learners_b = c("glm"),
                      learners_hz = c("glm"),
                      learners_u = c("glm"),
                      learners_ubar = c("glm"),
                      learners_v = c("glm"),
                      learners_vbar = c("glm"),
                      learners_cens = "glm",
                      tmle = FALSE) {
    
    checkmate::assertDataFrame(data[, c(A, S, W, Z, M, Y)])
    checkmate::assertNumber(folds, lower = 1, upper = nrow(data) - 1)
    
    if (!is.null(S)) {
        ans <- transported(data, A, S, W, Z, M, Y, family, folds, partial_tmle,
                           bounds, learners_g, learners_e, learners_c, learners_b,
                           learners_hz, learners_u, learners_ubar, learners_v, learners_vbar)
        return(ans)
    }
    
    ate_calc(data, A, W, Z, M, Y, cens, family, folds, partial_tmle,
                    bounds, learners_g, learners_e, learners_b,
                    learners_hz, learners_u, learners_ubar, learners_v, learners_vbar, learners_cens,
                    tmle)
}

ate_calc <- function(data, A, W, Z, M, Y, cens,
                            family, folds = 1, partial_tmle, bounds,
                            learners_g = "glm",
                            learners_e = "glm",
                            learners_b = "glm",
                            learners_hz = "glm",
                            learners_u = "glm",
                            learners_ubar = "glm",
                            learners_v = "glm",
                            learners_vbar = "glm",
                            learners_cens = "glm",
                     tmle) {
    # new code
    require(fastDummies)
    require(data.table)
    outcome <- Y
    covar <- W
    trt <- A
    # end new code
    
    npsem <- Npsem$new(A = A, W = W, Z = Z, M = M, Y = Y, cens = cens)
    folds <- make_folds(data, folds)
    
    bounds <- scale_y(data[[npsem$Y]], family, bounds)
    data[[npsem$Y]] <- Y <- bounds$y
    Y <- ifelse(is.na(Y), -999, Y)
    A <- data[[npsem$A]]
    
    gg <- g(data, npsem, folds, learners_g)
    gg <- apply(gg, 2, function(x) pmax(pmin(x, 1 - 0.001), 0.001))
    Hs <- matrix(nrow = nrow(data), ncol = 3)
    if (!is.null(cens)) {
        prob_obs <- pObs(data, npsem, folds, learners_cens)
    }

        tmle_folds <- 1
        trt_obs <- data[[npsem$A]]
        folded <- origami::make_folds(data, V = tmle_folds)
        if (tmle_folds == 1) {
            folded[[1]]$training_set <- folded[[1]]$validation_set
        }
        
        obs <- !is.na(data[[outcome]])
        trt_factor <- dummy_cols(data, trt, remove_first_dummy = FALSE, remove_selected_columns = TRUE)
        lvls <- setdiff(names(trt_factor), names(data))
        
        for (i in 1:tmle_folds) {
            train <- as.data.table(data[folded[[i]]$training_set, ])
            #train_trt_factor <- as.data.table(trt_factor[folded[[i]]$training_set, ])
            #valid <- as.data.table(data[folded[[i]]$validation_set, ])
            valid_trt_factor <- as.data.table(trt_factor[folded[[i]]$validation_set, ])
            
            #valids <- vector("list", length(lvls) + 1)
            #names(valids) <- c("A", lvls)
            #valids[["A"]] <- valid_trt_factor
            #for (lvl in lvls) {
            #    other <- setdiff(lvls, lvl)
            #    
            #    valids[[lvl]] <- data.table::copy(valid_trt_factor)
            #    valids[[lvl]][[lvl]] <- 1
            #    valids[[lvl]][, (other) := lapply(.SD, function(x) rep(0, length(x))), .SDcols = other]
            #}
            
            Qs <- b(data, npsem, family, folds, learners_b, tmle = TRUE)
            
            Qs <- list(ifelse(trt_obs == 0, Qs[, 1], Qs[, 2]),
                       Qs[, 1],
                       Qs[, 2])
            
            names(Qs) <- c("A", lvls)
            
            Qs <- lapply(Qs, function(x) pmax(pmin(x, 1 - 0.0001), 0.0001))
            
            Cs <- matrix(nrow = nrow(data), ncol = 1, data = 1)
            prob_observed <- matrix(nrow = nrow(data), ncol = 1, data = 1)
            prob_obs2 <- ifelse(trt_obs == 0, prob_obs[,1], prob_obs[,2])
            if (!all(obs)) {
                prob_observed[folded[[i]]$validation_set, 1] <- prob_obs2
                
                Cs[folded[[i]]$validation_set, 1] <-
                    as.numeric(obs[folded[[i]]$validation_set]) / prob_observed[folded[[i]]$validation_set, 1]
            }
            
            Gs <- matrix(nrow = nrow(data), ncol = length(lvls))
            Hs <- matrix(nrow = nrow(data), ncol = length(lvls) + 1)
            colnames(Hs) <- c("A", lvls)
            colnames(Gs) <- c(lvls)
            gg2 <- ifelse(trt_obs == 0, gg[,1], gg[,2])
            for (target in lvls[2:length(lvls)]) {
                Gs[folded[[i]]$validation_set, target] <- gg2
                Hs[folded[[i]]$validation_set, target] <-
                    valid_trt_factor[[target]] / Gs[folded[[i]]$validation_set, target]
            }
            
            Gs[folded[[i]]$validation_set, lvls[1]] <- 1 - rowSums(Gs, na.rm = T)
            Hs[folded[[i]]$validation_set, lvls[1]] <-
                valid_trt_factor[[lvls[1]]] / Gs[folded[[i]]$validation_set, lvls[1]]
            
            Hs[, "A"] <- rowSums(Hs[, lvls] * as.matrix(trt_factor[, lvls]))
            # calculate tmle
            Qeps <- matrix(nrow = nrow(data), ncol = length(lvls) + 1)
            colnames(Qeps) <- c("A", lvls)
            
            tmle_data <- data.frame(y = train[[outcome]],
                                    Q_A = Qs[["A"]],
                                    H_A = Hs[, "A"]*Cs[, 1])
            
            fluc <- glm(y ~ -1 + offset(qlogis(Q_A)) + H_A, data = tmle_data[obs, ], family = binomial)
            eps <- coef(fluc)
            
            for (lvl in c("A", lvls)) {
                Qeps[folded[[i]]$validation_set, lvl] <- plogis(qlogis(Qs[[lvl]]) + eps*Hs[, lvl]*Cs[, 1])
            }
        }
        
        y <- ifelse(is.na(data[[outcome]]), -999, data[[outcome]])
        psis <- apply(Qeps, 2, mean)
        eics <- lapply(c("A", lvls), function(lvl) Hs[, lvl] * Cs[, 1] * (y - Qeps[, lvl]) + Qeps[, lvl])
        names(eics) <- c("A", lvls)
        
        ses <- lapply(c("A", lvls), function(lvl) sqrt(var(eics[[lvl]]) / nrow(data)))
        names(ses) <- c("A", lvls)
        
        tmle_res <- list(psi = psis,
                         std.error = ses,
                         ic = eics,
                         g = gg,
                         prob_observed = prob_observed)
        
        val_list <- list(
            results = ans,
            eif11 = eifs$`11`, 
            eif10 = eifs$`10`, 
            eif00 = eifs$`00`,
            tmle_res = tmle_res
        )
}
