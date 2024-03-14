tmle_mlr3 <- function(data, trt, covar, outcome,
                      g_learners = "glm", Q_learners = "glm", c_learners = "glm",
                      outcome_type = c("binomial", "continuous"), 
                      super_efficient = FALSE,
                      .mlr3superlearner_folds = 10) {
    tmle_folds <- 1
    trt_obs <- data[[npsem$A]]
    folded <- origami::make_folds(data, V = tmle_folds)
    if (tmle_folds == 1) {
        folded[[1]]$training_set <- folded[[1]]$validation_set
    }
    
    obs <- !is.na(data[["Y"]])
    
    covar <- W
    trt <- "A"
    trt_factor <- dummy_cols(data, trt, remove_first_dummy = FALSE, remove_selected_columns = TRUE)
    lvls <- setdiff(names(trt_factor), names(data))
    
    for (i in 1:tmle_folds) {
        train <- as.data.table(data[folded[[i]]$training_set, ])
        train_trt_factor <- as.data.table(trt_factor[folded[[i]]$training_set, ])
        valid <- as.data.table(data[folded[[i]]$validation_set, ])
        valid_trt_factor <- as.data.table(trt_factor[folded[[i]]$validation_set, ])
        
        valids <- vector("list", length(lvls) + 1)
        names(valids) <- c("A", lvls)
        valids[["A"]] <- valid_trt_factor
        for (lvl in lvls) {
            other <- setdiff(lvls, lvl)
            
            valids[[lvl]] <- data.table::copy(valid_trt_factor)
            valids[[lvl]][[lvl]] <- 1
            valids[[lvl]][, (other) := lapply(.SD, function(x) rep(0, length(x))), .SDcols = other]
        }
        
        Qs <- b(data, npsem, family, folds, learners_b, tmle = TRUE)
        
        Qs <- list(A = ifelse(trt_obs == 0, Qs[, 1], Qs[, 2]),
                   A_0 = Qs[, 1],
                   A_1 = Qs[, 2])
        
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
        for (target in lvls[2:length(lvls)]) {
            Gs[folded[[i]]$validation_set, target] <- gg[, 2]
            Hs[folded[[i]]$validation_set, target] <-
                valid_trt_factor[[target]] / g2[folded[[i]]$validation_set, target]
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
    eics <- lapply(c("A", lvls), function(lvl) Hs[, lvl] * Cs[, 1] * (y - Qeps[, lvl]) + Qeps[, lvl] - psis[lvl])
    names(eics) <- c("A", lvls)
        
        eics <- lapply(eics, rescale_y_continuous)
        psis <- sapply(psis, rescale_y_continuous)
    }
    
    ses <- lapply(c("A", lvls), function(lvl) sqrt(var(eics[[lvl]]) / nrow(data)))
    names(ses) <- c("A", lvls)
    
    list(psi = psis,
         std.error = ses,
         ic = eics,
         g = g,
         prob_observed = prob_observed)
}
