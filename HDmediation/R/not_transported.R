not_transported <- function(data, A, W, Z, M, Y, cens,
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
    ee <- e(data, npsem, folds, learners_e)
    ee <- apply(ee, 2, function(x) pmax(pmin(x, 1 - 0.001), 0.001))
    bb <- b(data, npsem, family, folds, learners_b)
    hz <- h_z(data, npsem, folds, learners_hz)

    if (!is.null(cens)) {
        prob_obs <- pObs(data, npsem, folds, learners_cens)
    }

    thetas <- ipws <- eifs <- list()
    vvbar <- matrix(nrow = nrow(data), ncol = 3)
    colnames(vvbar) <- c("00", "10", "11")
    for (param in list(c(1, 1), c(1, 0), c(0, 0))) {
        aprime <- param[1]
        astar <- param[2]

        hm <- h_m(hz, gg, ee, aprime, astar)

        if (!is.null(cens)) {
            obs <- data[[npsem$cens]]
            ipcw_ap <- obs / prob_obs[, gl("P(delta=1|A={aprime},Z,M,W)")]
            ipcw_as <- obs / prob_obs[, gl("P(delta=1|A={astar},Z,M,W)")]
        } else {
            ipcw_as <- ipcw_ap <- 1
        }

        ipwy <- ((A == aprime) / gg[, gl("g({aprime}|w)")])*ipcw_ap
        
        if (partial_tmle) {
            fit <- glm(Y ~ 1, offset = qlogis(bb[, gl("b({aprime},Z,M,W)")]), family = "binomial",
                       subset = A == aprime, weights = ipwy * hm #/ mean(ipwy * hm)
                       )
            bb[, gl("b({aprime},Z,M,W)")] <- plogis(coef(fit) + qlogis(bb[, gl("b({aprime},Z,M,W)")]))
        }

        uu <- u(data, npsem, bb, hm, aprime, folds, learners_u)
        uubar <- ubar(data, npsem, uu, aprime, folds, learners_ubar)
        vv <- v(data, npsem, bb, hz, aprime, folds, learners_v)
        vvbar[, paste(param, collapse = "")] <- vbar(data, npsem, vv, astar, folds, learners_vbar)

        # EIF calculation
        eify_weight <- ipwy * hm / mean(ipwy * hm)

        eify <- eify_weight * (Y - bb[, gl("b({aprime},Z,M,W)")])
        # eify <- ipwy * hm * (Y - bb[, gl("b({aprime},Z,M,W)")])
        
        ipwz <- ((A == aprime) / gg[, gl("g({aprime}|w)")])*ipcw_ap
        eifz_weight <- ipwz / mean(ipwz)
        eifz <- eifz_weight * (uu[, 1] - uubar[, 1])
        # eifz <- ipwz  * (uu[, 1] - uubar[, 1])
        
        ipwm <- ((A == astar) / gg[, gl("g({astar}|w)")])*ipcw_as
        eifm_weight <- ipwm / mean(ipwm)
        eifm <- eifm_weight * (vv[, 1] - vvbar[, paste(param, collapse = "")])
        # eifm <- ipwm  * (vv[, 1] - vvbar[, paste(param, collapse = "")])

        eif <- rescale_y(eify + eifz + eifm + vvbar[, paste(param, collapse = "")], bounds$bounds)
        theta <- mean(eif)

        thetas <- c(thetas, list(theta))
        ipws <- c(ipws, list(mean(ipwy * hm / mean(ipwy * hm) * Y)))
        eifs <- c(eifs, list(eif))
    }

    names(eifs) <- c("11", "10", "00")
    names(thetas) <- c("11", "10", "00")
    names(ipws) <- c("11", "10", "00")

    ans <- data.frame(total = thetas$`11` - thetas$`00`,
                      indirect = thetas$`11` - thetas$`10`,
                      direct = thetas$`10` - thetas$`00`,
                      gcomp_total = mean(vvbar[, "11"] - vvbar[, "00"]),
                      gcomp_indirect = mean(vvbar[, "11"] - vvbar[, "10"]),
                      gcomp_direct = mean(vvbar[, "10"] - vvbar[, "00"]), 
                      ipw_total = ipws$`11` - ipws$`00`,
                      ipw_indirect = ipws$`11` - ipws$`10`, 
                      ipw_direct = ipws$`10` - ipws$`00`)

    ans$var_total <- var(eifs$`11` - eifs$`00`)
    ans$var_indirect <- var(eifs$`11` - eifs$`10`)
    ans$var_direct <- var(eifs$`10` - eifs$`00`)
    
    ci_total <- ans$total + c(-1, 1) * qnorm(0.975) * sqrt(ans$var_total / nrow(data))
    ci_indirect <- ans$indirect + c(-1, 1) * qnorm(0.975) * sqrt(ans$var_indirect / nrow(data))
    ci_direct <- ans$direct + c(-1, 1) * qnorm(0.975) * sqrt(ans$var_direct / nrow(data))

    ans$ci_total_low <- ci_total[1]
    ans$ci_total_high <- ci_total[2]
    ans$ci_indirect_low <- ci_indirect[1]
    ans$ci_indirect_high <- ci_indirect[2]
    ans$ci_direct_low <- ci_direct[1]
    ans$ci_direct_high <- ci_direct[2]
    
    # TMLE code
    if(tmle == TRUE)
    {
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
        for (target in lvls[2:length(lvls)]) {
            Gs[folded[[i]]$validation_set, target] <- gg[,2]
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
    else
    {
        val_list <- list(
            results = ans,
            eif11 = eifs$`11`, 
            eif10 = eifs$`10`, 
            eif00 = eifs$`00`
            ) 
    }
}
